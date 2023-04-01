#import bpy
#import bmesh
#import mathutils
from typing import List, Optional, Tuple
from enum import Enum
TH = 0.00001
TH2 = TH**2
DE = "-3.2782554626464844e-07 -2.281630039215088 -0.13208013772964478:3.5762786865234375e-07 -1.5119284391403198 1.7138760089874268|-1.3113021850585938e-06 2.910860061645508 -5.257140159606934:3.0994415283203125e-06 5.783379077911377 1.631919264793396|-3.2782554626464844e-07 -2.281630039215088 -0.13208013772964478:-1.3113021850585938e-06 2.9108572006225586 -5.257140159606934:5.784740447998047 6.579355239868164 -2.426814079284668|3.5762786865234375e-07 -1.5119284391403198 1.7138760089874268:3.0994415283203125e-06 5.783379077911377 1.631919264793396:2.732125997543335 10.099564552307129 1.4211957454681396|5.784740447998047 6.579355239868164 -2.426813840866089:2.732125997543335 10.099563598632812 1.4211965799331665"

def same_coords(c1: mathutils.Vector, c2: mathutils.Vector) -> bool:
    return (c1-c2).length_squared < TH2

class StepRes(Enum):
    FINISHED = 1
    NOT_FINISHED = 2
    PART_FINISHED = 3

class BigPoint:
    def __init__(self, i: int):
        self.points: List[Point] = []
        self.i = i
        self.count = 0

    def add_point(self, point: "Point"):
        self.points.append(point)
        self.count += 1


class Point:
    def __init__(self, i: int, spline: "Spline", bpoint: "BigPoint"):
        self.bpoint: BigPoint = bpoint
        self.i: int = i
        self.prev_seg: Optional[Segment] = None
        self.post_seg: Optional[Segment] = None
        self.spline: Spline = spline


class Spline:
    def __init__(self, glist: "GlobalList"):
        self.points: List[Point] = []
        self.segments: List[Segment] = []
        self.glist = glist
        self.glist.add_spline(self)

    def add_point(self, coords: mathutils.Vector):
        i = 0
        added = False
        count = self.glist.get_count()
        bpoint = BigPoint(-1) # placeholder
        while not added:
            if i == count:
                bpoint = self.glist.create_bpoint(coords)
                added = True
            else:
                if same_coords(coords, self.glist.get_coords(i)):
                    bpoint = self.glist.get_bpoint(i)
                    added = True
            i += 1
        point = Point(i-1, self, bpoint)
        bpoint.add_point(point)
        self.points.append(point)
        p_num = len(self.points)
        if p_num > 1:
            seg = Segment(self.points[p_num - 2], point, self.glist)
            self.segments.append(seg)

    def round_spline(self):
        seg = Segment(self.points[-1], self.points[0], self.glist)
        self.segments.append(seg)


class Segment:
    def __init__(self, p1: Point, p2: Point, glist: "GlobalList"):
        self.quads: List[Quad] = []
        self.p1 = p1
        self.p2 = p2
        self.p1.post_seg = self
        self.p2.prev_seg = self
        self.finished = False
        glist.add_segment(self)

    def add_quad(self, quad: "Quad"):
        self.quads.append(quad)
        if len(self.quads) == 2:
            self.finished = True

class Quad:
    @staticmethod
    def step_left(segments: List[Segment], s: Segment):
        p1 = s.p1
        prev_seg = p1.prev_seg
        if prev_seg is None:
            return None
        if prev_seg not in segments:
            return None
        return prev_seg

    @staticmethod
    def step_right(segments: List[Segment], s: Segment):
        p2 = s.p2
        post_seg = p2.post_seg
        if post_seg is None:
            return None
        if post_seg not in segments:
            return None
        return post_seg

    @staticmethod
    def move_along_spline(segments: List[Segment], s0: Segment):
        edg: List[Segment] = [s0]
        s = s0
        bps: List[BigPoint] = []
        while s is not None:
            s = Quad.step_right(segments, s)
            if s == s0:
                raise ValueError("cycling!")
            if s is not None:
                edg.append(s)
                bps.append(s.p1.bpoint)
        s = s0
        while s is not None:
            s = Quad.step_left(segments, s)
            if s == s0:
                raise ValueError("cycling!")
            if s is not None:
                edg = [s] + edg
                bps = [s.p2.bpoint] + bps
        return edg, bps

    @staticmethod
    def verify_segment_connects_point(s: Segment, p: Point) -> bool:
        for point in p.bpoint.points:
            if point in (s.p1, s.p2):
                return True
        return False

    @staticmethod
    def get_neighbours(segments: List[Segment], edg: List[Segment]) -> List[Segment]:
        edges = (edg[0], edg[-1])
        ret: List[Optional[Segment]] = [None, None]
        right_added = False
        left_added = False
        points = [edges[0].p1, edges[1].p2]
        for i in range(2):
            e = edges[i]
            found = False
            n = 0
            while not found:
                s = segments[n]
                if s == e:
                    found = True
                    if not left_added:
                        prev_seg = segments[n-1]
                        if prev_seg not in edg:
                            if Quad.verify_segment_connects_point(prev_seg, points[i]):
                                ret[i] = prev_seg
                                left_added = True
                    if not right_added:
                        nn = n + 1
                        if nn == len(segments):
                            nn = 0
                        next_seg = segments[nn]
                        if segments[nn] not in edg:
                            if Quad.verify_segment_connects_point(next_seg, points[i]):
                                ret[i] = segments[nn]
                                right_added = True
                n += 1
        ret2 = [seg for seg in ret if seg is not None]
        assert len(ret2) == 2
        return ret2

    @staticmethod
    def steps_dividing_edges_right(p: Point, n: int):
        for _ in range(n):
            p_post = p.post_seg
            if p_post is None:
                return None
            p = p_post.p2
        return p.bpoint

    @staticmethod
    def steps_dividing_edges_left(p: Point, n: int):
        for _ in range(n):
            p_prev = p.prev_seg
            if p_prev is None:
                return None
            p = p_prev.p1
        return p.bpoint

    @staticmethod
    def verify_dividing_edges(bp1: BigPoint, bp2: BigPoint, n: int):
        for p in bp1.points:
            s1 = Quad.steps_dividing_edges_right(p, n)
            if s1 is not None:
                if s1 == bp2:
                    return True
            s2 = Quad.steps_dividing_edges_left(p, n)
            if s2 is not None:
                if s2 == bp2:
                    return True
        return False

    @staticmethod
    def get_dir_from_right(edg: List[Segment], segments: List[Segment]):
        bp = edg[-1].p2.bpoint
        for p in bp.points:
            post = p.post_seg
            if post is not None:
                if post in segments:
                    return True
        return False

    @staticmethod
    def get_dir_from_left(edg: List[Segment], segments:List[Segment]):
        bp = edg[0].p1.bpoint
        for p in bp.points:
            prev = p.prev_seg
            if prev is not None:
                if prev in segments:
                    return True
        return False

    @staticmethod
    def verify_div(b1: List[BigPoint], b2: List[BigPoint], dir1: bool, dir2: bool, n: int):
        for i in range(len(b1)):
            if dir1 != dir2:
                if Quad.verify_dividing_edges(b1[i], b2[i], n):
                    return True
            else:
                if Quad.verify_dividing_edges(b1[i], b2[-i-1], n):
                    return True
        return False

    @staticmethod
    def verify_and_init(segments: List[Segment], glist: "GlobalList") -> bool:
        le = len(segments)
        if le < 4:
            return False
        if (le % 2) == 1:
            return False
        for segment in segments:
            if segment.finished:
                return False
        edg1, b1 = Quad.move_along_spline(segments, segments[0])
        ret = Quad.get_neighbours(segments, edg1)
        edg2, b2 = Quad.move_along_spline(segments, ret[1])
        edg4, b4 = Quad.move_along_spline(segments, ret[0])
        if len(edg2) != len(edg4):
            return False
        if le - len(edg2)*2 != len(edg1)*2:
            return False
        ret = Quad.get_neighbours(segments, edg2)
        if ret[0] not in edg1:
            seg = ret[0]
        else:
            seg = ret[1]
        edg3, b3 = Quad.move_along_spline(segments, seg)
        if le != len(edg1) + len(edg2) + len(edg3) + len(edg4):
            return False
        dir1 = True
        dir2 = Quad.get_dir_from_right(edg1, segments)
        dir4 = Quad.get_dir_from_left(edg1, segments)
        if dir2:
            dir3 = Quad.get_dir_from_right(edg2, segments)
        else:
            dir3 = not Quad.get_dir_from_left(edg2, segments)
        if Quad.verify_div(b1, b3, dir1, dir3, len(edg2)) or Quad.verify_div(b2, b4, dir2, dir4, len(edg1)):
            return False
        Quad(edg1, edg2, edg3, edg4, dir1, dir2, dir3, dir4, glist)
        return True

    def add_quad_to_edge(self, edg: List[Segment]):
        for segment in edg:
            segment.add_quad(self)

    def __init__(self,
                 edge1: List[Segment],
                 edge2: List[Segment],
                 edge3: List[Segment],
                 edge4: List[Segment],
                 dir1: bool,
                 dir2: bool,
                 dir3: bool,
                 dir4: bool,
                 glist: "GlobalList"):
        self.edges = [edge1, edge2, edge3, edge4]
        self.directions = [dir1, dir2, dir3, dir4]
        for edge in self.edges:
            self.add_quad_to_edge(edge)
        glist.add_quad(self)

    '''def progress_segment(self,
                         nedges: int,
                         vert: bmesh.types.BMVert,
                         end_co: List[float]) -> Tuple[List[bmesh.types.BMEdge], bmesh.types.BMVert]:
        selected_edges: List[bmesh.types.BMEdge] = []
        for edge in vert.link_edges:
            selected_edges = [edge]
            prev_vert = vert
            now_edge: bmesh.types.BMEdge = edge
            next_vert = now_edge.other_vert(prev_vert)
            for _ in range(nedges - 1):
                link_edges : List[bmesh.types.BMEdge] = next_vert.link_edges
                assert len(link_edges) == 2
                edg0 = link_edges[0]
                next_edge = edg0 if edg0 != now_edge else link_edges[1] 
                selected_edges.append(next_edge)
                now_edge = next_edge
                next_vert = now_edge.other_vert(prev_vert)
            if same_coords(next_vert.co, end_co):
                return selected_edges, next_vert
        raise ValueError("segment on mesh not found")'''

    '''def select_quad_edge(self, edge: List[Segment], nedges: int, glist: "GlobalList", bm: bmesh.types.BMesh):
        selected_edges: List[bmesh.types.BMEdge] = []
        vert = glist.search_vert(bm, self.edges[0][0].get_p1())
        for segment in edge:
            res = self.progress_segment(nedges, vert, glist.get_coords(segment.get_p2().get_i()))
            selected_edges.extend(res[0])
            vert = res[1]
        return selected_edges'''

    '''def render_quad(self, bm: bmesh.types.BMesh, nedges: int, glist: "GlobalList"):
        edges = self.select_quad_edge(self.edges[0], nedges, glist, bm)
        edges.extend(self.select_quad_edge(self.edges[2], nedges, glist, bm))
        bmesh.ops.grid_fill(bm, edges=edges)'''

class GlobalList:
    def __init__(self):
        self.reduced_points: List[mathutils.Vector] = []
        self.big_points: List[BigPoint] = []
        self.count = 0
        self.splines: List[Spline] = []
        self.quads: List[Quad] = []
        self.segments: List[Segment] = []

    def create_bpoint(self, coords: mathutils.Vector):
        self.reduced_points.append(coords)
        bpoint = BigPoint(self.count)
        self.big_points.append(bpoint)
        self.count += 1
        return bpoint

    def get_count(self):
        return self.count

    def get_coords(self, i: int):
        return self.reduced_points[i]

    def get_bpoint(self, i: int):
        return self.big_points[i]

    def add_spline(self, spline: Spline):
        self.splines.append(spline)

    def get_splines(self):
        return self.splines

    def add_segment(self, segment: Segment):
        self.segments.append(segment)

    def add_quad(self, quad: Quad):
        self.quads.append(quad)

    def work_with_segment(self, segment: Segment):
        segments_verified = [segment]
        bpoints_verified: List[BigPoint] = []
        counter_splines = 1
        prev_spline = segment.p1.spline
        initial_spline = prev_spline
        initial_bpoint = segment.p1.bpoint
        self.edge_step(initial_bpoint, segment.p2, segments_verified, bpoints_verified, counter_splines, prev_spline, initial_spline, True)
        segment.finished = True

    def edge_step(self,
                  initial_bpoint: BigPoint,
                  point: Point,
                  segments_verified: List[Segment],
                  bpoints_verified: List[BigPoint],
                  counter_splines: int,
                  prev_spline: Spline,
                  initial_spline: Spline,
                  first: bool) -> StepRes:
        bpoint = point.bpoint
        spline = point.spline
        if spline != prev_spline:
            counter_splines += 1
            first = False
            if counter_splines > 5:
                return StepRes.NOT_FINISHED
            if counter_splines == 5 and spline != initial_spline:
                return StepRes.NOT_FINISHED
        if initial_bpoint == bpoint:
            if counter_splines < 4:
                return StepRes.NOT_FINISHED
            if Quad.verify_and_init(segments_verified, self):
                if segments_verified[0].finished:
                    return StepRes.FINISHED
                return StepRes.PART_FINISHED
        if bpoint in bpoints_verified:
            return StepRes.NOT_FINISHED
        bpoints_verified.append(bpoint)
        for point in bpoint.points:
            segment1 = point.prev_seg
            res = self.step_segment(initial_bpoint,
                                    segment1,
                                    True,
                                    segments_verified,
                                    bpoints_verified,
                                    counter_splines,
                                    spline,
                                    initial_spline,
                                    first)
            if res in (StepRes.FINISHED, StepRes.PART_FINISHED):
                bpoints_verified.pop()
                return res
            segment2 = point.post_seg
            res = self.step_segment(initial_bpoint,
                                    segment2,
                                    False,
                                    segments_verified,
                                    bpoints_verified,
                                    counter_splines,
                                    spline,
                                    initial_spline,
                                    first)
            if res in (StepRes.FINISHED, StepRes.PART_FINISHED):
                bpoints_verified.pop()
                return res
        bpoints_verified.pop()
        return StepRes.NOT_FINISHED

    def step_segment(self,
                     initial_bpoint: BigPoint,
                     segment: Optional[Segment],
                     is_p1: bool,
                     sv: List[Segment],
                     bv: List[BigPoint],
                     counter_splines: int,
                     spline: Spline,
                     initial_spline: Spline,
                     first: bool) -> StepRes:
        if segment is not None:
            if not segment in sv:
                if not segment.finished:
                    sv.append(segment)
                    point = segment.p1 if is_p1 else segment.p2
                    res = self.edge_step(initial_bpoint, point, sv, bv, counter_splines, spline, initial_spline, first)
                    sv.pop()
                    match res:
                        case StepRes.FINISHED: 
                            return StepRes.FINISHED
                        case StepRes.PART_FINISHED:
                            if not first:
                                return StepRes.PART_FINISHED
                        case StepRes.NOT_FINISHED:
                            return StepRes.NOT_FINISHED
        return StepRes.NOT_FINISHED

    def add_quads(self):
        for segment in self.segments:
            if not segment.finished:
                self.work_with_segment(segment)

    '''def search_vert(self, bm: bmesh.types.BMesh, point: Point) -> bmesh.types.BMVert:
        for v in bm.verts:
            if same_coords(v.co, self.reduced_points[point.get_i()]):
                return v
        raise ValueError("vert not found!")'''

    '''def render_quads(self, bm: bmesh.types.BMesh, nedges:int):
        for quad in self.quads:
            quad.render_quad(bm, nedges, self)'''


'''def work_with_mesh(glist: GlobalList, active: bpy.types.Object, nedges: int):
    rotation = active.rotation_euler
    location = active.location
    scale = active.scale
    mesh = bpy.data.meshes.new_from_object(active)

    bm = bmesh.new()
    bm.from_mesh(mesh)

    bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=TH)

    glist.render_quads(bm, nedges)

    bm.to_mesh(mesh)
    bm.free()
    obj = bpy.data.objects.new("gridfill", mesh)
    obj.rotation_euler = rotation
    obj.location = location
    obj.scale = scale
    bpy.context.collection.objects.link(obj)'''

glist = GlobalList()

'''active = bpy.context.active_object
cur = active.to_spline(bpy.context.evaluated_depsgraph_get())
nedges = cur.resolution_u
splines = cur.splines
mat = active.matrix_world'''
nedges = 16

'''for s in splines:
    spline = Spline(glist)
    for p in s.bezier_points:
        spline.add_point(mat @ p.co)
    if s.use_cyclic_u:
        spline.round_spline'''
'''ss = []
for s in splines:
    cos = []
    for p in s.bezier_points:
        cos.append(" ".join(f"{c}" for c in p.co))
    ss.append(":".join(cos))
DE = "|".join(ss)
raise ValueError(DE)
'''

for s in DE.split("|"):
    spline = Spline(glist)
    for p in s.split(":"):
        co = [float(aa) for aa in p.split(" ")]
        spline.add_point(co)
    #if s.use_cyclic_u:
        #spline.round_spline
for spline in glist.splines:
    print("spline")
    for point in spline.points:
        print(point.red_index)
glist.add_quads()
for quad in glist.quads:
    print("quad!!!")
#work_with_mesh(glist, active, nedges)
