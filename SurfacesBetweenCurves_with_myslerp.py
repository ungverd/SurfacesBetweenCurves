bl_info = {
    "name": "Create Surfaces Between Curves",
    "blender": (3, 4, 0),
    "category": "Object",
}


import bpy
import mathutils
from typing import List, Optional, Tuple, Callable, Literal
from enum import Enum
import math

TH = 0.00001
TH2 = TH**2

def quad_edges_to_normal(co_a1: mathutils.Vector,
                         co_a2: mathutils.Vector,
                         co_b1: mathutils.Vector,
                         co_b2: mathutils.Vector):
    diff_a = co_a2 - co_a1
    diff_b = co_b2 - co_b1
    diff_a.normalize()
    diff_b.normalize()
    no = diff_a + diff_b
    no.normalize()
    return no

def quad_verts_to_barycentric_tri(co_a: mathutils.Vector,
                                  co_b: mathutils.Vector,
                                  co_a_next: mathutils.Vector,
                                  co_b_next: mathutils.Vector,
                                  co_a_prev: mathutils.Vector | None,
                                  co_b_prev: mathutils.Vector | None,
                                  is_flip: bool):
    tri = [co_a, co_b]
    no = quad_edges_to_normal(co_a, co_a_next, co_b, co_b_next)
    if co_a_prev is not None and co_b_prev is not None:
        no_t = quad_edges_to_normal(co_a_prev, co_a, co_b_prev, co_b)
        no += no_t
        no.normalize()
    if is_flip:
        no = -no
    d = (co_a - co_b).length
    no *= d
    tri.append((co_a + co_b)/2)
    tri[2] += no
    return tri

def normal_tri_v3(v1: mathutils.Vector,
                  v2: mathutils.Vector,
                  v3: mathutils.Vector):
    n1 = v1 - v2
    n2 = v2 - v3
    n = n1.cross(n2)
    n.normalize()
    return n

def ortho_basis_v3v3_v3(n: mathutils.Vector):
    eps = 1.192092896e-07
    f = n[0] * n[0] + n[1] * n[1]
    if f > eps:
        d = 1 / (f ** 0.5)
        r_n1 = [n[1] * d, -n[0] * d, 0]
        r_n2 = [-n[2] * r_n1[1],
                n[2] * r_n1[0],
                n[0] * r_n1[1] - n[1] * r_n1[0]]
    else:
        r_n1 = [-1 if n[2] < 0 else 1, 0, 0]
        r_n2 = [0, 1, 0]
    return r_n1, r_n2

def axis_dominant_v3_to_m3(normal: mathutils.Vector):
    r_n1, r_n2 = ortho_basis_v3v3_v3(normal)
    r_mat = mathutils.Matrix((r_n1, r_n2, normal))
    r_mat.transpose()
    return r_mat

def cross_tri_v2(v1: mathutils.Vector,
                 v2: mathutils.Vector,
                 v3: mathutils.Vector):
    return (v1[0] - v2[0]) * (v2[1] - v3[1]) + (v1[1] - v2[1]) * (v3[0] - v2[0])

def barycentric_weights_v2(v1: mathutils.Vector,
                           v2: mathutils.Vector,
                           v3: mathutils.Vector,
                           co: mathutils.Vector):
    w = mathutils.Vector((cross_tri_v2(v2, v3, co),
                          cross_tri_v2(v3, v1, co),
                          cross_tri_v2(v1, v2, co)))
    wtot = w[0] + w[1] + w[2]
    if wtot != 0:
        w *= (1 / wtot)
    else:
        w = mathutils.Vector((1/3, 1/3, 1/3))
    return w

def interp_v3_v3v3v3(v1: mathutils.Vector,
                     v2: mathutils.Vector,
                     v3: mathutils.Vector,
                     w: mathutils.Vector):
    return mathutils.Vector((
        v1[0] * w[0] + v2[0] * w[1] + v3[0] * w[2],
        v1[1] * w[0] + v2[1] * w[1] + v3[1] * w[2],
        v1[2] * w[0] + v2[2] * w[1] + v3[2] * w[2]))

def interp_v3_v3v3(a: mathutils.Vector,
                   b: mathutils.Vector,
                   t: float):
    s = 1 - t
    return s * a + t * b

def transform_point_by_tri_v3(pt_src: mathutils.Vector,
                              tri_tar_p1: mathutils.Vector,
                              tri_tar_p2: mathutils.Vector,
                              tri_tar_p3: mathutils.Vector,
                              tri_src_p1: mathutils.Vector,
                              tri_src_p2: mathutils.Vector,
                              tri_src_p3: mathutils.Vector):
    no_tar = normal_tri_v3(tri_tar_p1, tri_tar_p2, tri_tar_p3)
    no_src = normal_tri_v3(tri_src_p1, tri_src_p2, tri_src_p3)
    mat_src = axis_dominant_v3_to_m3(no_src)
    pt_src_xy = pt_src @ mat_src
    tri_xy_src = (tri_src_p1 @ mat_src,
                  tri_src_p2 @ mat_src,
                  tri_src_p3 @ mat_src)
    w_src = barycentric_weights_v2(tri_xy_src[0], tri_xy_src[1], tri_xy_src[2], pt_src_xy)
    pt_tar = interp_v3_v3v3v3(tri_tar_p1, tri_tar_p2, tri_tar_p3, w_src)
    area_tar = math.sqrt(mathutils.geometry.area_tri(tri_tar_p1, tri_tar_p2, tri_tar_p3))
    area_src = math.sqrt(mathutils.geometry.area_tri(tri_xy_src[0], tri_xy_src[1], tri_xy_src[2]))
    z_ofs_src = pt_src_xy[2] - tri_xy_src[0][2]
    return pt_tar + no_tar * (z_ofs_src / area_src) * area_tar

def XY(x: int, y: int, xtot: int):
    return x + y * xtot

def grid_fill(verts1: List[mathutils.Vector],
              verts2: List[mathutils.Vector],
              railverts1: List[mathutils.Vector],
              railverts2: List[mathutils.Vector]):
    xtot = len(verts1)
    ytot = len(railverts1)
    v0 = mathutils.Vector((0, 0, 0))
    v_grid = verts1 + ([v0] * xtot * (ytot-2)) + verts2
    for i in range(1, ytot-1):
        v_grid[xtot*i] = railverts1[i]
        v_grid[xtot*(i+1) - 1] = railverts2[i]
    tri_a = quad_verts_to_barycentric_tri(
        v_grid[XY(0, 0, xtot)],
        v_grid[XY(xtot-1, 0, xtot)],
        v_grid[XY(0, 1, xtot)],
        v_grid[XY(xtot-1, 1, xtot)],
        None,
        None,
        False
    )
    tri_b = quad_verts_to_barycentric_tri(
        v_grid[XY(0, ytot-1, xtot)],
        v_grid[XY(xtot-1, ytot-1, xtot)],
        v_grid[XY(0, ytot-2, xtot)],
        v_grid[XY(xtot-1, ytot-2, xtot)],
        None,
        None,
        True
    )
    for y in range(1, ytot-1):
        tri_t = quad_verts_to_barycentric_tri(
            v_grid[XY(0, y, xtot)],
            v_grid[XY(xtot-1, y, xtot)],
            v_grid[XY(0, y+1, xtot)],
            v_grid[XY(xtot-1, y+1, xtot)],
            v_grid[XY(0, y-1, xtot)],
            v_grid[XY(xtot-1, y-1, xtot)],
            False
        )
        for x in range(1, xtot-1):
            co_a = transform_point_by_tri_v3(v_grid[x],
                                             tri_t[0],
                                             tri_t[1],
                                             tri_t[2],
                                             tri_a[0],
                                             tri_a[1],
                                             tri_a[2])
            co_b = transform_point_by_tri_v3(v_grid[(xtot * ytot) + (x - xtot)],
                                             tri_t[0],
                                             tri_t[1],
                                             tri_t[2],
                                             tri_b[0],
                                             tri_b[1],
                                             tri_b[2])
            co = interp_v3_v3v3(co_a, co_b, y / (ytot - 1))
            v_grid[(y * xtot) + x] = co
    return v_grid

#*******************************************************************************************

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
    def __init__(self,
                 i: int,
                 spline: "Spline",
                 bpoint: "BigPoint",
                 handle_left: mathutils.Vector,
                 handle_right: mathutils.Vector):
        self.bpoint: BigPoint = bpoint
        self.i: int = i
        self.prev_seg: Optional[Segment] = None
        self.post_seg: Optional[Segment] = None
        self.spline: Spline = spline
        self.handle_left = handle_left
        self.handle_right = handle_right


class Spline:
    def __init__(self, glist: "GlobalList"):
        self.points: List[Point] = []
        self.segments: List[Segment] = []
        self.glist = glist
        self.glist.add_spline(self)

    def add_point(self, coords: mathutils.Vector, handle_left: mathutils.Vector, handle_right: mathutils.Vector):
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
        point = Point(i-1, self, bpoint, handle_left - coords, handle_right - coords)
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
        self.quads: List[Quad | BigQuad] = []
        self.p1 = p1
        self.p2 = p2
        self.p1.post_seg = self
        self.p2.prev_seg = self
        self.finished = False
        glist.add_segment(self)

    def add_quad(self, quad: "Quad | BigQuad"):
        self.quads.append(quad)
        if len(self.quads) == 2:
            self.finished = True

    @staticmethod
    def iterate_direct(segs: "List[Segment]"):
        yield segs[0].p1
        for seg in segs:
            yield seg.p2

    @staticmethod
    def iterate_reversed(segs: "List[Segment]"):
        yield segs[-1].p2
        for seg in reversed(segs):
            yield seg.p1

def get_y_normalized(p: mathutils.Vector, x: mathutils.Vector):
    return (p - p.project(x)).normalized()

def calc_basis(init: mathutils.Vector,
               edge: mathutils.Vector,
               point: mathutils.Vector,
               handle1: mathutils.Vector,
               handle2: Optional[mathutils.Vector] = None):
    v0 = mathutils.Vector((0, 0, 0))
    x_vec = (edge - init).normalized()
    z_vec = x_vec.cross(handle1)
    if z_vec.length_squared > TH2:
        y_vec = get_y_normalized(handle1, x_vec)
        z_vec.normalize()
    else:
        p_vec = point - init
        z_vec = x_vec.cross(p_vec)
        if (z_vec.length_squared > TH2):
            y_vec = get_y_normalized(p_vec, x_vec)
            z_vec.normalize()
        else:
            if (handle2 is None):
                y_vec = v0
                z_vec = v0
            else:
                z_vec = x_vec.cross(handle2)
                if (z_vec.length_squared > TH2):
                    z_vec.normalize()
                    y_vec = get_y_normalized(handle2, x_vec)
                else:
                    y_vec = v0
                    z_vec = v0
    return mathutils.Matrix((x_vec, y_vec, z_vec))

def get_coords_from_vec_and_basis_matrix(v: mathutils.Vector, m: mathutils.Matrix) -> mathutils.Vector:
    res = m @ v
    assert isinstance(res, mathutils.Vector)
    return res

def get_vec_from_coords_and_basis_matrix(v: mathutils.Vector, m: mathutils.Matrix) -> mathutils.Vector:
    return v @ m

def my_slerp(p1: mathutils.Vector, p2: mathutils.Vector, factor: float):
    # slerp returns unit vector, we interpolate length separately
    return p1.slerp(p2, factor) * (p1.length * (1 - factor) + p2.length * factor)

def make_collinear(v1: mathutils.Vector, v2: mathutils.Vector):
    dif = (v1 - v2).normalized()
    new_v1 = dif * v1.length
    new_v2 = -dif * v2.length
    return new_v1, new_v2


class MiddlePoint:
    def __init__(self, co: mathutils.Vector):
        self.co = co
        self.handle_up: mathutils.Vector
        self.handle_down: mathutils.Vector
        self.handle_left: mathutils.Vector
        self.handle_right: mathutils.Vector


class BigQuad:
    def __init__(self,
                 edg1: List[Segment],
                 edg2: List[Segment],
                 edg3: List[Segment],
                 edg4: List[Segment],
                 dir1: bool,
                 dir2: bool,
                 dir3: bool,
                 dir4: bool,
                 glist: "GlobalList"):
        self.v_grid: List[mathutils.Vector]
        self.xtot = len(edg1) + 1
        self.ytot = len(edg2) + 1
        self.edges = [edg1, edg2, edg3, edg4]
        self.dirs = [dir1, dir2, dir3, dir4]
        self.points_grid: List[List[MiddlePoint]] = []
        self.corner_handles: List[Tuple[mathutils.Vector, mathutils.Vector]] = []
        for edg in self.edges:
            for segment in edg:
                segment.add_quad(self)
        glist.add_big_quad(self)

    @staticmethod
    def extract_first_handle(edge: List[Segment]):
        return edge[0].p1.handle_right

    @staticmethod
    def extract_last_handle(edge: List[Segment]):
        return edge[-1].p2.handle_left

    def extract_corner_handles(self):
        # edges[0]
        if self.dirs[1]:
            handle_right_0 = BigQuad.extract_first_handle(self.edges[1])
        else:
            handle_right_0 = BigQuad.extract_last_handle(self.edges[1])
        if self.dirs[3]:
            handle_left_0 = BigQuad.extract_first_handle(self.edges[3])
        else:
            handle_left_0 = BigQuad.extract_last_handle(self.edges[3])
        self.corner_handles.append((handle_left_0, handle_right_0))

        # edges[1]
        handle_left_1 = BigQuad.extract_last_handle(self.edges[0])
        if self.dirs[2]:
            handle_right_1 = BigQuad.extract_last_handle(self.edges[2])
        else:
            handle_right_1 = BigQuad.extract_first_handle(self.edges[2])
        if not self.dirs[1]:
            handle_left_1, handle_right_1 = handle_right_1, handle_left_1
        self.corner_handles.append((handle_left_1, handle_right_1))

        # edges[2]
        if self.dirs[1]:
            handle_right_2 = BigQuad.extract_last_handle(self.edges[1])
        else:
            handle_right_2 = BigQuad.extract_first_handle(self.edges[1])
        if self.dirs[3]:
            handle_left_2 = BigQuad.extract_last_handle(self.edges[3])
        else:
            handle_left_2 = BigQuad.extract_first_handle(self.edges[3])
        if not self.dirs[2]:
            handle_left_2, handle_right_2 = handle_right_2, handle_left_2
        self.corner_handles.append((handle_left_2, handle_right_2))

        # edges[3]
        handle_left_3 = BigQuad.extract_first_handle(self.edges[0])
        if self.dirs[2]:
            handle_right_3 = BigQuad.extract_first_handle(self.edges[2])
        else:
            handle_right_3 = BigQuad.extract_last_handle(self.edges[2])
        if not self.dirs[3]:
            handle_left_3, handle_right_3 = handle_right_3, handle_left_3
        self.corner_handles.append((handle_left_3, handle_right_3))

    @staticmethod
    def extract_coords(edg: List[Segment], dir: bool, glist: "GlobalList"):
        if dir:
            res = [glist.get_coords(p.i) for p in Segment.iterate_direct(edg)]
        else:
            res = [glist.get_coords(p.i) for p in Segment.iterate_reversed(edg)]
        return res

    def extract_handles(self, edg: List[Segment], i: int):
        handles: List[Optional[mathutils.Vector]] = [None] * (len(edg) + 1)
        handles[0] = self.corner_handles[i][0]
        handles[-1] = self.corner_handles[i][1]
        for j, p in enumerate(Segment.iterate_direct(edg)):
            if j not in (0, len(edg)):
                if len(p.bpoint.points) != 1:
                    point_handles: List[mathutils.Vector] = []
                    for pp in p.bpoint.points:
                        if pp != p:
                            if pp.post_seg and not pp.prev_seg:
                                point_handles.append(pp.handle_left)
                            elif pp.prev_seg and not pp.post_seg:
                                point_handles.append(pp.handle_right)
                    if len(point_handles) != 0:
                        if len(point_handles) == 1:
                            handles[j] = point_handles[0]
                        else:
                            total = point_handles[0]
                            for h in point_handles[1:]:
                                total += h
                            handles[j] = total/len(point_handles)
        return handles

    @staticmethod
    def populate_handles(handles: List[Optional[mathutils.Vector]]) -> List[mathutils.Vector]:
        i = 0
        count = 1
        while i < (len(handles) - 1):
            p1 = handles[i]
            assert isinstance(p1, mathutils.Vector)
            count = 1
            p2 = handles[i + count]
            while p2 is None:
                count += 1
                p2 = handles[i + count]
            for j in range(i + 1, i + count):
                factor = j / count
                # slerp returns unit vector, we interpolate length separately
                handles[j] = my_slerp(p1, p2, factor)
            i += count
        res_handles = [h for h in handles if h is not None]
        assert len(res_handles) == len(handles)
        return res_handles

    def calc_init_coords(self,
                         i: int,
                         right: bool,
                         horisontal: bool,
                         h1: bool,
                         h: List[mathutils.Vector]):
        tot = self.ytot if horisontal else self.xtot
        second_coord = 0 if h1 else tot - 1
        first_coord1 = i - 1 if right else i + 1
        first_coord2 = i + 1 if right else i - 1
        coords1 = [first_coord1, second_coord]
        coords2 = [first_coord2, second_coord]
        coords3 = [i, second_coord]
        if not horisontal:
            for l in (coords1, coords2, coords3):
                l.reverse()
        for l in (coords1, coords2, coords3):
            l.append(self.xtot)
        matrix = calc_basis(self.v_grid[XY(*coords1)],
                            self.v_grid[XY(*coords2)],
                            self.v_grid[XY(*coords3)],
                            h[first_coord1],
                            h[i])
        return get_coords_from_vec_and_basis_matrix(h[i], matrix)

    def handle_fin(self,
                   i: int,
                   j: int,
                   vh1: List[mathutils.Vector],
                   vh2: List[mathutils.Vector],
                   coords1: mathutils.Vector,
                   coords2: mathutils.Vector,
                   right: bool,
                   horisontal: bool):
        xtot = self.xtot if horisontal else self.ytot
        ytot = self.ytot if horisontal else self.xtot
        coord1 = i - 1 if right else i + 1
        coord2 = i + 1 if right else i - 1
        c1 = [coord1, j]
        c2 = [coord2, j]
        c3 = [i, j]
        h_side = my_slerp(vh1[j], -vh2[j], (i - 1)/(xtot - 1)) if right else my_slerp(-vh1[j], vh2[j], (i + 1)/(xtot - 1))
        if not horisontal:
            for c in (c1, c2, c3):
                c.reverse()
        for c in (c1, c2, c3):
            c.append(self.xtot)
        dest_m = calc_basis(self.v_grid[XY(*c1)],
                            self.v_grid[XY(*c2)],
                            self.v_grid[XY(*c3)],
                            h_side)
        coords_fin = my_slerp(coords1, coords2, j/(ytot - 1))
        return get_vec_from_coords_and_basis_matrix(coords_fin, dest_m)

    def subdivide(self, glist: "GlobalList", obj: bpy.types.ID):
        v1 = BigQuad.extract_coords(self.edges[0], True, glist)
        v2 = BigQuad.extract_coords(self.edges[2], self.dirs[2], glist)
        rv1 = BigQuad.extract_coords(self.edges[3], self.dirs[3], glist)
        rv2 = BigQuad.extract_coords(self.edges[1], self.dirs[1], glist)
        self.v_grid = grid_fill(v1, v2, rv1, rv2)
        for i in range(self.ytot):
            row: List[MiddlePoint] = []
            for j in range(self.xtot):
                row.append(MiddlePoint(self.v_grid[XY(j, i, self.xtot)]))
            self.points_grid.append(row)
        self.extract_corner_handles()
        cross_hs = [BigQuad.populate_handles(self.extract_handles(edg, i)) for i, edg in enumerate(self.edges)]
        for i in range(4):
            if not self.dirs[i]:
                cross_hs[i].reverse()
        cross_h1 = cross_hs[0]
        cross_h2 = cross_hs[2]
        cross_vh1 = cross_hs[3]
        cross_vh2 = cross_hs[1]
        h_right_left = [
            (
                [p.handle_right for p in Segment.iterate_direct(edge)] if d else 
                [p.handle_left for p in Segment.iterate_reversed(edge)],
                [p.handle_left for p in Segment.iterate_direct(edge)] if d else
                [p.handle_right for p in Segment.iterate_reversed(edge)]
            ) for edge, d in zip(self.edges, self.dirs)]
        h1_right = h_right_left[0][0]
        h1_left = h_right_left[0][1]
        h2_right = h_right_left[2][0]
        h2_left= h_right_left[2][1]
        vh1_down = h_right_left[3][0]
        vh1_up = h_right_left[3][1]
        vh2_down = h_right_left[1][0]
        vh2_up = h_right_left[1][1]
        horisontal = True
        for i in range(1, self.xtot - 1):
            coords_right1 = self.calc_init_coords(i, True, horisontal, True, h1_right)
            coords_right2 = self.calc_init_coords(i, True, horisontal, False, h2_right)
            coords_left1 = self.calc_init_coords(i, False, horisontal, True, h1_left)
            coords_left2 = self.calc_init_coords(i, False, horisontal, False, h2_left)
            for j in range(1, self.ytot - 1):
                res_right = self.handle_fin(i, j, cross_vh1, cross_vh2, coords_right1, coords_right2, True, horisontal)
                res_left = self.handle_fin(i, j, cross_vh1, cross_vh2, coords_left1, coords_left2, False, horisontal)
                collinear_right, collinear_left = make_collinear(res_right, res_left)
                self.points_grid[j][i].handle_right = collinear_right
                self.points_grid[j][i].handle_left = collinear_left
        vertical = False
        for j in range(1, self.ytot - 1):
            coords_down1 = self.calc_init_coords(j, True, vertical, True, vh1_down)
            coords_down2 = self.calc_init_coords(j, True, vertical, False, vh2_down)
            coords_up1 = self.calc_init_coords(j, False, vertical, True, vh1_up)
            coords_up2 = self.calc_init_coords(j, False, vertical, False, vh2_up)
            for i in range(1, self.xtot - 1):
                res_down = self.handle_fin(j, i, cross_h1, cross_h2, coords_down1, coords_down2, True, vertical)
                res_up = self.handle_fin(j, i, cross_h1, cross_h2, coords_up1, coords_up2, False, vertical)
                collinear_down, collinear_up = make_collinear(res_down, res_up)
                self.points_grid[j][i].handle_down = collinear_down
                self.points_grid[j][i].handle_up = collinear_up
        for i in range(self.xtot):
            self.points_grid[0][i].handle_down = cross_h1[i]
            self.points_grid[0][i].handle_up = -cross_h1[i]
            self.points_grid[self.ytot - 1][i].handle_down = -cross_h2[i]
            self.points_grid[self.ytot - 1][i].handle_up = cross_h2[i]
        for i in range(self.ytot):
            self.points_grid[i][0].handle_right = cross_vh1[i]
            self.points_grid[i][0].handle_left = -cross_vh1[i]
            self.points_grid[i][self.xtot - 1].handle_right = -cross_vh2[i]
            self.points_grid[i][self.xtot - 1].handle_left = cross_vh2[i]
        for i in range(1, self.xtot - 1):
            spline = obj.data.splines.new("BEZIER")
            spline.bezier_points.add(self.ytot - 1)
            for j, point in zip(range(self.ytot), spline.bezier_points):
                p = self.points_grid[j][i]
                point.co = p.co
                point.handle_left = p.handle_up + p.co
                point.handle_right = p.handle_down + p.co
                point.handle_left_type = "ALIGNED"
                point.handle_right_type = "ALIGNED"
        for j in range(1, self.ytot - 1):
            spline = obj.data.splines.new("BEZIER")
            spline.bezier_points.add(self.xtot - 1)
            for i, point in zip(range(self.xtot), spline.bezier_points):
                p = self.points_grid[j][i]
                point.co = p.co
                point.handle_left = p.handle_left + p.co
                point.handle_right = p.handle_right + p.co
                point.handle_left_type = "ALIGNED"
                point.handle_right_type = "ALIGNED"


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
    def turn(p: Point):
        for pp in p.bpoint.points:
            if pp != p:
                yield pp

    @staticmethod
    def little_quads_step_left(p: Point):
        if p.prev_seg:
            return p.prev_seg.p1
        return None

    @staticmethod
    def little_quads_step_right(p: Point):
        if p.post_seg:
            return p.post_seg.p2
        return None

    @staticmethod
    def step(p: Point,
            func: Callable[[Point], Optional[Point]],
            s1: int,
            s2: int,
            bp_s1: BigPoint,
            bp_s2: BigPoint) -> bool | Point:
        res = func(p)
        if res:
            for pp in Quad.turn(res):
                for s, bp in zip((s1, s2), (bp_s1, bp_s2)):
                    if s != 0:
                        for f in (Quad.steps_dividing_edges_left, Quad.steps_dividing_edges_right):
                            if f(pp, s) == bp:
                                return True
            return res
        return False

    @staticmethod
    def steps(p: Point,
              n: int,
              end_bps_s1: List[BigPoint],
              end_bps_s2: List[BigPoint],
              s1: int,
              s2: int,
              func: Callable[[Point], Optional[Point]]) -> bool:
        res = p
        for i in range(n):
            if i != (n - 1) or (s1 != 0 and s2 != 0):
                res = Quad.step(res, func, s1, s2, end_bps_s1[i], end_bps_s2[i])
                if res in (True, False):
                    return res
        return False

    @staticmethod
    def verify_little_quads(p1: Point,
                            end_bps_s1: List[BigPoint],
                            end_bps_s2: List[BigPoint],
                            s1: int,
                            s2: int,
                            n: int):
        for p in Quad.turn(p1):
            for func in (Quad.little_quads_step_left, Quad.little_quads_step_right):
                if Quad.steps(p, n, end_bps_s1, end_bps_s2, s1, s2, func):
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
            if dir1 == dir2:
                if Quad.verify_dividing_edges(b1[i], b2[i], n):
                    return True
            else:
                if Quad.verify_dividing_edges(b1[i], b2[-i-1], n):
                    return True
        return False

    @staticmethod
    def extract_points(edge: List[Segment], dir: bool) -> List[Point]:
        points: List[Point] = []
        if edge:
            if dir:
                points = list(Segment.iterate_direct(edge))
            else:
                points = list(Segment.iterate_reversed(edge))
        return points

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
        dir4 = not Quad.get_dir_from_left(edg1, segments)
        if dir2:
            dir3 = not Quad.get_dir_from_right(edg2, segments)
        else:
            dir3 = Quad.get_dir_from_left(edg2, segments)
        if le == 4:
            Quad(edg1[0], edg2[0], edg3[0], edg4[0], dir1, dir2, dir3, dir4, glist)
            return True
        if Quad.verify_div(b1, b3, dir1, dir3, len(edg2)) or Quad.verify_div(b2, b4, dir2, dir4, len(edg1)):
            return False
        points1 = Quad.extract_points(edg1, True)
        points2 = Quad.extract_points(edg2, dir2)
        bpoints2 = [p.bpoint for p in points2]
        points3 = Quad.extract_points(edg3, dir3)
        points4 = Quad.extract_points(edg4, dir4)
        bpoints4 = [p.bpoint for p in points4]
        le_e1 = len(edg1)
        n = len(edg2)
        for i, p in enumerate(points1):
            if Quad.verify_little_quads(p, bpoints4[1:], bpoints2[1:], i, le_e1 - i, n):
                return False
        if Quad.steps(points1[0], n, bpoints2[1:], bpoints2[1:], le_e1, 0, Quad.little_quads_step_left):
            return False
        if Quad.steps(points1[-1], n, bpoints4[1:], bpoints4[1:], le_e1, 0, Quad.little_quads_step_right):
            return False
        bpoints2.reverse()
        bpoints4.reverse()
        for i, p in enumerate(points3):
            if Quad.verify_little_quads(p, bpoints4[1:], bpoints2[1:], i, le_e1 - i, n):
                return False
        if Quad.steps(points3[0], n, bpoints2[1:], bpoints2[1:], le_e1, 0,
                      Quad.little_quads_step_right if dir3 else Quad.little_quads_step_left):
            return False
        if Quad.steps(points3[-1], n, bpoints4[1:], bpoints4[1:], le_e1, 0,
                      Quad.little_quads_step_left if dir3 else Quad.little_quads_step_right):
            return False
        BigQuad(edg1, edg2, edg3, edg4, dir1, dir2, dir3, dir4, glist)
        return True

    def __init__(self,
                 seg1: Segment,
                 seg2: Segment,
                 seg3: Segment,
                 seg4: Segment,
                 dir1: bool,
                 dir2: bool,
                 dir3: bool,
                 dir4: bool,
                 glist: "GlobalList"):
        self.segments = [seg1, seg2, seg3, seg4]
        self.directions = [dir1, dir2, dir3, dir4]
        for seg in self.segments:
            seg.add_quad(self)
        glist.add_quad(self)


class GlobalList:
    def __init__(self):
        self.reduced_points: List[mathutils.Vector] = []
        self.big_points: List[BigPoint] = []
        self.count = 0
        self.splines: List[Spline] = []
        self.quads: List[Quad] = []
        self.big_quads: List[BigQuad] = []
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

    def add_big_quad(self, big_quad: BigQuad):
        self.big_quads.append(big_quad)

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

    def subdivide_quads(self, obj: bpy.types.ID):
        for bq in self.big_quads:
            bq.subdivide(self, obj)


def work_with_mesh(glist: GlobalList, active: bpy.types.Object, nedges: int):
    '''rotation = active.rotation_euler
    location = active.location
    scale = active.scale

    obj = bpy.data.objects.new(f"{active.name}.gridfill", mesh)
    obj.rotation_euler = rotation
    obj.location = location
    obj.scale = scale
    bpy.context.collection.objects.link(obj)'''


class CreateSurfacesBetweenCurves(bpy.types.Operator):
    """My Creating Surfaces Between Curves Script"""      # Use this as a tooltip for menu items and buttons.
    bl_idname = "object.create_surfs"        # Unique identifier for buttons and menu items to reference.
    bl_label = "Create surfaces"         # Display name in the interface.
    bl_options = {'REGISTER', 'UNDO'}  # Enable undo for the operator.

    def execute(self, context):        # execute() is called when running the operator.

        glist = GlobalList()

        active = context.active_object
        cur = active.data
        #nedges = cur.resolution_u
        splines = cur.splines

        for s in splines:
            spline = Spline(glist)
            for p in s.bezier_points:
                spline.add_point(p.co, p.handle_left, p.handle_right)
            if s.use_cyclic_u:
                spline.round_spline()
        glist.add_quads()
        obj = active.copy()
        obj.data = active.data.copy()
        glist.subdivide_quads(obj)
        context.collection.objects.link(obj)

        #work_with_mesh(glist, active, nedges)

        return {'FINISHED'}            # Lets Blender know the operator finished successfully.

def menu_func(self, context):
    self.layout.operator(CreateSurfacesBetweenCurves.bl_idname)

def register():
    bpy.utils.register_class(CreateSurfacesBetweenCurves)
    bpy.types.VIEW3D_MT_object.append(menu_func)  # Adds the new operator to an existing menu.

def unregister():
    bpy.utils.unregister_class(CreateSurfacesBetweenCurves)


# This allows you to run the script directly from Blender's Text editor
# to test the add-on without having to install it.
if __name__ == "__main__":
    register()

