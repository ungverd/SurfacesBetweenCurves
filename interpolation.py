import mathutils
from typing import List

def quad_edges_to_normal(co_a1: mathutils.Vector,
                         co_a2: mathutils.Vector,
                         co_b1: mathutils.Vector,
                         co_b2: mathutils.Vector):
    diff_a = co_a2 - co_a1
    diff_b = co_b2 - co_b1
    diff_a.normalize()
    diff_b.normalize()
    return diff_a + diff_b

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
            True
        )
        for x in range(1, xtot-1):

