import mathutils

def quad_edges_to_normal(co_a1: mathutils.Vector,
                         co_a2: mathutils.Vector,
                         co_b1: mathutils.Vector,
                         co_b2: mathutils.Vector):
    diff_a = co_a2 - co_a1
    diff_b = co_b2 - co_b1
    diff_a = diff_a.normalized()
    diff_b = diff_b.normalized()
    return diff_a + diff_b

def quad_verts_to_barycentric_tri(co_a: mathutils.Vector,
                                  co_b: mathutils.Vector,
                                  co_a_next: mathutils.Vector,
                                  co_b_next: mathutils.Vector,
                                  co_a_prev: mathutils.Vector,
                                  co_b_prev: mathutils.Vector,
                                  flip: bool):
    tri = [co_a, co_b]
    no = quad_edges_to_normal(co_a, co_a_next, co_b, co_b_next)
    if co_a_prev:
        no_t = quad_edges_to_normal(co_a_prev, co_a, co_b_prev, co_b)
        no += no_t
        no = no.normalized()
    if is_flip:
        no = -no
    d = (co_a - co_b).length
    no *= d
    tri.append((co_a + co_b)/2)
    tri[2] += no
    return tri

