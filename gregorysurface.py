import bpy
import bmesh
import mathutils
import math
from typing import List
# from enum import Enum
import time
from copy import deepcopy

N = 64 # number of subdivisions
N1 = N + 1

TH = 0.00001
TH2 = TH**2

def calculate_faces():
    faces: List[List[int]] = []
    for i in range(N-2):
        for j in range(N-2):
            face = [i * (N - 1) + j]
            face.append(i * (N - 1) + j + 1)
            face.append((i + 1) * (N - 1) + j + 1)
            face.append((i + 1) * (N - 1) + j)
            faces.append(face)
    return faces

'''def calculate_faces():
    faces: List[List[int]] = []
    for i in range(N):
        for j in range(N):
            face = [i * (N + 1) + j]
            face.append(i * (N + 1) + j + 1)
            face.append((i + 1) * (N + 1) + j + 1)
            face.append((i + 1) * (N + 1) + j)
            faces.append(face)
    return faces'''

FACES = calculate_faces()

try:
    #raise(ModuleNotFoundError)
    import numpy as np
    import numpy.typing as npt
    t = time.time()
    J_S = np.arange(4)                               # shape (4)
    COMBS = np.array([math.comb(3, i) for i in J_S]) # shape (4)
    V = np.linspace(0, 1, N1) # shape      (N+1)
    U = np.expand_dims(V, 1)  # shape (N+1, 1  )
    V0 = V[1:N] # shape      (N-1)
    U0 = U[1:N] # shape (N-1, 1  )
    V02 = np.expand_dims(np.stack((V0, 1-V0), 1), axis=(1,3)) # shape      (N-1, 1, 2, 1)
    U02 = np.expand_dims(V02, axis=4)                         # shape (N-1, 1  , 2, 1, 1)
    UV_DIV = U02 + V02                                        # shape (N-1, N-1, 2, 2, 1)

    '''BERNS_NP = np.expand_dims(COMBS * U**J_S * (1-U)**(3-J_S), axis=1) # shape      (N+1, 1, 4)
    BERNS_NP_2 = np.expand_dims(BERNS_NP, axis=3)                      # shape (N+1, 1  , 4, 1)
    BERNS2 = np.expand_dims(BERNS_NP*BERNS_NP_2, axis=4)               # shape (N+1, N+1, 4, 4, 1)'''

    BERNS_NP = np.expand_dims(COMBS * U0**J_S * (1-U0)**(3-J_S), axis=1) # shape      (N-1, 1, 4)
    BERNS_NP_2 = np.expand_dims(BERNS_NP, axis=3)                        # shape (N-1, 1  , 4, 1)
    BERNS2 = np.expand_dims(BERNS_NP*BERNS_NP_2, axis=4)                 # shape (N-1, N-1, 4, 4, 1)

    # numpy magic to get products of all combinations of Bernstein coefficients

    DIMS1 = np.empty((N-1, N-1, 1, 2, 3))
    DIMS2 = np.empty((N-1, N-1, 4, 1, 3))
    #DIMS3 = np.empty((1, N-1, 4, 4, 3))
    #DIMS4 = np.empty((N+1, 1, 4, 4, 3))

    def calc_control_points_np(input1: List[List[mathutils.Vector]],
                               input2: List[List[mathutils.Vector]]):
        full = np.array(input1)                     # shape           (4, 4, 3)
        cp_u = full[1:3, 1:3, :]                    # shape           (2, 2, 3)
        cp_v = np.array(input2)                     # shape           (2, 2, 3)
        cp_central = (cp_u*V02 + cp_v*U02) / UV_DIV # shape (N-1, N-1, 2, 2, 3)

        top, _ = np.broadcast_arrays(full[0, 1:3], DIMS1)        # shape (N-1, N-1, 1, 2, 3)
        bottom, _ = np.broadcast_arrays(full[3, 1:3], DIMS1)     # shape (N-1, N-1, 1, 2, 3)
        res1 = np.concatenate((top, cp_central, bottom), axis=2) # shape (N-1, N-1, 4, 2, 3)
        left, _ = np.broadcast_arrays(np.expand_dims(full[:, 0], axis=1), DIMS2)  # shape (N-1, N-1, 4, 1, 3)
        right, _ = np.broadcast_arrays(np.expand_dims(full[:, 3], axis=1), DIMS2) # shape (N-1, N-1, 4, 1, 3)
        res2 = np.concatenate((left, res1, right), axis=3)  # shape (N-1, N-1, 4, 4, 3)
        '''side1, _ = np.broadcast_arrays(full, DIMS3)         # shape (1  , N-1, 4, 4, 3)
        res3 = np.concatenate((side1, res2, side1), axis=0) # shape (N+1, N-1, 4, 4, 3)
        side2, _ = np.broadcast_arrays(full, DIMS4)         # shape (N+1, 1  , 4, 4, 3)
        return np.concatenate((side2, res3, side2), axis=1) # shape (N+1, N+1, 4, 4, 3)'''
        return res2

    def calc_gregory_surf(kk: List[List[mathutils.Vector]],
                          kk1: List[List[mathutils.Vector]],
                          mesh: bpy.types.Mesh):
        control_points = calc_control_points_np(kk, kk1)
        #res: npt.NDArray[np.float32] = np.sum((BERNS2 * control_points), (2, 3)).reshape(N1 * N1, 3)
        res: npt.NDArray[np.float32] = np.sum((BERNS2 * control_points), (2, 3)).reshape(N-1 * N-1, 3)
        # numpy representation of the formula p(u,v) = sum_i_from_0_to_3(sum_j_from_0_to_3( k(i,j)*B(i,u)*B(j,v) ))
        mesh.from_pydata(res, [], FACES)

except (ModuleNotFoundError, ImportError):
    t = time.time()
    def bernstein(i: int, u: float):
        return math.comb(3, i) * u**i * (1-u)**(3-i)

    BERNS = [[bernstein(i, k/N) for i in range(4)] for k in range(N1)]
    def berns2_el(i: int, j: int, nu: int, nv:int):
        return BERNS[nu][i] * BERNS[nv][j]

    B2 = [[[[berns2_el(i,j,nu,nv) for j in range(4)] for i in range(4)] for nv in range(N1)] for nu in range(N1)]

    def calc_control_points(kk: List[List[mathutils.Vector]],
                            kk1: List[List[mathutils.Vector]],
                            nu: int,
                            nv: int):
        u = nu / N
        v = nv / N
        cp = deepcopy(kk)
        cp[1][1] = (u * kk1[0][0] + v * kk[1][1]) / (u + v)
        cp[2][1] = ((1 - u) * kk1[1][0] + v * kk[2][1]) / (1 - u + v)
        cp[1][2] = (u * kk1[0][1] + (1 - v) * kk[1][2]) / (1 - v + u)
        cp[2][2] = ((1 - u) * kk1[1][1] + (1 - v) * kk[2][2]) / (2 - u - v)
        return cp

    def calc_point(nu: int, nv: int, kk: List[List[mathutils.Vector]], kk1: List[List[mathutils.Vector]]):
        el: List[float] = []
        control_points = calc_control_points(kk, kk1, nu, nv)
        for d in range(3):
            el.append(sum(sum(control_points[i][j][d] * B2[nu][nv][i][j] for j in range(4)) for i in range(4)))
        return mathutils.Vector(el)

    def calc_gregory_surf(kk: List[List[mathutils.Vector]],
                          kk1: List[List[mathutils.Vector]],
                          mesh: bpy.types.Mesh):
        res = [calc_point(nu, nv, kk, kk1) for nv in range(1, N) for nu in range(1, N)]
        mesh.from_pydata(res, [], FACES)

def are_coplanar(v1: mathutils.Vector, v2: mathutils.Vector, v3: mathutils.Vector):
    return abs(mathutils.Matrix((v1, v2, v3)).determinant()) < TH

def calc_basis(edge: mathutils.Vector, point: mathutils.Vector, handle: mathutils.Vector):
    x_vec = edge
    y_vec = point - point.project(edge)
    z_vec = edge.cross(point)
    return (x_vec, y_vec, z_vec)

def make_coplanar(v1: mathutils.Vector, v2: mathutils.Vector, v3: mathutils.Vector):
    v4 = v2 - v3
    normal = v1.cross(v4)
    new_v2 = v2 - v2.project(normal)
    new_v3 = v3 - v3.project(normal)
    return (new_v2, new_v3)

class BezierCurve:
    def __init__(self,
                 p0: mathutils.Vector,
                 p1: mathutils.Vector,
                 p2: mathutils.Vector,
                 p3: mathutils.Vector):
        self.p0 = p0
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.c0 = p1 - p0
        self.c1 = p2 - p1
        self.c2 = p3 - p2

    def get_a0(self, curve: "BezierCurve"):
        return curve.c0 if same_coords(curve.p0, self.p0) else (-curve.c2)

    def get_a3(self, curve: "BezierCurve"):
        return curve.c0 if same_coords(curve.p0, self.p3) else (-curve.c2)

    def get_b0(self, curve: "BezierCurve"):
        return curve.c2 if same_coords(curve.p3, self.p0) else (-curve.c0)

    def get_b3(self, curve: "BezierCurve"):
        return curve.c2 if same_coords(curve.p3, self.p3) else (-curve.c0)

    def calc_Chiyokura_cp(self,
                          curve_a0: "BezierCurve",
                          curve_a3: "BezierCurve",
                          curve_b0: "BezierCurve",
                          curve_b3: "BezierCurve") -> bool:
        a0 = self.get_a0(curve_a0)
        a3 = self.get_a3(curve_a3)
        b0 = self.get_b0(curve_b0)
        b3 = self.get_b3(curve_b3)
        c0 = self.c0
        c1 = self.c1
        c2 = self.c2
        for v in (a0, a3, b0, b3, c0, c2):
            if v.length < TH:
                return False
        if not are_coplanar(a0, b0, c0) or not are_coplanar(a3, c2, b3):
            return False
        if b0.cross(c0).length < TH or b3.cross(c2).length < TH:
            return False
        k0, h0 = get_coefs(b0, c0, a0)
        k1, h1 = get_coefs(b3, c2, a3)
        b1 = 2/3 * b0 + 1/3 * b3
        b2 = 1/3 * b0 + 2/3 * b3
        a1 = 1/3 * (k1 - k0) * b0 + k0 * b1 + 2/3 * h0 * c1 + 1/3 * h1 * c0
        a2 = k1 * b2 - 1/3 * (k1 - k0) * b3 + 1/3 * h0 * c2 + 2/3 * h1 * c1
        return True


def get_coefs(b: mathutils.Vector, c: mathutils.Vector, a: mathutils.Vector):
    if b[0] < TH and c[0] < TH:
        return coefs_use_axis(b, c, a, 1, 2)
    if b[1] < TH and c[1] < TH:
        return coefs_use_axis(b, c, a, 0, 2)
    return coefs_use_axis(b, c, a, 0, 1)

def coefs_use_axis(b: mathutils.Vector, c: mathutils.Vector, a: mathutils.Vector, i: int, j: int):
    h = (a[i]*b[j] - a[j]*b[i]) / (c[i]*b[j] - c[j]*b[i])
    if b[i] > TH:
        k = (a[i] - h * c[i]) / b[i]
    else:
        k = (a[j] - h * c[j]) / b[j]
    return (k, h)


def get_loc(name: str) -> mathutils.Vector:
    return bpy.data.objects[name].location

def same_coords(c1: mathutils.Vector, c2: mathutils.Vector) -> bool:
    return (c1-c2).length_squared < TH2

kk = [[get_loc(f"k{j}{i}")for i in range(4)] for j in range(4)]

kk1 = deepcopy((kk)[1:3])
kk1 = [row[1:3] for row in kk1]
kk1[0][0] += mathutils.Vector((0.1, 0.05, 0.02))
kk1[0][1] += mathutils.Vector((0.1, 0.05, 0.02))
kk1[1][0] += mathutils.Vector((0.1, 0.05, 0.02))
kk1[1][1] += mathutils.Vector((0.1, 0.05, 0.02))

'''kk:np.ndarray = np.array([
               [np.array([-1.        ,  0.        ,  2.15774703]),
                np.array([-1.        ,  1.48469925,  2.15774703]),
                np.array([ 1.        , -1.45185316,  2.15774703]),
                np.array([1.        , 0.        , 2.15774703])],
               [np.array([-0.5       ,  0.5       ,  2.15774703]),
                np.array([-0.54151225, -0.01831961,  1.25261807]),
                np.array([0.20227045, 0.15929604, 1.22609735]),
                np.array([ 0.65001631, -0.07832319,  1.00388753])],
               [np.array([-1.88212955, -0.26279324, -0.45193356]),
                np.array([-0.56988615, -0.42863035,  0.45224738]),
                np.array([ 0.05436653, -0.26342535,  0.4620285 ]),
                np.array([ 0.        ,  0.        , -0.93721968])],
               [np.array([-1.        ,  0.        , -0.93721968]),
                np.array([-0.93310571,  0.9774524 , -1.01902032]),
                np.array([ 1.        , -1.        , -0.93721968]),
                np.array([ 1.        ,  0.        , -0.93721968])]])'''

'''kk = [
               [[-1.        ,  0.        ,  2.15774703],
                [-1.        ,  1.48469925,  2.15774703],
                [ 1.        , -1.45185316,  2.15774703],
                [1.        , 0.        , 2.15774703]],
               [[-0.5       ,  0.5       ,  2.15774703],
                [-0.54151225, -0.01831961,  1.25261807],
                [0.20227045, 0.15929604, 1.22609735],
                [ 0.65001631, -0.07832319,  1.00388753]],
               [[-1.88212955, -0.26279324, -0.45193356],
                [-0.56988615, -0.42863035,  0.45224738],
                [ 0.05436653, -0.26342535,  0.4620285 ],
                [ 0.        ,  0.        , -0.93721968]],
               [[-1.        ,  0.        , -0.93721968],
                [-0.93310571,  0.9774524 , -1.01902032],
                [ 1.        , -1.        , -0.93721968],
                [ 1.        ,  0.        , -0.93721968]]]'''
mesh = bpy.data.meshes.new(name="New Object Mesh")
calc_gregory_surf(kk, kk1, mesh)
obj = bpy.data.objects.new("MyObject", mesh)
bpy.context.collection.objects.link(obj)