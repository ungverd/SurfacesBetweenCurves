import bpy
import bmesh
import mathutils
import math
from typing import List
# from enum import Enum
import time

N = 64 # number of subdivisions
N1 = N + 1

TH = 0.00001
TH2 = TH**2

def bernstein(i: int, u: float):
    return math.comb(3, i) * u**i * (1-u)**(3-i)

BERNS = [[bernstein(i, k/N) for i in range(4)] for k in range(N1)]

def calculate_faces():
    faces: List[List[int]] = []
    for i in range(N):
        for j in range(N):
            face = [i * N1 + j]
            face.append(i * N1 + j + 1)
            face.append((i + 1) * N1 + j + 1)
            face.append((i + 1) * N1 + j)
            faces.append(face)
    return faces

FACES = calculate_faces()

try:
    #raise(ModuleNotFoundError)
    import numpy as np
    import numpy.typing as npt
    t = time.time()
    BERNS_NP = np.array(BERNS)
    BERNS2 = np.multiply(*np.meshgrid(BERNS_NP, BERNS_NP)).reshape(N1*4, N1, 4).transpose(1, 0, 2).reshape(N1, N1, 4, 4)
    # numpy magic to get products of all combinations of Bernstein coefficients

    def calc_bezier_surf(control_points: List[List[mathutils.Vector]], mesh: bpy.types.Mesh):
        cp = np.array(control_points)
        cp_transposed = np.expand_dims(np.expand_dims(cp.transpose(2, 0, 1), 1), 1)
        # numpy magic to be able to multiply control point vectors with Bernstein coefficients
        res: npt.NDArray[np.float32] = np.sum((BERNS2 * cp_transposed), (3, 4)).transpose(2, 1, 0)
        # numpy representation of the formula p(u,v) = sum_i_from_0_to_3(sum_j_from_0_to_3( k(i,j)*B(i,u)*B(j,v) ))
        res = res.reshape(N1 * N1, 3)
        mesh.from_pydata(res, [], FACES)

except (ModuleNotFoundError, ImportError):
    t = time.time()
    def berns2_el(i: int, j: int, nu: int, nv:int):
        return BERNS[nu][i] * BERNS[nv][j]

    B2 = [[[[berns2_el(i,j,nu,nv) for i in range(4)] for j in range(4)] for nu in range(N1)] for nv in range(N1)]

    def calc_point(nu: int, nv: int, control_points: List[List[mathutils.Vector]]):
        el: List[float] = []
        for d in range(3):
            el.append(sum(sum(control_points[j][i][d] * B2[nv][nu][j][i] for j in range(4)) for i in range(4)))
        return mathutils.Vector(el)

    def calc_bezier_surf(control_points: List[List[mathutils.Vector]], mesh: bpy.types.Mesh):
        res = [calc_point(nu, nv, control_points) for nu in range(N1) for nv in range(N1)]
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
                          curve_b3: "BezierCurve"):
        a0 = self.get_a0(curve_a0)
        a3 = self.get_a3(curve_a3)
        b0 = self.get_b0(curve_b0)
        b3 = self.get_b3(curve_b3)
        c0 = self.c0
        c1 = self.c1
        c2 = self.c2


def get_loc(name: str) -> mathutils.Vector:
    return bpy.data.objects[name].location


def same_coords(c1: mathutils.Vector, c2: mathutils.Vector) -> bool:
    return (c1-c2).length_squared < TH2

kk = [[get_loc(f"k{j}{i}")for i in range(4)] for j in range(4)]

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
calc_bezier_surf(kk, mesh)
obj = bpy.data.objects.new("MyObject", mesh)
bpy.context.collection.objects.link(obj)