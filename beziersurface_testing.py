import bpy
import bmesh
import mathutils
import math
import numpy as np
from typing import List, Optional, Tuple
from enum import Enum
import time
t = time.time()
TH = 0.00001
TH2 = TH**2
N = 32 #number of subdivisions
N1 = N + 1

def bernstein(i: int, u: float):
    return math.comb(3, i) * u**i * (1-u)**(3-i)

BERNS = np.array([[bernstein(i, k/N) for i in range(4)] for k in range(N1)])
BERNS2 = np.multiply(*np.meshgrid(BERNS, BERNS)).reshape(N1*4, N1, 4).transpose(1, 0, 2).reshape(N1, N1, 4, 4)

def bbel(i,j,nu,nv):
    return BERNS[nu][i] * BERNS[nv][j]

B2 = [[[[bbel(i,j,nu,nv) for i in range(4)] for j in range(4)] for nu in range(N1)] for nv in range(N1)]

def same_coords(c1: mathutils.Vector, c2: mathutils.Vector) -> bool:
    return (c1-c2).length_squared < TH2

def get_loc(name: str):
    return np.array(bpy.data.objects[name].location)

#kk = np.array([[get_loc(f"k{j}{i}")for i in range(4)] for j in range(4)])

kk:np.ndarray = np.array([
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
                np.array([ 1.        ,  0.        , -0.93721968])]])

kk_transposed = np.expand_dims(np.expand_dims(kk.transpose(2,0,1),1),1)

def sumk(nu, nv):
    el = []
    for d in range(3):
        el.append(sum(sum(kk[j][i][d] * B2[nv][nu][j][i] for j in range(4)) for i in range(4)))
    return el

res2 = [[sumk(nu, nv) for nu in range(N1)] for nv in range(N1)]
print(np.array(res2))

print("!!!!!!!!!!!!!")
res = np.sum((BERNS2 * kk_transposed),(3,4)).transpose(2,1,0)
print(res)
#print(time.time() - t)