import numpy as np
import numba


n = 124 
A = np.random.rand(n,n,n)
B = np.random.rand(n,n,n)
# @numba.njit
def add(A,B,n):
    for i in range(n):
        for k in range(n):
            for j in range(n):
               A[j,k,i] = B[i,j,k]**2 / 10 + B[i,j,k]**(-2) + B[i,j,k]

import time

start = time.time()
for i in range(10):
    add(A,B,n)
print(time.time()-start)