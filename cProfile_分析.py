#!usr/bin/env python
#-*- coding:utf-8 _*-
"""
@author:tianchang
@file: cProfile_分析.py
@time: 2019/11/23 19:02-02s
"""

import cProfile
import gaoti

import  numpy as np
N_x_grid=100*100
N_y_grid=100*100
T = np.zeros((N_x_grid,N_y_grid))

result= np.zeros(N_x_grid*N_y_grid)
for i in range(N_x_grid*N_y_grid):
        result[i] = N_x_grid**5+12354
def aa(N_x_grid,N_y_grid,result):
    for i in range (0, N_x_grid):
        for j in range (0, N_y_grid):
            T[i][j] = result[i * N_y_grid + j]
    return T



cProfile.run('T = result.reshape (N_x_grid, N_y_grid)')
cProfile.run('T = aa(N_x_grid,N_y_grid,result)')


'''
   Ordered by: standard name
   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.011    0.011    0.011    0.011 <string>:1(<module>)
        1    0.001    0.001    0.012    0.012 {built-in method builtins.exec}
        1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
        1    0.000    0.000    0.000    0.000 {method 'reshape' of 'numpy.ndarray' objects}
         4 function calls in 69.435 seconds
   Ordered by: standard name
   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.000    0.000   69.435   69.435 <string>:1(<module>)
        1   69.435   69.435   69.435   69.435 cProfile_分析.py:20(aa)
        1    0.000    0.000   69.435   69.435 {built-in method builtins.exec}
        1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
'''