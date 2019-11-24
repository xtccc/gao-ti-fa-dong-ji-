#!usr/bin/env python
#-*- coding:utf-8 _*-
"""
@author:tianchang
@file: 计算分析.py
@time: 2019/11/23 19:18-10s
"""
import gaoti_tanli
import cProfile, pstats, io
from pstats import SortKey
pr = cProfile.Profile()
pr.enable()


# ... do something ...
gaoti_tanli.solve(20,20,5,0.6949,0.7515)
#108s

pr.disable()
s = io.StringIO()
sortby = SortKey.CUMULATIVE
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print(s.getvalue())


'''

74985 function calls in 108.015 seconds
   Ordered by: cumulative time
   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1   25.538   25.538  108.015  108.015 C:\Users\tianchang\膏体\gaoti_tanli.py:5(solve)
     2499    0.011    0.000   82.477    0.033 <__array_function__ internals>:2(solve)
     2499    0.021    0.000   82.464    0.033 {built-in method numpy.core._multiarray_umath.implement_array_function}
     2499   82.321    0.033   82.443    0.033 C:\Users\tianchang\AppData\Local\Programs\Python\Python38\lib\site-packages\numpy\linalg\linalg.py:327(solve)

'''