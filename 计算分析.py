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