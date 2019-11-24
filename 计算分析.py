#!usr/bin/env python
#-*- coding:utf-8 _*-
"""
@author:tianchang
@file: 计算分析.py
@time: 2019/11/23 19:18-10s
"""
import gaoti
import cProfile, pstats, io
from pstats import SortKey
pr = cProfile.Profile()
pr.enable()


# ... do something ...
gaoti.solve()

pr.disable()
s = io.StringIO()
sortby = SortKey.CUMULATIVE
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print(s.getvalue())