#!usr/bin/env python
#-*- coding:utf-8 _*-
"""
@author:tianchang
@file: cProfile_分析.py
@time: 2019/11/23 19:02-02s
"""

import cProfile
import gaoti_tanli

cProfile.run('gaoti_tanli.solve(20,20,5,0.6949,0.7515)')
'''
从分析报告结果中我们可以得到很多信息：

    整个过程一共有197个函数调用被监控，其中192个是原生调用（即不涉及递归调用）总共执行的时间为0.002秒结果列表中是按照标准名称进行排序，也就是按照字符串的打印方式（数字也当作字符串）在列表中：
        ncalls表示函数调用的次数（有两个数值表示有递归调用，总调用次数/原生调用次数）tottime是函数内部调用时间（不包括他自己调用的其他函数的时间）percall等于 tottime/ncallscumtime累积调用时间，与tottime相反，它包含了自己内部调用函数的时间最后一列，文件名，行号，函数名
'''