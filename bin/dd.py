#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2021-09-03 13:33:48
# @Author  : Muxiaoxiong 
# @email   : xiongweinie@foxmail.com


def combinations(source: list,n:int)->list:
    '''从一个元素不重复的列表里选出n个元素
    :参数 source:元素不重复的列表
    :参数 n: 要选出的元素数量，正整数，小于等于列表source的长度
    '''
    # 如果n正好等于列表的长度，那么只有一种组合
    # 就是把所有元素都选出来
    if len(source)==n:
        return [source]
    # 如果n是1，那么列表中的每个元素都是一个组合
    if n == 1:
        ans = []
        for i in source:
            ans.append([i])
        return ans
    # 下面处理n小于列表长度的情况
    ans = []
    # 从列表里选n个元素出来，可以理解为先把列表里的第0个元素拿出来放进组合
    # 然后从剩下的元素里选出n-1个
    for each_list in combinations(source[1:],n-1):
        ans.append([source[0]]+each_list)
    # 还可以直接从剩下的元素里选出n个
    for each_list in combinations(source[1:],n):
        ans.append(each_list)
    return ans

print(combinations(["A",'B','C','D','E','F'],3))
