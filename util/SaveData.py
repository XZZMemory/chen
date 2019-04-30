# -*- coding:utf-8 -*-
from . import RnWLists

def saveIndividual(path, POPUL, flag):
    '''
    将解的信道分配情况和功率分配情况写入TXT， flag为0覆写，flag为1续写
    '''
    if flag == 0:
        RnWLists.clearTxt(path)
    
    for x in POPUL:
        RnWLists.writeTxt(path, x.genec, 1)  #追加写入信道分配
        RnWLists.writeTxt(path, x.genep, 1) #追加写入功率分配