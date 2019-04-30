# -*- coding:utf-8 -*-
'''
Created on 2017年3月8日

@author: Graaf.S.Angell
'''

'''

用来生成测试用例的工具函数，包括：
1.花式随机生成pico基站坐标位置
2.花式随机生成用户坐标位置
3.生成的位置画出来

'''


import numpy as np
import math
import random
from math import sqrt
import matplotlib.pyplot as plt

pi, sin, cos = np.pi, np.sin, np.cos

def CulculateDistance(x1, y1,x2,y2):
    '''
    计算两点距离的小函数
    '''
    d =sqrt((x1 - x2)**2 + (y1 - y2)**2)
    return d

def LayPicoBS(MacroRange, PicoRange, num, flag):
    '''
    MacroRange, PicoRange分别为宏基站和微基站的覆盖半径
    设置微基站的位置，num是数量，flag是标记，0为边缘放置，1为逐渐远离中心放置，2为随机放置，并且尽量避免重叠
    '''
    PicoSetX = []
    PicoSetY = []  #要返回的坐标X,和Y的数组
    #用极坐标设定比较方便
    AverageTheta = pi*2/num  #求出平均间隔的角度
    length = 0  #设置距离初始为0   
    if flag == 0:  #边缘放置
        theta = 0   #极角的角度从0开始转一圈
        length = MacroRange - PicoRange #极径长度为最远
        for i in range(num):  #循环pico基站个数次
            PicoSetX.append(sin(theta)*length)  #算横坐标
            PicoSetY.append(cos(theta)*length) #算纵坐标，加入纵坐标数组
            theta+=AverageTheta  #计算过一次之后极角增加
        return (PicoSetX, PicoSetY)  #返回的元组，包含两个数组

    if flag == 1:    #按照距离排列
        length = 1*PicoRange
        AverageLength = (MacroRange-2*PicoRange)/(num+1) #极径的平均增加长度
        for i in range(num+1):
            length += AverageLength  #极径逐渐增加
            theta = pi*(math.log(i+1)*0.9+0.09*i)  #这个挺复杂的，就是为了尽量减少重叠区域
            PicoSetX.append(sin(theta)*length)  #算横坐标
            PicoSetY.append(cos(theta)*length) #算纵坐标，加入纵坐标数组
            
        PicoSetX.pop(0) #第一个的重复面积最大，删掉
        PicoSetY.pop(0)
        return (PicoSetX, PicoSetY)  #返回的元组，包含两个数组
    
    if flag == 2: #随机生成，并且尽量避免重叠
        length = random.random()*(MacroRange-2*PicoRange)+PicoRange  #随机生成极径
        theta = random.random()*2*pi #随机生成极角
        PicoSetX.append(sin(theta)*length)  #算横坐标
        PicoSetY.append(cos(theta)*length) #算纵坐标，加入纵坐标数组
        laid = 0
        while(laid<num-1):
            length = random.random()*(MacroRange-2*PicoRange)+PicoRange  #随机生成极径
            theta = random.random()*2*pi #随机生成极角
            SetX = sin(theta)*length
            SetY = cos(theta)*length
            IfTooNear = 0  #是否太近
            for i in range(len(PicoSetX)):
                if CulculateDistance(SetX, SetY, PicoSetX[i], PicoSetY[i])<PicoRange*2:  #如果有重叠
                    IfTooNear = 1
                    break
            if IfTooNear == 0:  #没有重叠才加入
                laid+=1
                PicoSetX.append(SetX)
                PicoSetY.append(SetY)
        
        return(PicoSetX, PicoSetY)

def LayUsers(MacroRange, PicoRange, num, flag, BSSet, pro = 0.75):
    '''
    MacroRange宏基站范围, PicoRange微基站范围, num放置的用户数量, BSSet为基站坐标
    flag标记为0均匀分布，标记为1，一定比率集中于微基站范围
    
    用户分配策略，如果flag为0，均匀分布，根据覆盖面积比计算用户是否在pico基站中
    如果flag为1，则按照比率决定用户在pico基站中或者pico基站外
    在pico基站中则轮流放置
    '''
    UserSetX = []
    UserSetY = []
    laid = 0
    
    if flag == 0:
        probability = PicoRange**2*len(BSSet)/(MacroRange**2)#计算面积比
    else:
        probability = pro
    
    point = 0  #标记
    while(laid<num):
            if random.random()<=probability:  #意味着这个用户要分配进pico基站中
                length = random.random()*PicoRange*1  #随机生成极径
                theta = random.random()*2*pi #随机生成极角
                SetX = sin(theta)*length + BSSet[0][point] #算横坐标
                SetY = cos(theta)*length + BSSet[1][point] #算纵坐标(在第point个pico基站范围内)   
                point = point+1 if point<len(BSSet) else 0  #循环标记在哪个pico基站里出现
                UserSetX.append(SetX)
                UserSetY.append(SetY)
                laid+=1
            else:  #放置在非pico基站范围
                IfInPico = 1  #标记，是否在pico基站范围内                
                while IfInPico == 1:  #如果在某个pico基站范围内，则标记为1，否则为0
                    IfInPico = 0
                    length = random.random()*MacroRange*1  #随机生成极径
                    theta = random.random()*2*pi #随机生成极角
                    SetX = sin(theta)*length  #算横坐标
                    SetY = cos(theta)*length  #算纵坐标(在第point个pico基站范围内)   
                    for i in range(len(BSSet[0])):  #循环pico基站个数次
                        if CulculateDistance(BSSet[0][i], BSSet[1][i], SetX, SetY)<PicoRange:  #如果用户位置在pico基站范围内
                            IfInPico = 1  #标记置1
                UserSetX.append(SetX)
                UserSetY.append(SetY)
                laid+=1
    return (UserSetX, UserSetY)
def LayUsersReal(MacroRange, PicoRange, num, flag, BSSet, pro = 0.75):
    '''
    此为简单随机
    
    MacroRange宏基站范围, PicoRange微基站范围, num放置的用户数量, BSSet为基站坐标
    flag标记为0均匀分布，标记为1，一定比率集中于微基站范围
    
    用户分配策略，如果flag为0，均匀分布，根据覆盖面积比计算用户是否在pico基站中
    如果flag为1，则按照比率决定用户在pico基站中或者pico基站外
    在pico基站中则轮流放置
    '''                                               
    UserSetX = []
    UserSetY = []
    laid = 0
    
    while(laid<num):                    
            length = random.random()*MacroRange*1  #随机生成极径
            theta = random.random()*2*pi #随机生成极角                    
            SetX = sin(theta)*length #算横坐标
            SetY = cos(theta)*length #算纵坐标             
            if flag == 0:
                laid+=1
                UserSetX.append(SetX)  #算横坐标
                UserSetY.append(SetY) #算纵坐标，加入纵坐标数组
            else:
                IfInPico = 0  #标记，是否在pico基站范围内
                for i in range(len(BSSet[0])):  #循环pico基站个数次
                    if CulculateDistance(BSSet[0][i], BSSet[1][i], SetX, SetY)<PicoRange:  #如果用户位置在pico基站范围内
                        IfInPico = 1  #标记置1
                        break
                if ((IfInPico == 0)and(random.random()>pro))or((IfInPico == 1)and(random.random()>(1-pro))):  
                    #如果在PICO基站范围内则有75%几率保留，不在pico基站范围内则有25%几率保留
                    laid += 1
                    UserSetX.append(SetX)  
                    UserSetY.append(SetY) #加入纵坐标数组                    
    
    return(UserSetX,UserSetY)


def DrawCircle(radius, xShifting = 0, yShifting = 0, color='r'):
    '''
    画一个圆圈，radius半径 xShifting，yShifting圆心偏移量，color颜色
    不要问我为什么，就这么弄就能画出来
    '''
    theta = np.linspace(0, 2*pi, 360)
    x = sin(theta)*radius+xShifting
    y = cos(theta)*radius + yShifting

    plt.plot(x,y,color)

def DrawCase(macrorange, picorange, picotum, usertum):
    '''
    画出图像
    '''
    DrawCircle(macrorange,0,0,'b')
    
    for i in range(len(picotum[0])):
        DrawCircle(picorange, picotum[0][i], picotum[1][i], 'b')
        
    for i in range(len(usertum[0])):
        plt.plot(usertum[0][i], usertum[1][i], marker = 'o')
        
    plt.show()


    