# -*- coding:utf-8 -*-
from _create_case import LayPicoBS, LayUsers, LayUsersReal
from util import RnWLists


MacroRange = 500
PicoRange = 100
PicoNum = 3
UserNum = 40
flag = 0



def CreatTxt(path, PicoNum, UserNum, Class, N):
    '''
    建立TXT文件, 路径，pico基站数int，用户数int，类别string, 编号int
    '''
    result = path + '/'+str(PicoNum)+'_'+str(UserNum)+'_'+Class+'_'+str(N)+'.txt'
    return result

def CreatCase(path, PicoNum, UserNum, Num):
    '''
    批量创建测试用例, 路径，pico基站数，用户数，创建个数
    如此创建的结果，用户分布按照密集型分布，基站分布分为两种，边缘或处于用户密集处
    所以两种基站分布方式场景用户的位置是一样的，保证计算结果在同一场景下进行
    '''
    
    for i in range(Num):
        PathR = CreatTxt(path, PicoNum, UserNum, 'R', i)  #边缘分布用例
        PathN = CreatTxt(path, PicoNum, UserNum, 'N', i) #集中分布用例
        print(PathR)
        
        tun1 = LayPicoBS(MacroRange, PicoRange, PicoNum, 2)  #集中分布基站位置
        tun2 = LayPicoBS(MacroRange, PicoRange, PicoNum, 0)  #边缘分布挤占位置
        tun3 = LayUsers(MacroRange, PicoRange, UserNum, 1, tun1, 0.5) #集中分布用户位置

        DataR = [tun2[0], tun2[1], tun3[0], tun3[1]]  #边缘分布数据
        DataN = [tun1[0], tun1[1], tun3[0], tun3[1]] #集中分布数据
        
        RnWLists.writeTxt(PathR, DataR, 1)
        RnWLists.writeTxt(PathN, DataN, 1)

def CreatCase2(path, PicoNum, UserNum, Num):
    '''
    新的创建样例方法，根据同样的用户分配在边缘分布3-6个基站，不用集中分布了——————————————17.11.05
    
    批量创建测试用例, 路径，pico基站数，用户数，创建个数
    如此创建的结果，用户分布按照密集型分布，基站分布分为两种，边缘或处于用户密集处
    所以两种基站分布方式场景用户的位置是一样的，保证计算结果在同一场景下进行
    '''
    
    for i in range(Num):

        
        tun1 = LayPicoBS(MacroRange, PicoRange, PicoNum, 2)# 没什么用，为了给tun3填参数
        
        tun3 = LayUsers(MacroRange, PicoRange, UserNum, 0, tun1, 0.5) #均匀分布用户位置
        tun3 = LayUsersReal(MacroRange, PicoRange, UserNum, 0, tun1, 0.5)
        
        for j in range(7):  #循环4次，分别创建基站数为3，4，5，6个的边缘分布样例
            # PicoNumNow = PicoNum + j #当前的pico基站数量
            PicoNumNow = j + 1
            PathR = CreatTxt(path, PicoNumNow, UserNum, 'R', i)  #创建边缘分布用例名称
            print(PathR)
              
            tun2 = LayPicoBS(MacroRange, PicoRange, PicoNumNow, 0)  #边缘分布基站位置
            DataR = [tun2[0], tun2[1], tun3[0], tun3[1]]  #边缘分布数据
        
            RnWLists.writeTxt(PathR, DataR, 1) #写入随机分布位置



path = './data/case'

# CreatCase(path, PicoNum, UserNum, 16)

CreatCase2(path, PicoNum, UserNum, 10)

'''
tum1 = CreatTestCase.LayPicoBS(MacroRange, PicoRange, PicoNum, 2)  #返回pico基站位置元组

tum2 = CreatTestCase.LayUsers(MacroRange, PicoRange, UserNum, 1, tum1, 0.5) #返回用户位置元组

CreatTestCase.DrawCase(MacroRange, PicoRange, tum1, tum2)  #绘制位置

print("如果确认保存输入ok，输入其他不保存")
getinput = input()  #获取输入

if getinput == 'ok' :
    
    data = []
    data.append(tum1[0])
    data.append(tum1[1])
    data.append(tum2[0])
    data.append(tum2[1])
    
    path = r'E:\nsgaii\3_20_N.txt'
    
    RnWLists.writeTxt(path, data, 1)
    print('数据已经保存')
else :
    print('数据被丢弃')
'''