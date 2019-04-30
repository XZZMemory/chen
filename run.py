# -*- coding:utf-8 -*-
import multiprocessing
import time

import _create_case
import solve
from util import RnWLists


def GetDM(Sets):
    '''
    获取距离矩阵
    输入的是坐标，一共4行，分别为基站X，基站Y，用户X，用户Y
    '''
    DM = []  #要返回的DM矩阵
    li = []
    Piconum = len(Sets[0])
    Usernum = len(Sets[2])
    
    for i in range(Usernum): #计算跟macro基站的距离
        li.append(_create_case.CulculateDistance(0, 0, Sets[2][i], Sets[3][i]))
    
    DM.append(li)
    for i  in range(Piconum):
        li = []
        for j in range(Usernum):
            li.append(_create_case.CulculateDistance(Sets[0][i], Sets[1][i], Sets[2][j], Sets[3][j]))
        
        DM.append(li)
    
    return DM
    


macroNumber = 1 #macro BS数量
infofmbs = [64,500,20]  #macro基站信息[信道数，覆盖范围（米），功率（瓦）]
infofpbs = [64,100,1]  #pico基站信息[信道数，覆盖范围（米），功率（瓦）]
maxchannal = 3   #每个用户可用的最大信道数量
nmub = 2   #每个用户可连接的最大基站数量
Qmax = 10000    #功率等级
B = 150       #信道带宽
Alpha = 10  #接收功率阈值，即接收功率大于等于alpha时传输成功  (????)
m = 100   #种群个数
g = 1000  #迭代次数
pm  = 0.05 #变异率

def cutData(data, longth):
    '''
    截取部分用户个体,data为坐标列表，longth为截取的用户数量
    '''
    userX = data[2][:longth]  #截取部分用户坐标，从第一个到第length个（大概）
    userY = data[3][:longth]
    newData = [data[0], data[1], userX, userY]  #[基站坐标x，基站坐标y，截取的部分用户坐标x，截取的部分用户坐标y]
    
    return newData

def MakeParamas(m, g, CaseRootPath, N, flag = 1, CutNumber = 60, PicoNum = 3, UserNum = 70,  Times = 8):

    begin = time.clock()
    ParamasList = [] #参数列表
    Class = 'R'
    CasePath = create_txt(CaseRootPath, PicoNum, UserNum, Class, N)
    for i in range(Times):
        data = RnWLists.readTxt(CasePath, 1)  #读取数据
        if flag == 0:
            data[0][0] = []  #这两个是用来清空pico坐标，形成对照组的
            data[0][1] = []
            CaseClass = 'O'
            
        Stes = cutData(data[0], CutNumber)  #切取部分用户数
        DM = GetDM(Stes)   #计算DM矩阵
        
        BSList = solve.SetBS(1, PicoNum, infofmbs, infofpbs, DM)
        
        Parama = [m, BSList, g,  CutNumber]
        ParamasList.append(Parama)  #加入参数列表
    
    end = time.clock()
    use = int(end-begin)
    minutes = int(use/60)
    secends = use%60
    # print("构造参数列表用时%d分%d秒" %(minutes, secends))
        
    return ParamasList

def solve_ga(Paramas):

    m = Paramas[0]
    BSList = Paramas[1]
    g = Paramas[2]
    userNumber = Paramas[3]
    Results = []
    begin = 25  #速率最小值
    end = 275  #速率最大值
    for LowestRate in range(begin, end, 25):
        Result = solve.ga(m, BSList, g, userNumber, LowestRate =LowestRate / 100)
        Results.append(Result)
    return Results

def run_ga(m, g, CaseRootPath, ResultRootPath, flag=1, ProcessNumber = 8, begin = 20, end = 60, PicoNumber = 3, CaseNumber = 0, times = 8):

    TotalTimerBegin = time.clock()   
    
    if end<begin:
        print('错误，end不能小于begin')
        exit(0)
        
    if flag == 0:
        Clas = 'O'
    elif flag == 1 :
        Clas = 'R'
    elif flag == 2:
        Clas = 'Z'
    else:
        Clas = 'N'
        
    pool = multiprocessing.Pool(ProcessNumber)
    
    for i in range(begin, end+1, 10):  #按照用户数量循环计算,以及间隔数量
        TimerBegin = time.clock()
        print("本次循环中有%d个用户" %(i))
        Paramas = MakeParamas(m, g, CaseRootPath, CaseNumber, flag=flag, CutNumber=i, PicoNum = PicoNumber, Times = times)
        CasePath = create_txt(CaseRootPath, PicoNum=PicoNumber, UserNum = i, Class = Clas, N = CaseNumber)
        
        Result = pool.map(solve_ga, Paramas) #并行计算
        
        print("result 长度为")
        print(len(Result))
        
        SavePath = create_txt(ResultRootPath, PicoNumber, i, Clas, CaseNumber) #生成保存位置
        RnWLists.writeTxt(SavePath, Result, 1)  #循环存入数据
        print("写入数据中，%d次重复数据" %(times))
        print("当前样例为,  %s"%(CasePath))
        print("写入位置为, %s"%(SavePath))

        
        
        TimerEnd = time.clock()
        use = int(TimerEnd-TimerBegin)
        minutes = int(use/60)
        secends = use%60
        print("%d个用户，用时%d分%d秒" %(i, minutes, secends))
        


    #SavePath = create_txt(ResultRootPath, 'save', 3, 60, Clas, CaseNumber)
    
    #RnWLists.writeTxt(SavePath, result, 1)  #将结果写入
    
    TotalTimerEnd = time.clock()
    use = int(TotalTimerEnd-TotalTimerBegin)
    minutes = int(use/60)
    secends = use%60
    print("当前计算结束，用时%d分%d秒，结果写入%s" %(minutes, secends, SavePath))

def create_txt(path, PicoNum, UserNum, Class, N):
    '''
    建立TXT文件, 路径，pico基站数int，用户数int，类别string, 编号int
    '''
    result = path + '/'+str(PicoNum)+'_'+str(UserNum)+'_'+Class+'_'+str(N)+'.txt'
    return result

def ga_convergence(m, g, CasePath, ConvergenceSavePath):

    data = RnWLists.readTxt(CasePath, 1)
    DM = GetDM(data[0])
    
    picoNumber = len(DM)-1
    userNumber = len(DM[0])
    
    begin = time.clock()
    BSList = solve.SetBS(macroNumber, picoNumber, infofmbs, infofpbs, DM)

    
    print("执行遗传算法初始化")
    Result = solve.ga_convergence(m = m, BSList = BSList, g = g, userNumber = userNumber)
    RnWLists.writeTxt(ConvergenceSavePath, Result, 1)
    print("执行遗传算法完成")


    end = time.clock()
    use = int(end-begin)
    minutes = int(use/60)
    secends = use%60
    print("程序执行用时%d分%d秒" %(minutes, secends))

flag = 1

path_1 = './data/case/6_70_R_6.txt'
path_2 = './data/case/6_40_R_2.txt'
CasePath = path_1
ConvergenceSavePath = './data/convergence/conv.txt'
RootPath = './data'
CaseRootPath = './data/case'
ResultRootPath = './data/result/data'


if __name__ == '__main__' :
    for i in range(50):
        ga_convergence(m, g, CasePath, ConvergenceSavePath) 

    # for j in range(4):
    #     run_ga(m, g, CaseRootPath, ResultRootPath, flag=1, ProcessNumber=10, begin=40, end=70, PicoNumber=j+3, CaseNumber=6, times=10)
    #run_ga(m, g, CaseRootPath, ResultRootPath, flag=1, ProcessNumber=5, begin=40, end=40, PicoNumber=6,CaseNumber=2, times=5)