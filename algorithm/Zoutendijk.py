# -*- coding:utf-8 -*-
from scipy.optimize import linprog
import numpy as np

#import threading

#lock = threading.Lock()

'''
输入部分
变量个数，变量上界变量下界 list
目标函数 f(x)------------不需要
当前偏导数list
当前可行解list
约束函数（矩阵）Ax ≥b Ex = e  （A, b, E, e）
----------------注意是"大于等于", "≥", ">=", "大於等於", 所以务必转换成大于等于形式

约束函数（矩阵）Ax ≤b Ex = e  （A, b, E, e）
----------------17.06.28更新放弃了pulp改为scipy之后就是"小于等于", "≤", "<=", "小於等於", 所以务必转换成小于等于形式(最后还是用的大于等于)

'''
class ZoutendijkProblem:
    '''
    问题类
    '''
    def __init__(self, feasibleX, deltaF, matrixA = None, vectorb = None, matrixE = None, vectore = None, xUpBounds = None, xLowBounds = None, sense = 'Minimize'):
        '''
        传入参数, feasibleX 可行解list，deltaF 偏导数list
        matrixA 矩阵A，vectorb 向量b  AX ≥ b
        matrixE 矩阵E，vectore 向量e，Ex = e
        xUpBounds x上界，xLowBounds x下界
        sense 默认Minimize，或者取Maximize
        '''
        self.feasibleX = feasibleX
        self.deltaF = deltaF
        self.matrixA = matrixA
        self.vectorb = vectorb
        self.matrixE = matrixE
        self.vectore = vectore
        self.xUpBounds = xUpBounds
        self.xLowBounds = xLowBounds
        self.sense = sense
        
        self.matrixA1 = []    #拆解成为的矩阵
        self.vectorb1 = []    #和向量
        self.matrixA2 = []    #拆解后
        self.vectorb2 = []    #在给定点A1*x  = b1，A2*x≥b2
        
        
        self.xNum = len(feasibleX)  #x个数
        
        self.vectord = [0 for i in range(self.xNum)]  #线性规划解
        self.dUpBounds = [1 for i in range(self.xNum)]
        self.dLowBounds = [-1 for i in range(self.xNum)]  #d的上下界
        
        self. ifKKT = 1
        self.MaxLambda = float("inf")        
    
    
    def VectorDot(self, a, b):
        '''
        计算内积的小函数
        '''
        n = len(a)
        '''
        if n != len(b):  #检测长度不一致错误
            print("ERROR, vector's length dose not match")
            exit(0)
        '''
        sum = 0
        for i in range(n):
            sum += a[i] * b[i]
            
        return sum
            
    
    def DepartMatrixA(self):
        '''
        只有当Ax≥b条件存在时才能执行-------------------------------切记-
        分解矩阵A，使其成为A1和A2，其中A1是当前已经达到边界条件的参数矩阵，A2是当前未达到边界条件的参数矩阵，同理分解b1和b2
        '''
        if self.matrixA == None or self.vectorb == None:
            print("ERROR, bound Ax≥b is not exit")
            exit(0)
        n = len(self.matrixA)
        if n != len(self.vectorb):  #检测左端的表达式与右端的值数量是否匹配
            print("ERROR, length of Vector b dose not match the lenght of MatrixA")
            exit(0)
            
        for i in range(n):  #分解矩阵和向量
            if abs(self.VectorDot(self.feasibleX, self.matrixA[i]) - self.vectorb[i])<0.0000001:  #判断是否到达边界
                self.matrixA1.append(self.matrixA[i])
                self.vectorb1.append(self.vectorb[i])
                #print(self.matrixA[i])
                #print("加入A1")
            else:
                self.matrixA2.append(self.matrixA[i])
                self.vectorb2.append(self.vectorb[i])
                #print(self.matrixA[i])
                #print("加入A2")
                
    def GetLPBound(self):
        '''
        获取LP问题的上下界，|d|<1，d>=0或者d<=0
        '''
        if self.xUpBounds != None:
            for i in range(self.xNum):
                if self.feasibleX[i] == self.xUpBounds[i]:
                    self.dUpBounds[i] = 0   #如果x已经抵达上界，则对应的可行方向d一定小于0，所以d的上界是0

        if self.xLowBounds != None:
            for i in range(self.xNum):
                if self.feasibleX[i] == self.xLowBounds[i]:
                    self.dLowBounds[i] = 0   #如果x已经抵达下界，则对应的可行方向d一定大于0，所以d的下界是0

    def LinearProgramming(self):
        '''
        解决线性规划问题
        self.deltaF为目标函数参数（线性的）
        dLowBounds为d下界
        dUpBounds为d上界
        上下界由GetLPBound计算好
        
        '''
        
        c = np.array(self.deltaF)
        if self.sense == "Maximize" :
            c = -c
        if len(self.matrixA1) != 0:
            a = np.mat(self.matrixA1)
            a = -a
        else:
            a = None
            
        if len(self.vectorb1) != 0: 
            b = np.array([0 for i in range(len(self.vectorb1))])
            b = -b
        else:
            b = None
        
    
        boundslist = []
        for i in range(self.xNum):  #生成d上下界
            boundslist.append((self.dLowBounds[i], self.dUpBounds[i]))
        bounds = tuple(boundslist)   #d上下界转元组
        #print(bounds)
        
        #lock.acquire()  #线程锁似乎没什么用
        res = linprog(c, A_ub = a, b_ub = b, bounds = bounds)  #计算线性规划
        #lock.release()
        #print(res)
        vectord = res['x']
        for i in range(self.xNum):   #numpy array 转list
            self.vectord[i] =round(vectord[i],13)

    def GetLambda(self):
        '''
        计算最大步长lambda max
        '''
        n = len(self.matrixA2)
        MaxLambda = float("inf")
        for i in range(n):
            bhat = self.vectorb2[i] - self.VectorDot(self.matrixA2[i], self.feasibleX)
            dhat = self.VectorDot(self.matrixA2[i], self.vectord)
            if dhat < 0 and bhat/dhat<MaxLambda:  #取最小
                MaxLambda = bhat/dhat
        
        #print("位置1处MaxLambda = %f" %(MaxLambda))
                
        if self.xLowBounds != None:  #x下界条件存在
            for i in range(self.xNum):
                bhat = self.xLowBounds[i]-self.feasibleX[i]
                if abs(self.vectord[i])<0.00001 or abs(bhat)<0.00001:
                    self.vectord[i] = 0
                dhat = self.vectord[i]
                if bhat != 0 and dhat < 0 and bhat/dhat<MaxLambda:
                    MaxLambda = bhat/dhat
        #print("位置2处MaxLambda = %f" %(MaxLambda))

        if self.xUpBounds != None:  #x上界条件
            for i in range(self.xNum):
                bhat = self.feasibleX[i]-self.xUpBounds[i]
                if abs(self.vectord[i])<0.00001 or abs(bhat)<0.00001:
                    self.vectord[i] = 0
                dhat = -self.vectord[i]
                if bhat != 0 and dhat < 0 and bhat/dhat<MaxLambda:
                    MaxLambda = bhat/dhat
                    
        #print("位置3处MaxLambda = %f" %(MaxLambda))
      
        self.MaxLambda = round(MaxLambda, 13)
        
        
    def solve(self):
        '''
        对问题求解
        '''
        if self.matrixA != None and self.vectorb != None:   #如果不等式条件存在，分解矩阵
            self.DepartMatrixA()
            
        self.GetLPBound()   #获取d的上下限
        self.LinearProgramming()  #线性规划
        self.GetLambda()  #获取lambda
        return (self.MaxLambda, self.vectord)
    
    def resolve(self):
        '''
        需要修改为多条件限制，以后再说  2017.06.07
        '''
        if self.MaxLambda == float("inf"):
            print("strat to resolve")
            for i in range(len(self.matrixA)):
                print("resolve times ---------------" + str(i))
                switchA = self.matrixA.pop(0)
                switchb = self.vectorb.pop(0)
                self.matrixE = [] if self.matrixE == None else self.matrixE
                self.vectore = [] if self.vectore == None else self.vectore
                self.matrixE.append(switchA)
                self.vectore.append(switchb)
                NewSolve = self.solve()
                #print(NewSolve)
                if NewSolve[0] != float("inf"):
                    print("over")
                    break
                switchA = self.matrixE.pop(-1)
                switchb = self.vectore.pop(-1)
                self.matrixA.append(switchA)
                self.vectorb.append(switchb)            
            
    def BoundarySolve(self):
        '''
        如果找不到解就改变条件
        AX>b的条件加入EX = e中，尝试求出不同的线性规划解
        '''
        self.solve()
        if self.MaxLambda == float("inf"):
            print("strat to resolve")
            for i in range(len(self.matrixA)):
                print("resolve times ---------------" + str(i+1))
                switchA = self.matrixA.pop(0)
                switchb = self.vectorb.pop(0)
                self.matrixE = [] if self.matrixE == None else self.matrixE
                self.vectore = [] if self.vectore == None else self.vectore
                self.matrixE.append(switchA)
                self.vectore.append(switchb)
                NewSolve = self.solve()
                #print(NewSolve)
                if NewSolve[0] != float("inf"):
                    break
                switchA = self.matrixE.pop(-1)
                switchb = self.vectore.pop(-1)
                self.matrixA.append(switchA)
                self.vectorb.append(switchb)
        #else:
            #problemSolve = self.solve()
            
        return (self.MaxLambda, self.vectord)
    
def test():
    '''
    测试1
    FeasibleX = [0, 0]
    DeltaF = [-2, -4]
    MatrixA = [[-2, 1], [-1, -1]]
    vectorB = [-1, -2]
    xLowBound = [0, 0]
    
    测试2
    FeasibleX = [0, 2]
    DeltaF = [0, 16]
    MatrixA = [[1, 1], [15, 10]]
    vectorB = [1, 12]
    xLowBound = [0, 0]
    
    测试3
    FeasibleX = [0, 1.2]
    DeltaF = [0, 9.6]
    MatrixA = [[1, 1], [15, 10]]
    vectorB = [1, 12]
    xLowBound = [0, 0]
    
    测试4
    FeasibleX = [0.4 , 0.6]
    DeltaF = [4/5, 24/5]
    MatrixA = [[1, 1], [15, 10]]
    vectorB = [1, 12]
    xLowBound = [0, 0]
    
    测试5
    FeasibleX = [0.8 , 0.2]
    DeltaF = [8/5, 8/5]
    MatrixA = [[1, 1], [15, 10]]
    vectorB = [1, 12]
    xLowBound = [0, 0]
    '''
    FeasibleX = [0, 2]
    DeltaF = [0, 16]
    MatrixA = [[1, 1], [15, 10]]
    vectorB = [1, 12]
    xLowBound = [0, 0]
    z = ZoutendijkProblem(feasibleX=FeasibleX, deltaF=DeltaF, matrixA=MatrixA, vectorb=vectorB, xLowBounds=xLowBound)
    print(z.solve())


    
#test()


            
    

    
    