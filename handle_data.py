# -*- coding:utf-8 -*-
import NSGA_II.RnWLists
import matplotlib

def CreatTxt(path, PicoNum, UserNum, Class, N):
    '''
    建立TXT文件, 路径，pico基站数int，用户数int，类别string, 编号int
    '''
    result = path + '\\'+str(PicoNum)+'_'+str(UserNum)+'_'+Class+'_'+str(N)+'.txt'
    return result




def getFirst2Line(path, path2):  
    '''
    读取前两行写入
    '''
    a = NSGA_II.RnWLists.readTxt(path, 0)
    print(len(a))
    b = [[],[]]
    for i in range(len(a)):
        b[0] += a[i][0]   #读取前两行，速率和功率
        b[1] += a[i][1]
        
    print(len(b[0]))
    c = [b[0], b[1]]
    
    NSGA_II.RnWLists.writeTxt(path2, c, 1)
    

def HandleFirst2Line(path,path2,n):
    '''
    筛选优秀的解写入 ，n为个数
    '''
    a = NSGA_II.RnWLists.readTxt(path, 0)

    b = [[],[]]
    for i in range(len(a)):
        b[0] += a[i][0]   #读取前两行，速率和功率
        b[1] += a[i][1]
        
    c = [[],[]]
    
    for i in range(n):
        c[0].append(b[0][i])
        c[1].append(b[1][i])
    
    for i in range(1,len(b[0])):
        for j in range(len(c[0])):
            if b[0][i]>c[0][j] and b[1][i]<c[1][j]:  #判断支配
                c[0][j]=b[0][i]
                c[1][j]=b[1][i]
                break
        if len(c[0])>n:
            c[0].pop
            c[1].pop
            
        
    d = [c[0], c[1]]
    
    NSGA_II.RnWLists.writeTxt(path2, d, 1)

def HandleFirst2LineBest(path,path2):
    '''
    筛选第一支配层，不限个数
    '''
    a = NSGA_II.RnWLists.readTxt(path, 0)

    b = [[],[]]
    for i in range(len(a)):
        b[0] += a[i][0]   #读取前两行，速率和功率
        b[1] += a[i][1]
        
    c = [[],[]]
    
    for i in range(1,len(b[0])):
        for j in range(len(c[0])):
            if b[0][i]>c[0][j] and b[1][i]<c[1][j]:  #判断支配
                c[0][j]=b[0][i]
                c[1][j]=b[1][i]
                break
            
        
    d = [c[0], c[1]]
    
    NSGA_II.RnWLists.writeTxt(path2, d, 1)

def HandleOne(path):
    '''
    给一个文件算平均数
    '''
    a = NSGA_II.RnWLists.readTxt(path, 0)
    b = [[],[]]
    for i in range(len(a)):
        b[0] += a[i][0]   #读取前两行，速率和功率
        b[1] += a[i][1]
    c = [0,0]
    c[0] = sum(b[0])/len(b[0])
    c[1] = sum(b[1])/len(b[1])
    return c

def HandelMany(RootPath, path2, PicoNum, Begin, End, Class, N):
    '''
    计算很多文件
    '''
    a = [[], []]
    for i in range(Begin, End+1):
        path = CreatTxt(RootPath, PicoNum, i, Class, N)
        b = HandleOne(path)
        a[0].append(b[0])
        a[1].append(b[1])
    
    NSGA_II.RnWLists.writeTxt(path2, a, 1)
        

def HandelMuchMany(RootPath, path2, PicoNum, Begin, End, Class, BeginN = 7, EndN = 15):
    '''
    计算不同样例
    '''
    
    c = [[],[]]
    n = EndN-BeginN+1
    for i in range(Begin, End+1,5):
        a = [[], []]
        for j in range(BeginN, EndN+1):
            path = CreatTxt(RootPath, PicoNum, i, Class, j)
            b = HandleOne(path)
            a[0].append(b[0])
            a[1].append(b[1])
        c[0].append(sum(a[0])/n)
        c[1].append(sum(a[1])/n)
    
    NSGA_II.RnWLists.writeTxt(path2, c, 1)
    
path = r'E:\nsgaii\Results\SingleTarget\2\6_60_R_17.txt'    #读取数据位置
path2 = r'E:\nsgaii\Results\contrastNSGA.txt'                              #写入数据位置
RootPath = r'E:\nsgaii\Results\OandR\6_60_R_7'  

#HandleFirst2Line(path, path2, 40)
getFirst2Line(path, path2)
#HandelMany(RootPath, path2, PicoNum = 3, Begin = 20, End = 60, Class = 'O', N = 11)
#HandelMuchMany(RootPath, path2, PicoNum = 6, Begin = 20, End = 60, Class = 'R', BeginN=0, EndN=15)