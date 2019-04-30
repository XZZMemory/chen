# -*- coding:utf-8 -*-
import os

def clearTxt(path):
    '''
    清空文件内容
    '''
    try:
        f = open(path, 'w')
        f.close
    
    except IOError:
        print("fail to open file") #打开文件出现错误

def writeTxt(path, data, flag):
    
    '''
    输入要写入的路径 path 和要存入的二维数组 data
    flag为0则覆写，1（或其他）则为追加
    把一段数据（一个二维数组）按行写入一个TXT文件中
    每行写入一个一维数组，用逗号‘，’来分隔
    '''
    try:
        if flag == 0:
            f = open(path, 'w' )
        else:
            f = open(path, 'a')
        for line in data:
            for number in line:
                f.write(str(number))
                f.write(',') #用逗号作为间隔
            f.write('\n')
        f.write('\n')  #最后留一个间隔符，用来表示二维数组之间的间隔
        f.close()
    except IOError:
        print("fail to open file") #打开文件出现错误
        
def readTxt(path, flag):
    '''
    flag表示是否转化为整形，flag为0不转换，为1则按行转换，为2则按矩阵转换
    （取决于同行或同矩阵内是否有浮点数）
    从TXT文件读取文件中的所有二维数组
    path为地址
    返回的data是一个三维数组，其中包括多个读取到的二维数组
    '''
    
    try:
        f = open(path,'r')
        data = []  #完整数据，三维数组，用来储存所有数据
        mat = []  # 单独一个二维数组
        isint = 1
        while True:
            line = f.readline()  #读取一行
            if (not line):  #如果没有行了，表示文件到结尾
                break
            elif (line =='\n' ): #如果是分隔符，表示这个矩阵读取完成，当前矩阵加入data
                if (isint == 1)and(flag ==2):
                    mat = [[int(y) for y in x if y] for x in mat if x] #如果满足条件并且需要转整形
                data.append(mat) 
                mat = []
                continue
            else:
                floatline = line.split(',' )#按照逗号分隔开获取的字符串
                floatline.pop()#删除掉最后一个\n字符      
                if('.' in line):
                    floatline = [float(x) for x in floatline if x]#把每一段字符串，转化为浮点数
                    isint = 0  #isint置为0则该矩阵不能整体转为整形 
                elif(flag == 1):
                    floatline = [int(x) for x in floatline if x]
                elif(flag == 2):
                    floatline = [float(x) for x in floatline if x]#把每一段字符串，转化为浮点数             
                mat.append(floatline)
        
        f.close()
        return data
    except IOError:
        print("fail to open file") #打开文件出现错误
    
