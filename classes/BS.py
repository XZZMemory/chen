class BS:  #还需要实现获取可连接用户功能（通过DM距离矩阵）

    def __init__(self,chaNum,coverArea,power,DL):  #信道数，覆盖面积，功率，距离列表
        self.chaNum = chaNum
        self.coverArea = coverArea
        self.power = power
        self.DL = DL
        self.userUnderCover = []
        self.getUser(DL)
    def getUser(self,DL):
        for i in range(len(DL)):
            if DL[i] < self.coverArea:  #如果用户距离小于覆盖范围
                self.userUnderCover.append(i)  #那么就把用户编号加入到可连接列表