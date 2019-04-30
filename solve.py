# -*- coding:utf-8 -*-
import copy
import math
import random
import time

from algorithm import Zoutendijk


class BS:  #还需要实现获取可连接用户功能（通过DM距离矩阵）

    def __init__(self,chaNum,coverArea,power,DL):  #信道数，覆盖面积，功率，距离列表
        self.chaNum = chaNum
        self.coverArea = coverArea
        self.power = power
        self.DL = DL   #还是需要存一下距离列表来计算干扰距离
        self.userUnderCover = []  #user under cover在覆盖范围内的用户组
        self.getUser(DL)
    def getUser(self,DL):
        for i in range(len(DL)):   #对于所有用户
            if DL[i] < self.coverArea:  #如果用户距离小于覆盖范围
                self.userUnderCover.append(i)  #那么就把用户编号加入到可连接列表
                
                
                
def SetBS (macronum, piconum, infofmbs, infofpbs, DM):

    BSList = []
    for i in range(macronum):    #循环macro基站个数次
        m  = BS(infofmbs[0],infofmbs[1],infofmbs[2],DM[len(BSList)])  #建立macr属性的BS对象 根据现有的基站数（len（BSList））来确定该取DM的哪一行
        BSList.append(m)  #加入到基站列表
    if len(DM)>1:  #如果DM矩阵长度大于1，表示不只有一个pico基站--------------------------------------------------------17.11.07
        for i in range(piconum):  #同上，pico基站个数次
            p  = BS(infofpbs[0],infofpbs[1],infofpbs[2],DM[len(BSList)])
            BSList.append(p)
    
    return BSList

class Individual:
    '''
    individual  传入参数BS为基站类数组
    '''
    
    def __init__(self, BS, userNum, maxc, Qmax, Alpha, B):
        '''
        传入参数，BS为基站类数组，userNum是用户个数，maxc是每个用户可连接最大信道数，pm变异率
        '''
  
        self.BSNum = len(BS)   #基站数量
        self.BS = BS
        self.maxc = maxc  #用户可连接最大信道数
        #self.Qmax = Qmax
        self.Qmax = 1   #-----------------------------------------2017.06.08修改为实数功率等级，所以最大功率等级是1
        self.power = 0  #总功率
        self.rate = 0  #总速率
        self.userRate = 0 #用户平均速率
        self.rp = 0 # r/p，速率功率比值，一个用来衡量解的优劣性的一个指标
        self.userNum = userNum  #用户数量
        self.userList = [0 for i in range(userNum)]  #用户占用的信道数列表
        self.userRateList = [0 for i in range(userNum)]  #----------------------17.10.25新增，用来保存每个用户的速率
        self.minRate = B*0.075  #--------------------------------------需要保证的用户最低速率 17.10.25
        self.tump = self.InitialiseGene()  #生成初始化信道分配
        self.genec = self.tump[0]   #个体基因信道分配 
        self.genep = self.tump[1] #个体基因功率系数         
        self.gener = self.tump[2]  #个体每个信道上的速率
        self.pm = 0.05
        #self.cd = 0 #拥挤度
        self.Alpha = Alpha  #计算噪声参数
        self.B = B  #带宽
        
        self.Bias = 1000000  #防止干扰和噪声过小设置的偏移值，在计算SINR之后会约掉
        self.Bias2 = 100  #用来给zoutengdijk放大maxlambda
        
        self.ZoutendijkCount = 0    #进行下一次Zoutendijk寻优操作需要等待的时间，每次迭代减小1
        self.ZoutendijkFailedTime = 0 #Zoutendijk失败次数，每次失败都要增加等待时间的初始值
        
        self.rank = 0  #初始支配等级
        self.CrowdingDistance = 0 #拥挤度初始化
        self.SequenceNumber = -1 #当前个体的唯一序号
        self.DominateIndividuals = [] #当前支配的个体序号列表
        self.DominatedNumber = 0 #当前支配该个体的个体数
        
        #self.powerLimit = 0.0000001 #功率等级小于这个数就直接认为是0并且释放信道
        self.powerLimit = 0.00001 #功率等级小于这个数就直接认为是0并且释放信道
        
        self.geneIplusN = self.tump[3]
        self.geneSINR = self.tump[4]
        self.genePartialDerivative = self.tump[5]
        
        #以下为了测试
        self.genePartialR =[]
        self.geneI = []
        self.geneN = []
        self.geneTest()
        
        self.LowestRate = 0  #为了单目标优化设定的最低速率
        
    def geneTest(self):
        '''
        为了初始化基因测试数组的函数
        '''
        #以下为了测试
        templist = [[] for i in range(self.BSNum)]  #生成一个空的基因数组
        for i in range(self.BSNum):
            for j in range(self.BS[i].chaNum):
                templist[i].append(0)
                
        self.genePartialR = copy.deepcopy(templist)
        self.geneI = copy.deepcopy(templist)
        self.geneN = copy.deepcopy(templist)
        
    def InitialiseGene(self):  
        '''
        初始化信道分配和功率分配
        
        功率分配的从1开始有点太少了，应该计算出一个最小功率
        '''
        genec = []
        genep = []
        cnum = 0
        for i in range(self.BSNum):  #计算信道总数
            cnum = cnum + self.BS[i].chaNum
        
        p = self.userNum/(cnum*self.BSNum)    #依概率决定是否分配信道，概率为用户数/信道数
        
        for i in range(self.BSNum):  #循环基站个数次
            chan = []   #当前基站的信道分配
            powe = []  #当前基站的功率分配
            wholePower = self.Qmax  #总功率设置为Qmax功率等级
            for j in range (self.BS[i].chaNum):  #循环当前基站信道个数次
                if(random.random() > p)and(len(self.BS[i].userUnderCover)>0):  
                    #如果随机数大于P,并且可选用户列表不为空，则为当前信道分配用户
                    try:
                        user = random.choice(self.BS[i].userUnderCover)
                        chan.append(user)   #随机在可选用户群中选取一个用户分配给当前信道
                        self.userList[user] += 1   #该用户被分配的信道数+1
                        #wp = random.randint(2,int(wholePower/(self.BS[i].chaNum-j)*0.6))-1  #随机生成一个功率值(1到剩余功率的50%+1)
                        wp = random.random()*(wholePower/(self.BS[i].chaNum-j))  #--------------------------------------------2017.06.08修改功率等级为实数
                        #powe.append(random.randint(1,wp)) 
                        powe.append(wp)#将这个功率值分配给这个信道
                        wholePower = wholePower-wp    #可用功率值减少
                    except IndexError:
                        print('user序号:' + str(user))
                        print('userlist列表:' + str(self.userList))
                        print('错误，列表溢出，程序终止')
                        exit(0)
                else:
                    chan.append(-1)   #如果不分配信道，设置为-1
                    powe.append(0)   #不用的信道功率设置为0
            
            genec.append(chan)   #当前基站的信道分配加入基因中
            genep.append(powe)  #当前基站的功率分配加入基因中
        
        gener = []  #初始化gener
        geneSINR = []  #初始化SINR和Interfere，为了计算zoutendijk-------------------------------------------------------2017.06.27新增
        geneIplusN = []   #干扰和噪声的和
        genePartialDerivative = []  #pij偏导数
        for i in range(self.BSNum):  #复制跟基因同样大小的数组
            gener.append([0 for j in range(self.BS[i].chaNum)])
            geneSINR.append([0 for j in range(self.BS[i].chaNum)])
            geneIplusN.append([0 for j in range(self.BS[i].chaNum)])
            genePartialDerivative.append([0 for j in range(self.BS[i].chaNum)])

        tump = (genec,genep,gener,geneIplusN,geneSINR,genePartialDerivative)
        return(tump)
    
    def AddChannel(self,i):  #这个得改
        '''
        增加用户i的信道分配,i是用户在用户数组中的位置，表示第i个用户
        执行策略1和策略2，策略1算法复杂度较低，只寻找基站有空余信道并且用户i在这个基站覆盖范围下的情况
        ，分配一个空余信道给用户i
        策略2复杂度较高，在没有空余信道情况下才执行，剥夺一个多信道用户的信道，分配给i用户
        策略2也不能保证用户一定能分配到信道 best effort
        '''
        c = self.userList[i]  #获取用户i分配到的信道数
        
        for j in range (self.BSNum):  #策略1，在所有基站搜索是否i处于覆盖范围，并且信道有空余
            if (-1 in self.genec[j]) and (i in self.BS[j].userUnderCover):  
                #如果当前基站有空余信道（值为-1）并且用户在这个基站的覆盖范围
                k = self.genec[j].index(-1)    #就找到空余信道
                self.genec[j][k] = i        #把空余信道分配给用户i
                self.userList[i] += 1        #用户分配的信道数+1
                if self.Qmax<sum(self.genep[j]):
                    self.RevisePower(j)             
                #self.genep[j][k] = random.randint(0,int((self.Qmax-sum(self.genep[j]))*0.8))+1  #同时在剩余可用功率中分配给该信道一些功率
                self.genep[j][k] = random.random()*((self.Qmax-sum(self.genep[j]))/self.userNum) #-----------------------------------------2017.06.08修改功率等级为实数
                break  #问题解决，跳出循环
        
        if c == self.userList[i]:  #如果之前的循环没有成功分配信道，进入策略2
            #策略2的思想是，将用户根据分配到的信道数多少排序，从分配信道最多的用户调整信道给需要分配信道的用户
            mlist = []     #多信道用户列表
            m = max(self.userList)     #先获取获得信道最大数
            for n in range(m-1):   #m是信道分配数量的最大值，循环m-1次，循环到信道分配数为2为止
                for r in range(self.userNum):  #对于每个用户
                    if self.userList[r] == m:    #如果信道分配数等于m
                        mlist.append(r)    #加入多信道用户列表
                m -= 1    #循环一次之后，寻找信道分配数少1的用户
            mlen = len(mlist)   #可调整用户数量
            for j in range(self.BSNum):   #开始调整信道，对于每个基站
                if i in self.BS[j].userUnderCover:  #如果i用户在这个基站的覆盖范围内
                    for s in range(mlen):    #对于每一个在mlist里的用户
                        if mlist[s] in self.genec[j]:    #如果该用户也在这个基站的信道分配中
                            k = self.genec[j].index(mlist[s])  #定位这个信道的编号
                            self.genec[j][k] = i     #将这个信道转分配给用户i
                            self.userList[mlist[s]] -= 1  #被剥夺信道的用户分配信道数-1
                            self.userList[i] += 1    #用户i的信道数 +1
                            break    #结束本循环
                if c < self.userList[i]:    #如果在当前基站中，用户i获得了信道
                    break    #结束本循环
                        
    def DecChannel(self,i):
        '''
        减少用户i的信道分配
        2017.3.15 策略修改：寻找速率/功率比最低的那个信道删除
        '''
        rp = float('inf')
        bf = -1  #基站编号
        sf = -1 #信道编号
        for n in range (self.BSNum):
            for s in range(self.BS[n].chaNum):
                catch = [n,s,i]
                try:
                    if (self.genec[n][s] == i) and (self.gener[n][s]/(self.genep[n][s]/self.Qmax*self.BS[n].power)<rp):
                        #如果 信道分配给了用户i，并且信道的速率功率比小于rp，则标记该信道
                        bf = n  #标记这个基站的
                        sf = s  #这个信道
                except ZeroDivisionError:
                    print('减少信道函数出现除以0错误')
                    print(catch)
                    print(self.genec[catch[0]])
                    print(self.genep[catch[0]])
                    exit(0)
        
        self.genec[bf][sf] = -1 #信道置空
        self.genep[bf][sf] = 0 #功率置0
        self.gener [bf][sf] = 0 #速率置为0
        self.userList[i] -= 1 #该用户分配信道数-1
                
                
    def RevisePower(self,i):
        '''
        修复功率，即基站i的总功率不能大于1000
        
        +1那里应该设置个最低功率，1有点太少了
        
        ----------------------------------------------------------------------------2017.06.08修改功率等级为实数导致这个函数重写
        
        '''
        
        for j in range(len(self.genep[i])):  #----------------------------------------2017.06.08新增
            if self.genep[i][j]/self.Qmax<self.powerLimit:  #如果功率非常非常小(小于最低限度)，则直接置为0
                self.genep[i][j] = 0
                self.genec[i][j] = -1  #同时这个信道标注为空闲
            
        a = sum(self.genep[i])  #该基站当前功率的总数，应该小于Qmax才对
        if a >= self.Qmax:    #如果a超过Qmax
            for j in range(len(self.genep[i])): 
                if self.genep[i][j] != 0:
                    self.genep[i][j] *= self.Qmax/a #按照比例减少
                    
            a = sum(self.genep[i])
            for j in range(len(self.genep[i])):
                if self.genep[i][j] != 0:
                    self.genep[i][j] += (self.Qmax-a)  #修正，因为实数计算会不可避免出现不能取整的情况，该操作用于取整，使得基站满载
                    break

    def ReviseUserRate(self, k):
        '''
        用来解决单个用户速率过低的问题----------------------------------------------------------------17.10.25新增
        单纯的提升用户分配到的信道功率
        '''          
        for i in range(self.BSNum,0):
            for j in range(self.BS[i].chaNum,0):
                if self.genec[i][j] == k:
                    self.genep[i][j] *= 2  #从后向前查找用户K分配的信道，功率提升50%（为了给pico基站）
                                   
    def Revise(self):  
        '''
        对个体进行修复,修复的目标是，每个用户至少一个信道（保证一定的带宽），每个用户至多maxc个信道
        每个基站总功率不超过1000
        '''
        self.userList = [0 for i in range(self.userNum)]  #每次修复之前都要重新计算用户占用信道数列表
        for i in range(self.BSNum):  #循环基站数次
            for j in range(self.BS[i].chaNum):  #循环信道数次
                k = self.genec[i][j]   #信道分配给了第K个用户
                if  k != -1:  #如果分配了
                    self.userList[k] += 1  #那么用户k占用的信道数+1
                    
        for i in range(self.BSNum): #先修复一次功率,以免出现溢出
            self.RevisePower(i)
        
        for i in range(self.userNum):   #先修复信道分配
            while(self.userList[i] > self.maxc):  #如果某用户分配了多余maxc条信道，则执行用户信道减少方法
                self.DecChannel(i)
        for i in range (self.userNum):#由于先执行了减少信道，这样增大了空余信道存在的概率，可以提高信道增加的速度
            if self.userList[i] ==0:   #如果某用户没有分配到信道，则执行信道增加方法
                self.AddChannel(i)
        for i in range(self.userNum):  #如果用户速率小于最低速率，那么就执行修复算法，提升pico基站信道的功率-----------------------------------17.10.25新增
            if self.userRateList[i]<self.minRate:
                self.ReviseUserRate(i)
                
        for i in range(self.BSNum): #修复功率
            self.RevisePower(i)
        
        #self.power = self.CalculateTotalPower()  #计算总功率
        #self.rate = self.CalculateTotalRate()  #计算总速率
        #self.rp = self.rate/self.power
        self.CalculateAll()  #计算总速率，总功率，平均速率和速率功率比值
        
    #def Gsnk(self,s,n,k):
        '''
        计算信道功率增益
        '''
        
    def CheckChannel(self):
        '''
        信道检查
        '''
        flag =0
        for i in range(self.userNum):
            if self.userList[i] == 0 or self.userList[i]>3:
                flag = 1
                
        if flag == 1:
            print("信道未修复")
        #else:
           # print("信道分配未发现问题")
    def CalculateInterfere(self, n, s):
        '''
        计算干扰, 第n个基站的第s个信道收到的干扰
        干扰公式  信道功率*距离^（-4）  加法叠加
        '''
        Interfere = 0  #干扰初始设置为0
        for i in range(self.BSNum):   #对于每个基站
            if (i != n)and(self.genec[i][s] != -1):   #如果不是当前基站n，并且要计算的信道功率不为0（即已经分配）
                k = self.genec[n][s]     #取要计算干扰信道链接的用户K
                Distence = self.BS[i].DL[k]  #计算用户K与要计算干扰的基站的距离
                Interfere = Interfere + self.genep[i][s]/self.Qmax*self.BS[i].power*self.Bias*(Distence**(-4))   #用公式计算干扰，加法叠加
        return(Interfere)
    
    def CalculateRate(self, n, s, flag = 0):
        '''
        计算速率，第n个基站的第s个信道的速率
        速率公式   带宽B*log2（(功率*距离^(-4))/(总功率*覆盖范围^（-4）)/alpha+干扰）
        '''
        channalNum = self.BS[n].chaNum  #获取当前基站信道总数
        Interfere = self.CalculateInterfere(n, s)  #计算干扰
        
        self.geneI[n][s] = Interfere  #---------------------------为了测试
        
        k = self.genec[n][s]
        Distence = self.BS[n].DL[k] #得到基站到用户的距离
        g = self.genep[n][s]/self.Qmax*self.BS[n].power*self.Bias*(Distence**(-4))   #计算当前信道功率（17.04.18修正为信号强度）
        Noise = self.Qmax*self.BS[n].power*self.Bias*(self.BS[n].coverArea**(-4))/self.Alpha  #计算噪声
        self.geneN[n][s] = Noise  #----------------------------------为了测试
        
        if Noise<= 0 :
            print("Error")
        IplusN = Interfere + Noise
        #if flag == 2:
            #print("Interefere = %f, Noise = %f, self.Qmax = %f, self.BS[n].power = %f, self.BS[n].coverArea**(-4) = %f"%(Interfere, Noise, self.Qmax, self.BS[n].power, self.BS[n].coverArea**(-4)))
        self.geneIplusN[n][s] = IplusN #-----------------------------------------------------------------------顺带储存了干扰噪声和（17.06.27新增）
        SINR = g/IplusN
        print("base: "+str(n)+" channel: "+str(s)+" sinr: "+str(SINR))
        self.geneSINR[n][s] = SINR
        '''
        d = 1+p/(Noise+Interfere)
        print(p)
        print(self.genep[n][s])
        print(self.Qmax)
        print(self.BS[n].power)
        print(Distence)
        print()
        
        print(Noise)
        print(Interfere)
        print()
        print()
        '''
        if flag == 0:   #---------------------------------------------------------------------------------17.06.28修正，如果信道不分配导致IplusN为0，所以予以修正
            #if SINR<0:
                #print("g = %f, Interfere = %f, Noise = %f" %(g*10000, Interfere*10000, Noise*10000))
            Rate = self.B/channalNum*math.log2(1+SINR) #计算速率
            
            self.gener[n][s] = Rate  #把计算好的速率存入gener中
            self.userRateList[k] += Rate  #------------------------------------------------------------------速率累加进用户速率 17.10.25
            
            return(Rate)
    def CalculateTotalPower(self):  
        '''
        计算个体的总功率
        '''
        TotalPower = 0
        for i in range(self.BSNum):  #对于每一个基站
            for j in range(self.BS[i].chaNum):  #的每一条信道
                if abs(self.genep[i][j])<self.powerLimit:   #如果这条信道功率非常接近0  -------------------------------------------------------2017.07.11修改防止误差
                    self.genep[i][j] = 0
                    self.genec[i][j] = -1
                else:  
                    TotalPower = TotalPower + self.genep[i][j]/self.Qmax*self.BS[i].power  #那么就计算功率累加至totalpower
        
        return(TotalPower)
    
    def CalculateTotalRate (self):  
        '''
        计算个体的总速率
        '''
        TotalRate = 0
        for i in range(self.userNum):   #---------------------------------------17.10.25新增，首先把用户速率列表置为0
            self.userRateList[i] = 0
        for i in range(self.BSNum):  #对于每一个基站
            for j in range(self.BS[i].chaNum):  #的每一条信道
                if self.genec[i][j] != -1:     #如果这条信道不是空置
                    rateij = self.CalculateRate(i, j)
                    TotalRate += rateij   #那么就计算速率累加至totalrate
                else:
                    self.gener[i][j] = 0 #空置信道速率为0
                    self.genep[i][j]= 0 #功率也置为0
                    #self.CalculateRate(i, j, 1) #否则计算该信道的干扰和噪声
                    
        return(TotalRate)
    
    def CalculatePartialDerivative(self):
        '''
        计算偏导数
        '''
        for i in range(self.BSNum):   #对于每一个基站
            for j in range(self.BS[i].chaNum):  #的每一个信道计算偏导数
                SumParts = 0  #求和部分初始化为0
                k = self.genec[i][j]
                gij = self.BS[i].power*self.Bias*(self.BS[i].DL[k])**(-4)  #公因子gik
                for m in range(self.BSNum):  #循环基站个数次
                    l = self.genec[m][j]  #取信道关联的用户
                    if l != -1:  #如果信道不为闲置
                        if self.geneIplusN[m][j] <= 0:    #用来检测IplusN大小的
                            #print("before self.geneIplusN[%d][%d] = %f" %(m, j, self.geneIplusN[m][j]))
                            self.CalculateRate(m, j, 2)#-------------------------------------------------------------------------------这里导致了一个print
                            #print("after self.geneIplusN[%d][%d] = %f" %(m, j, self.geneIplusN[m][j]))
                        if m == i:  #如果是基站i
                            SumParts += -1/self.geneIplusN[i][j]*(1+self.geneSINR[i][j])**(-2)  #gik/(I+N)
                        else: #否则
                            SumParts += self.genep[m][j]*self.Bias*self.BS[m].power*(self.BS[m].DL[k]**(-4))/(self.geneIplusN[m][j]**2)*(1+self.geneSINR[m][j])**(-2) #gmk*pmk*gik/(I+N)^2--------------------17.06.27除以0问题
                PartderviR = self.B/self.BS[i].chaNum/math.log(2)*gij*SumParts  #R偏导数
                self.genePartialR[i][j] = PartderviR
                
                #Partdervi = (PartderviR*self.power-(self.BS[i].power**0.5*0.1)*self.rate)/(self.power**2)  #计算偏导数  对power开方除以10是为了减少功率造成的影响，尽量优化速率
                Partdervi = (PartderviR*math.log2(1+self.power)-((1/(1+self.power)/math.log(2))*self.BS[i].power)*self.rate)/((math.log2(1+self.power))**2) #---------------------------17.11.04修改，因为R=log2（1+sinr） R/log2（1+P）才是线性的
                self.genePartialDerivative[i][j] = Partdervi  #存入值
    def RunningZoutendijk(self):
        '''
        进行Zoutendijk可行方向法
        '''
        self.CalculatePartialDerivative()  #先计算偏导数
        
        expand = self.Bias2  #扩张值，用来放大数据使得近似值的影响变小
        feasibleX = []
        deltaF = []
        matrixA = [[] for i in range(self.BSNum)]
        vectorb = [-1*expand for i in range(self.BSNum)]
        xUpBounds = []
        xLowBounds = []
        
        for i in range(self.BSNum):
            for j in range(self.BS[i].chaNum):
                if self.genec[i][j] != -1:   #只有分配出去的信道才进行计算
                    feasibleX.append(self.genep[i][j]*expand)
                    deltaF.append(self.genePartialDerivative[i][j])
                    for k in range(self.BSNum):
                        if i == k:
                            matrixA[k].append(-1)
                        else:
                            matrixA[k].append(0)
        n = len(feasibleX)

        xUpBounds = [1*expand for i in range(n)]
        xLowBounds = [0 for i in range(n)]    #x上下界
        Zouten = Zoutendijk.ZoutendijkProblem(feasibleX=feasibleX, deltaF=deltaF, matrixA=matrixA, vectorb=vectorb, xLowBounds=xLowBounds, xUpBounds=xUpBounds, sense='Maximize')
        solve = Zouten.BoundarySolve() #进行Zoutendijk
        MaxLambda = solve[0]
        vectord = solve[1]
        
        return (MaxLambda, vectord)
    
    def MatrixPlus(self, A1, A2, Pram = 1):
        '''
        矩阵相加函数  A1 = A1+Pram*A2*A3
        '''
        AR = []
        for i in range(len(A1)):
            aa = []
            for j in range(len(A1[i])):
                a = A1[i][j] + A2[i][j]*  Pram
                if a < 0:
                    print("MatrixPlus: %d,%d, A1 %f,A2 %f , Lambda %f导致小于0" %(i, j, A1[i][j], A2[i][j], Pram))
                aa.append(a)
            AR.append(aa)
                
        
        return AR
    
    def TestFunction(self):
        '''
        为了测试X = X+λ*d的函数性质，实验证明了是关于lambda的单调递增函数
        '''
        solve = self.RunningZoutendijk()
        MaxLambda = solve[0]
        print("MaxLambda = %f" %(MaxLambda))
        vectord = solve[1]
        
        Matrixd = []   
        t = 0 
        for i in range(self.BSNum):  #这一步是把vectord转换成可以直接和genep对应相加的矩阵
            temp = []
            for j in range(self.BS[i].chaNum):
                if self.genec[i][j] != -1:
                    temp.append(vectord[t])
                    t+=1
                else:
                    temp.append(0)
            Matrixd.append(temp)
            
        power = copy.deepcopy(self.genep)
        #PD = copy.deepcopy(self.genePartialDerivative)
        
        rplist = []
        genep = copy.deepcopy(self.genep) #用作测试，保存备份
        
        for i in range(20):
            self.genep = self.MatrixPlus(power, Matrixd, i*0.05*MaxLambda/self.Bias2)
            self.CalculateAll()
            rplist.append(self.rp)
            
        self.genep = genep   #恢复备份
        self.CalculateAll()
            
        #print("RPLIST:")
        #print(rplist)
        return rplist
        
        
    def  CheckChannelAndPower(self, flag = 0):
        '''
        用来检测功率跟信道分配不符的函数
        '''    
        for i in range(self.BSNum):
            for j in range(self.BS[i].chaNum):
                if self.genec[i][j] == -1 and self.genep[i][j] != 0:
                    #print("发现信道分配与功率不符，位置在基站%d，信道%d" %(i, j))
                    if flag == 0:
                        return 1
                    else:
                        self.genep[i][j] =0
                        #print("静默修复")
        
        return 0  
        
    def ZoutendijkReplace(self, times = 3):
        '''
        执行Zoutendijk替换原来的解
        '''
        self.CalculateAll()
        
        rate = self.rate
        power = self.power
        DominateValue = 0  #储存支配值
        DominateGenep = copy.deepcopy(self.genep) #储存最好的基因解
        DominateGenec = copy.deepcopy(self.genec)
        DominateFlag = 0 #是否产生过支配解
        genepower = copy.deepcopy(self.genep)
        rp = self.rp
        
        
        #rplists = []
        
        for k in range(times):  #执行Zoutendijk的次数
            solve = self.RunningZoutendijk()
            MaxLambda = solve[0]
            vectord = solve[1]
            
            
            #rplists.append(self.TestFunction())#----------------------------------------------用作测试Zoutendijk中能效比的单调性
            t = 0 #计数器
            for i in range(self.BSNum):
                for j in range(self.BS[i].chaNum):  #二重循环执行X = X+λ*d
                    if self.genec[i][j] != -1:
                        self.genep[i][j] = self.genep[i][j] + vectord[t]*0.8*MaxLambda/self.Bias2
                        if abs(self.genep[i][j])<self.powerLimit:  #强制归零
                            self.genep[i][j] = 0
                            self.genec[i][j] = -1
                        
                        if self.genep[i][j]< 0:
                            print("error, genec = %d, genep = %f, genepower = %f, vctord = %f, maxlambda = %f"%(self.genec[i][j],self.genep[i][j], genepower[i][j], vectord[t], MaxLambda))
                        t += 1
            self.CalculateAll()
        '''
        #验证能效比和λ取max的实验写入部分 
        rpsave = r'E:\nsgaii\indi\rplist.txt'   
        NSGA_II.RnWLists.writeTxt(rpsave, data = rplists, flag = 1)
        print("写入地址%s"%(rpsave))
        '''
        self.Revise()
        if self.rp>rp and self.rate>self.LowestRate:
            DominateFlag = 1
        else:
            self.genep = DominateGenep
            self.CalculateAll()
        return DominateFlag  
        
        
    def CalculateAll(self):
        '''
        一次计算所有需要计算的数值
        '''
        self.rate = self.CalculateTotalRate()  #算速率
        self.power = self.CalculateTotalPower()  #算功率
        self.userRate = self.rate/self.userNum  #算平均每个用户速率
        self.rp = self.rate/self.power #算速率功率比值


    def MutationOnceSingleTarget(self, TargetFlag = 0):  
        '''    
        个体变异，如果（0,1）随机数p小于变异率，则个体发生变异
        首先随机选择一条子基因，然后在子基因上确定变异的位置
        （其实就是在基因上选择变异位，因为基因有多条，所以称作子基因）
        变异规则：1，如果选择的信道是空信道，随机分配一个可行的用户，随机分配一个可行的功率，然后修复
                          2，如果信道非空，再进行一次判断
                              a，如果随机数p小于变异率/2（其实就是变异条件下的1/2概率，省事）就把信道置空，功率清零
                              b，否则就强行给这个信道随机换一个用户，然后再随机换一个功率
        '''
        p = random.random()  #首先生成一个（0,1）之间的随机数
        if p<self.pm:  #如果随机数小于变异率，则执行变异过程
            k = random.randint(0,self.BSNum-1)  #先随机选择子基因
            l = random.randint(0,len(self.genec[k])-1)  #再寻找子基因上的变异位
            if TargetFlag == 0:  #表示是信道分配
                if self.genec[k][l] == -1:    #如果信道是空的
                    if len(self.BS[k].userUnderCover) != 0:  #可选用户列表不为0才能操作
                        self.genec[k][l] = random.choice(self.BS[k].userUnderCover)  #在可选范围内随机选择一个用户分配给这个信道

                else:    #如果信道非空
                    if p<self.pm/2:   #这是个50%的概率
                        self.genec[k][l] = -1   #信道置空
                        self.genep[k][l] = 0   #功率置0
                        #self.Revise()   #修复
                    else:    #另外50%
                        self.genec[k][l] = random.choice(self.BS[k].userUnderCover)   #随机找一个可选用户
                        #self.genep[k][l] = random.randint(1,1000-sum(self.genep[k])) #随机分配一个可行功率
                        #self.Revise()   #修复
            else:
                if self.genec[k][l] == -1:    #如果信道是空的
                    self.MutationSingleTarget(TargetFlag = 1) #再变异一次功率
                else:
                    self.genep[k][l] = random.random()*(self.Qmax-sum(self.genep[k])-self.powerLimit)+self.powerLimit #再可选范围内随机分配一个功率
                    self.Revise()  #对基因执行修复
                

    def MutationOnce(self):  
        '''
        新版变异，变异操作由外部控制，个体自身不携带变异率，通过拥挤度来控制变异
        
        个体变异，如果（0,1）随机数p小于变异率，则个体发生变异
        首先随机选择一条子基因，然后在子基因上确定变异的位置
        （其实就是在基因上选择变异位，因为基因有多条，所以称作子基因）
        变异规则：1，如果选择的信道是空信道，随机分配一个可行的用户，随机分配一个可行的功率，然后修复
                          2，如果信道非空，再进行一次判断
                              a，如果随机数p小于变异率/2（其实就是变异条件下的1/2概率，省事）就把信道置空，功率清零
                              b，否则就强行给这个信道随机换一个用户，然后再随机换一个功率
        '''
        p = random.random()  #首先生成一个（0,1）之间的随机数
        k = random.randint(0,self.BSNum-1)  #先随机选择子基因
        l = random.randint(0,len(self.genec[k])-1)  #再寻找子基因上的变异位
        if self.genec[k][l] == -1:    #如果信道是空的
            if len(self.BS[k].userUnderCover) != 0:  #可选用户列表不为0才能操作
                #print(self.Qmax-sum(self.genep[k]))
                self.genec[k][l] = random.choice(self.BS[k].userUnderCover)  #在可选范围内随机选择一个用户分配给这个信道
                self.userList[self.genec[k][l]] += 1  #该用户占有的信道数 +1
                if self.Qmax-sum(self.genep[k])<self.powerLimit*100:  #-----------------------------------------------------------2017.07.24修改为0，如果功率溢出就修复
                    self.RevisePower(k)   #功率超了修复
                #self.genep[k][l] = random.randint(1,self.Qmax-sum(self.genep[k])) #再可选范围内随机分配一个功率
                self.genep[k][l] = random.random()*(self.Qmax-sum(self.genep[k])) #-----------------------------------------------------2017.06.08功率等级修改为实数

        else:    #如果信道非空
            if p<0.01:   #这是个20%的概率
                self.genec[k][l] = -1   #信道置空
                self.genep[k][l] = 0   #功率置0
            else:    #另外80%
                self.genec[k][l] = random.choice(self.BS[k].userUnderCover)   #随机找一个可选用户
                self.userList[self.genec[k][l]] += 1  #该用户占有的信道数 +1
                if self.Qmax - sum(self.genep[k])<self.powerLimit*100:
                    self.RevisePower(k)
                #self.genep[k][l] = random.randint(1,self.Qmax-sum(self.genep[k])) #随机分配一个可行功率
                self.genep[k][l] = random.random()*(self.Qmax-sum(self.genep[k])) #-----------------------------------------------------2017.06.08功率等级修改为实数
                
    def Mutation(self, n=1):
        '''
        多次变异
        '''
        for i in range(n):
            self.MutationOnce()
            
    def MutationSingleTarget(self, n=1, TargetFlag = 0):
        '''
        单目标多次变异
        '''
        for i in range(n):
            self.MutationOnceSingleTarget(TargetFlag = TargetFlag)

def Crossover(x, y, n, P = 0.95):
    '''
    交叉函数,输入两个个体X,Y，执行交叉，返回clist是两个子代个体
    对每个子基因执行单点交叉，全部执行完毕之后修复个体
    默认交叉率为0.95
    '''
    Timer = 0
    Begin = time.clock()
    xc = copy.deepcopy(x)  # 深拷贝父代个体x
    yc = copy.deepcopy(y) #深拷贝父代个体y

    if random.random()<=0.95:  #如果满足交叉率条件
        
        t = x.BSNum  #获取子基因个数
        for i in range(t):   #循环子基因个数次，每个基因进行交叉
            l = len(x.genec[i])  #获取当前子基因长度
            
            r = random.randint(1,l-1)  #获取一个随机数，即交叉点
            for j in range(r,l-1):  #交叉点后的进行交叉
                xc.genec[i][j] = y.genec[i][j]  #交叉信道
                xc.genep[i][j] = y.genep[i][j]  #交叉功率
                yc.genec[i][j] = x.genec[i][j]  #同上
                yc.genep[i][j] = x.genep[i][j]  #同上
        
        
        xc.Revise()  #修复
        yc.Revise()   
    End = time.clock()
    Timer = End-Begin
        
    clist = [x,y,xc,yc,Timer]   #两个新生个体
    return clist   #返回


def Crossover2(x, y, n, p = 0.95):
    '''
    重写的交叉函数，可以多点交叉
    
    交叉函数,输入两个个体X,Y，执行交叉，返回clist是两个子代个体
    对每个子基因执行单点交叉，全部执行完毕之后修复个体
    '''
    xc = copy.deepcopy(x)  # 深拷贝父代个体x
    yc = copy.deepcopy(y) #深拷贝父代个体y
    if random.random()<=0.95:
        t = x.BSNum  #获取子基因个数
        l = len(x.genec[0])  #获取当前子基因长度
        for i in range(t):   #循环子基因个数次，每个基因进行交叉
            List = [k for k in range(l)]  #基因位置编号
            slice = random.sample(List, n) #随机选取n个
            flag = 0
            for j in range(l):  #循环子基因个数次            
                if j in slice:    #如果i在选择的序列中
                    flag = 0 if flag == 1 else 1  #改变交叉标志
                if flag == 1:   #交叉标志是1的话
                    '''
                    temp = xc.genec[i][j]
                    xc.genec[i][j] = yc.genec[i][j]
                    yc.genec[i][j] = temp
                    
                    temp = xc.genep[i][j]
                    xc.genep[i][j] = yc.genep[i][j]
                    yc.genep[i][j] = temp
                    '''
                    xc.genec[i][j] = y.genec[i][j]
                    xc.genep[i][j] = y.genep[i][j]
                    
                    if xc.genec[i][j] == -1 and xc.genep[i][j] !=0:
                        print("交叉时产生的x个体发现信道和功率不符")
                        if y.genec[i][j] == -1 and y.genep[i][j] !=0:
                            print("交叉时复制系因为给xc的父代个体y的发现信道和功率不符")
                    
                    yc.genec[i][j] = x.genec[i][j]
                    yc.genep[i][j] = x.genep[i][j]
                    
                    if yc.genec[i][j] == -1 and yc.genep[i][j] !=0:
                        print("交叉时产生的y个体发现信道和功率不符")
                        if x.genec[i][j] == -1 and x.genep[i][j] !=0:
                            print("交叉时复制系因为给yc的父代个体x的发现信道和功率不符")
        xc.Revise()
        yc.Revise()
    clist = [x,y,xc,yc]   #两个新生个体
    return clist   #返回
pow=4.72
def pso(Population):
    '''
    对种群进行pso
    '''
    n = len(Population)
    count = 0
    m = 0
    print("对种群进行pso")
    # print(n)
    for i in range(n):
        if Population[i].CheckChannelAndPower(1) == 1:
            print("ZoutendijkReplace前检测到信道和功率不符")
            # PrintInfo(Population[i])

        # print("对个体执行Zoutendinjk")
        # ZoutendijkFlag = Population[i].ZoutendijkReplace(times = 10)

        if Population[
            i].ZoutendijkCount == 0:  # -------------------------------------------------------------------17.11.04退避算法在这里实现了

            ZoutendijkFlag = Population[i].ZoutendijkReplace(times=5)  # 如果计数结束，那么计算支配解
            m += 1
            count += ZoutendijkFlag  # 统计数量
            if ZoutendijkFlag == 0:  # 如果执行计算没有产生支配解
                Population[i].ZoutendijkFailedTime += 1  # 那么寻找失败次数加1
                Population[i].ZoutendijkCount = Population[i].ZoutendijkFailedTime  # 等待次数等于寻找失败次数，这样会越等越久
        else:
            Population[i].ZoutendijkCount -= 1  # 如果不执行，那么就等待次数减少1

        if Population[i].CheckChannelAndPower(1) == 1:
            print("ZoutendijkReplace后检测到信道和功率不符")
        Population[i].Revise()  # 修复个体
        Population[i].CheckChannel()
        # PrintInfo(Population[i])

    if m == 0:
        return -1  # 如果没执行计算就返回-1
    else:
        return count / m  # 返回支配比率


def PopulationZoutendijk(Population):
    '''
    对种群进行Zoutendijk
    '''
    n = len(Population)
    count = 0
    m = 0
    # print("对种群进行Zoutendijk")
    #print(n)
    for i in range(n):
        if Population[i].CheckChannelAndPower(1) == 1:
            print("ZoutendijkReplace前检测到信道和功率不符")     
            #PrintInfo(Population[i])  
             
        #print("对个体执行Zoutendinjk")    
        #ZoutendijkFlag = Population[i].ZoutendijkReplace(times = 10)     
        
        if Population[i].ZoutendijkCount ==0:   #-------------------------------------------------------------------17.11.04退避算法在这里实现了
            
            ZoutendijkFlag = Population[i].ZoutendijkReplace(times = 10) #如果计数结束，那么计算支配解
            m += 1
            count +=  ZoutendijkFlag #统计数量
            if ZoutendijkFlag == 0:  #如果执行计算没有产生支配解
                Population[i].ZoutendijkFailedTime += 1 #那么寻找失败次数加1
                Population[i].ZoutendijkCount = Population[i].ZoutendijkFailedTime #等待次数等于寻找失败次数，这样会越等越久
        else:
            Population[i].ZoutendijkCount -=1  #如果不执行，那么就等待次数减少1
        
        
        if Population[i].CheckChannelAndPower(1) == 1:
            print("ZoutendijkReplace后检测到信道和功率不符")    
            #PrintInfo(Population[i])
        Population[i].Revise()  #修复个体
        Population[i].CheckChannel()
        #PrintInfo(Population[i])
    
    if m == 0:
        return -1  #如果没执行计算就返回-1
    else:
        return count/m #返回支配比率
        
    #return count/n

def RunSubPopulation(SubPopulation):
    '''
    为了多进程执行写的函数
    '''
    flag = 1 #-------------------------------------------------------------------------------17.11.04手动开关Zoutendijk，0是开启，其他是不开启    
    PopulationNumber = len(SubPopulation)
    
    if flag == 0:
        begin = time.clock()
        m = PopulationZoutendijk(SubPopulation)  
        end = time.clock()   #对Zoutendijk计时
        if m!= -1:
            print("本次Zoutendijk个体个数为%d，用时%d秒，支配率为%f" %(PopulationNumber, end-begin, m))

    GeneLength = len(SubPopulation[0].genec[0])
    NewSubPopulation = []
    for j in range(int(len(SubPopulation)/2)):
        k = random.randint(0, len(SubPopulation)-1)
        X = SubPopulation.pop(k)
        k = random.randint(0, len(SubPopulation)-1)
        Y = SubPopulation.pop(k)   #随机选择XY交叉 
        
        #PrintInfo(X)
        #PrintInfo(Y)
        
        #n = 1
        n = random.randint(1, int(GeneLength/30))#交叉次数为基因长度除以10,取整
        Children = Crossover2(X, Y, n) #其中包括父代两个个体和子代两个个体一共四个[x,y,xc,yc]
        if Children[2].CheckChannelAndPower() == 1:
            print("交叉产生子代后检测到个体XC信道和功率不符")
            PrintInfo(Children[2])
        if Children[3].CheckChannelAndPower(1) == 1:
            print("交叉产生子代后检测到个体YC信道和功率不符")
            PrintInfo(Children[3])
        
        if random.random()<=(1-Children[0].CrowdingDistance/2)/11:  #依照拥挤度概率变异
            Children[0].Mutation(int(GeneLength/30))
            Children[0].Revise
            if Children[0].CheckChannelAndPower(1) == 1:
                print("父代个体X变异后检测到信道和功率不符")
                PrintInfo(Children[0])
        if random.random()<=(1-Children[1].CrowdingDistance/2)/11:
            Children[1].Mutation(int(GeneLength/30))  #对父代变异，因为子代拥挤度值无效
            Children[1].Revise
            if Children[1].CheckChannelAndPower(1) == 1:
                print("父代个体Y变异后检测到信道和功率不符")
                PrintInfo(Children[1])
        NewSubPopulation += Children
    
    
    return NewSubPopulation


def Crossover2SingleTarget(x, y, n, p = 0.95, TargetFlag = 0):
    '''
    重写的交叉函数，可以多点交叉
    
    交叉函数,输入两个个体X,Y，执行交叉，返回clist是两个子代个体
    对每个子基因执行单点交叉，全部执行完毕之后修复个体
    
    单目标，如果TargetFlag为0交叉信道，否则交叉功率
    '''
    xc = copy.deepcopy(x)  # 深拷贝父代个体x
    yc = copy.deepcopy(y) #深拷贝父代个体y
    if random.random()<=0.95:
        t = x.BSNum  #获取子基因个数
        l = len(x.genec[0])  #获取当前子基因长度
        for i in range(t):   #循环子基因个数次，每个基因进行交叉
            List = [k for k in range(l)]  #基因位置编号
            slice = random.sample(List, n) #随机选取n个
            flag = 0
            for j in range(l):  #循环子基因个数次            
                if j in slice:    #如果i在选择的序列中
                    flag = 0 if flag == 1 else 1  #改变交叉标志
                if flag == 1:   #交叉标志是1的话
                    '''
                    temp = xc.genec[i][j]
                    xc.genec[i][j] = yc.genec[i][j]
                    yc.genec[i][j] = temp
                    
                    temp = xc.genep[i][j]
                    xc.genep[i][j] = yc.genep[i][j]
                    yc.genep[i][j] = temp
                    '''
                    if TargetFlag ==0:   #单目标优化，如果flag为0，交叉信道，否则交叉功率
                        xc.genec[i][j] = y.genec[i][j]
                        yc.genec[i][j] = x.genec[i][j]
                    else:
                        xc.genep[i][j] = y.genep[i][j]
                        yc.genep[i][j] = x.genep[i][j]
                    
                    '''
                    if xc.genec[i][j] == -1 and xc.genep[i][j] !=0:
                        print("交叉时产生的x个体发现信道和功率不符")
                        if y.genec[i][j] == -1 and y.genep[i][j] !=0:
                            print("交叉时复制系因为给xc的父代个体y的发现信道和功率不符")
                            '''
                    
                    
                    
                    
                    '''
                    if yc.genec[i][j] == -1 and yc.genep[i][j] !=0:
                        print("交叉时产生的y个体发现信道和功率不符")
                        if x.genec[i][j] == -1 and x.genep[i][j] !=0:
                            print("交叉时复制系因为给yc的父代个体x的发现信道和功率不符")
                            '''
        if TargetFlag != 0: #如果优化的是功率，则修复
            xc.Revise()
            yc.Revise()
    clist = [x,y,xc,yc]   #两个新生个体
    return clist   #返回
   
def ChooseSingleTarget(Population, m):
    '''
    单目标选择
    '''
    for i in range(len(Population)):
        if Population[i].rate<Population[i].LowestRate:
            Population[i].power += Population[i].BS[0].power  #惩罚函数，如果速率不够就加功率
    Population.sort(key = lambda x:x.power)
    NewPop = []
    for i in range(m):
        NewPop.append(Population[i])
        Population[i].CalculateAll()
    return NewPop
def ga_rate_convergence(m, BSList, g, userNumber, LowestRate=0.5, maxchannal=3, Qmax=1, Alpha=10, B=20,
                             ConvergenceFlag=0):
    '''
    单目标优化遗传算法，采用Zoutendijk优化
    '''
    if ConvergenceFlag == 1:
        res = []  # 记录收敛性
    Population = []
    for i in range(m):  # 产生m个个体形成种群

        indi = Individual(BSList, userNumber, maxchannal, Qmax, Alpha, B)
        #    def __init__(self, BS, userNum, maxc, Qmax, Alpha, B):
        indi.CalculateTotalRate  # 计算一下速率，为gener赋值
        indi.Revise()
        indi.CalculateAll()
        indi.LowestRate = LowestRate * userNumber
        Population.append(indi)

    GeneLength = len(Population[0].genec[0])

    best = copy.deepcopy(Population[0])
    best.power = 26
    if ConvergenceFlag == 1:
        res.append(best.rp)
    print("GA收敛")
    for i in range(g):
        # print("开始第%d次迭代" % (i))
        NewPopulation = []
        for j in range(int(len(Population) / 2)):
            k = random.randint(0, len(Population) - 1)
            X = Population.pop(k)
            k = random.randint(0, len(Population) - 1)
            Y = Population.pop(k)  # 随机选择XY交叉

            Children = Crossover2SingleTarget(X, Y, 1, TargetFlag=1)  # 其中包括父代两个个体和子代两个个体一共四个[x,y,xc,yc]

            # Children[0].Mutation(6)
            Children[0].MutationSingleTarget(1, TargetFlag=1)
            Children[0].Revise()
            # Children[1].Mutation(6)  #对父代变异
            Children[1].MutationSingleTarget(1, TargetFlag=1)
            Children[1].Revise()

            NewPopulation += Children

        for j in range(len(NewPopulation)):
            NewPopulation[j].Revise()  # 修复个体
            NewPopulation[j].CalculateAll()

        Population = ChooseSingleTarget(NewPopulation, m)
        # PopulationZoutendijk(Population)


        for j in range(len(Population)):

            if Population[j].rate > Population[j].LowestRate and Population[j].power < best.power:
                best = copy.deepcopy(Population[j])
        if ConvergenceFlag == 1 and i%5==0:
            res.append(best.rp)

    # indiSavePath = './data/indi/SingleTarget.txt'
    # print("开始写入个体")
    # RnWLists.writeTxt(indiSavePath, best.genec, flag=0)
    # RnWLists.writeTxt(indiSavePath, best.genep, flag=1)

    if ConvergenceFlag == 1:
        return res
    else:
        return best



   
def random_search(m,BSList,g, userNumber, LowestRate=1.4,  maxchannal = 3, Qmax = 1, Alpha = 10, B = 20):
    Population = []
    for i in range(m):  # 产生m个个体形成种群

        indi = Individual(BSList, userNumber, maxchannal, Qmax, Alpha, B)
        #    def __init__(self, BS, userNum, maxc, Qmax, Alpha, B):
        indi.CalculateTotalRate()  # 计算一下速率，为gener赋值
        indi.Revise()
        indi.CalculateAll()
        indi.LowestRate = LowestRate * userNumber
        Population.append(indi)

    GeneLength = len(Population[0].genec[0])

    best = copy.deepcopy(Population[0])
    best.power = 26
    for j in range(len(Population)):
        if Population[j].rate > Population[j].LowestRate and Population[j].power < best.power:
            best = copy.deepcopy(Population[j])
    return best.power
    
def ga(m,BSList,g, userNumber, LowestRate=1.4,  maxchannal = 3, Qmax = 1, Alpha = 10, B = 20, ConvergenceFlag = 0):
    '''
    单目标优化遗传算法
    '''
    if ConvergenceFlag == 1:
        res = []  #记录收敛性
    Population = []
    for i in range (m):  #产生m个个体形成种群

        indi = Individual(BSList, userNumber,maxchannal, Qmax, Alpha, B)
        #    def __init__(self, BS, userNum, maxc, Qmax, Alpha, B):
        indi.CalculateTotalRate()  #计算一下速率，为gener赋值
        indi.Revise()
        indi.CalculateAll()
        indi.LowestRate = LowestRate*userNumber
        Population.append(indi)
    
    GeneLength = len(Population[0].genec[0])
    
    best = copy.deepcopy(Population[0])
    best.power = 26
    if ConvergenceFlag == 1:
        res.append(best.power)
    print("开始执行遗传算法")
    for i in range(g):
        # print("开始第%d次迭代"%(i))
        NewPopulation = []
        for j in range(int(len(Population)/2)):    
            k = random.randint(0, len(Population)-1)
            X = Population.pop(k)
            k = random.randint(0, len(Population)-1)
            Y = Population.pop(k)   #随机选择XY交叉 
            
            n = random.randint(1, int(GeneLength/10))#交叉次数为基因长度除以10,取整
            Children = Crossover2SingleTarget(X, Y, n) #其中包括父代两个个体和子代两个个体一共四个[x,y,xc,yc]
            
            Children[0].MutationSingleTarget(6)
            Children[0].Revise()
            Children[1].MutationSingleTarget(6)  #对父代变异
            Children[1].Revise()
    
            NewPopulation += Children
            
        if len(NewPopulation)<m:
            print("错误，当前NewPopulation大小为%d，目标种群大小为%d"%(len(NewPopulation),m))
         
        for j in range(len(Population)):
            Population[j].Revise()  #修复个体
            Population[j].CalculateAll()
         
            
        Population = ChooseSingleTarget(NewPopulation, m)
        for j in range(len(Population)):
            if Population[j].rate>Population[j].LowestRate and Population[j].power<best.power and Population[j].power > pow:
                best = copy.deepcopy(Population[j])
        if ConvergenceFlag == 1 and i % 5 == 0:
            res.append(best.power)
                
    if ConvergenceFlag == 1:
        return res             
    else:
        return best.power
pow=1.226
def ga_pso_rate_convergence(m, BSList, g, userNumber, LowestRate=0.5, maxchannal=3, Qmax=1, Alpha=10, B=20,
                             ConvergenceFlag=0):
    '''
    单目标优化遗传算法，采用Zoutendijk优化
    '''
    if ConvergenceFlag == 1:
        res = []  # 记录收敛性
    Population = []
    for i in range(m):  # 产生m个个体形成种群

        indi = Individual(BSList, userNumber, maxchannal, Qmax, Alpha, B)
        #    def __init__(self, BS, userNum, maxc, Qmax, Alpha, B):
        indi.CalculateTotalRate  # 计算一下速率，为gener赋值
        indi.Revise()
        indi.CalculateAll()
        indi.LowestRate = LowestRate * userNumber
        Population.append(indi)

    GeneLength = len(Population[0].genec[0])

    best = copy.deepcopy(Population[0])
    best.power = 26
    if ConvergenceFlag == 1:
        res.append(best.power)
    print("PSO收敛")
    for i in range(g):
        # print("开始第%d次迭代" % (i))
        NewPopulation = []
        for j in range(int(len(Population) / 2)):
            k = random.randint(0, len(Population) - 1)
            X = Population.pop(k)
            k = random.randint(0, len(Population) - 1)
            Y = Population.pop(k)  # 随机选择XY交叉

            n = random.randint(1, int(GeneLength / 10))  # 交叉次数为基因长度除以10,取整
            Children = Crossover2SingleTarget(X, Y, n, TargetFlag=1)  # 其中包括父代两个个体和子代两个个体一共四个[x,y,xc,yc]

            # Children[0].Mutation(6)
            Children[0].MutationSingleTarget(6, TargetFlag=1)
            Children[0].Revise()
            # Children[1].Mutation(6)  #对父代变异
            Children[1].MutationSingleTarget(6, TargetFlag=1)
            Children[1].Revise()

            NewPopulation += Children

        for j in range(len(NewPopulation)):
            NewPopulation[j].Revise()  # 修复个体
            NewPopulation[j].CalculateAll()

        Population = ChooseSingleTarget(NewPopulation, m)
        # PopulationZoutendijk(Population)


        for j in range(len(Population)):

            if Population[j].rate > Population[j].LowestRate and Population[j].power < best.power and Population[j].power > pow:
                best = copy.deepcopy(Population[j])
        if ConvergenceFlag == 1 and i%5==0:
            res.append(best.power)


    if ConvergenceFlag == 1:
        return res
    else:
        return best

def ga_pso(m, BSList, g, userNumber, LowestRate=1.4, maxchannal=3, Qmax=1, Alpha=10, B=20,
                             ConvergenceFlag=0):
    '''
    单目标优化遗传算法，采用Zoutendijk优化
    '''
    if ConvergenceFlag == 1:
        res = []  # 记录收敛性
    Population = []
    for i in range(m):  # 产生m个个体形成种群

        indi = Individual(BSList, userNumber, maxchannal, Qmax, Alpha, B)
        #    def __init__(self, BS, userNum, maxc, Qmax, Alpha, B):
        indi.CalculateTotalRate  # 计算一下速率，为gener赋值
        indi.Revise()
        indi.CalculateAll()
        indi.LowestRate = LowestRate * userNumber
        Population.append(indi)

    GeneLength = len(Population[0].genec[0])

    best = copy.deepcopy(Population[0])
    best.power = 26
    if ConvergenceFlag == 1:
        res.append(best.power)
    print("开始执行遗传算法+PSO")
    for i in range(g):
        # print("开始第%d次迭代" % (i))
        NewPopulation = []
        for j in range(int(len(Population) / 2)):
            k = random.randint(0, len(Population) - 1)
            X = Population.pop(k)
            k = random.randint(0, len(Population) - 1)
            Y = Population.pop(k)  # 随机选择XY交叉

            n = random.randint(1, int(GeneLength / 10))  # 交叉次数为基因长度除以10,取整
            Children = Crossover2SingleTarget(X, Y, n)  # 其中包括父代两个个体和子代两个个体一共四个[x,y,xc,yc]

            # Children[0].Mutation(6)
            Children[0].MutationSingleTarget(6)
            Children[0].Revise()
            # Children[1].Mutation(6)  #对父代变异
            Children[1].MutationSingleTarget(6)
            Children[1].Revise()

            NewPopulation += Children

        if len(NewPopulation) < m:
            print("错误，当前NewPopulation大小为%d，目标种群大小为%d" % (len(NewPopulation), m))


        for j in range(len(NewPopulation)):
            NewPopulation[j].Revise()  # 修复个体
            NewPopulation[j].CalculateAll()

        Population = ChooseSingleTarget(NewPopulation, m)

        NewPopulation = []
        for j in range(int(len(Population) / 2)):
            k = random.randint(0, len(Population) - 1)
            X = Population.pop(k)
            k = random.randint(0, len(Population) - 1)
            Y = Population.pop(k)  # 随机选择XY交叉

            n = random.randint(1, int(GeneLength / 10))  # 交叉次数为基因长度除以10,取整
            Children = Crossover2SingleTarget(X, Y, n, TargetFlag=1)  # 其中包括父代两个个体和子代两个个体一共四个[x,y,xc,yc]

            # Children[0].Mutation(6)
            Children[0].MutationSingleTarget(6, TargetFlag=1)
            Children[0].Revise()
            # Children[1].Mutation(6)  #对父代变异
            Children[1].MutationSingleTarget(6, TargetFlag=1)
            Children[1].Revise()

            NewPopulation += Children

        for j in range(len(NewPopulation)):
            NewPopulation[j].Revise()  # 修复个体
            NewPopulation[j].CalculateAll()

        Population = ChooseSingleTarget(NewPopulation, m)
        # PopulationZoutendijk(Population)


        for j in range(len(Population)):

            if Population[j].rate > Population[j].LowestRate and Population[i].power > pow and Population[j].power < best.power:
                best = copy.deepcopy(Population[j])
        if ConvergenceFlag == 1 and i%5==0:
            res.append(best.power)

    if ConvergenceFlag == 1:
        return res
    else:
        return best
pow=0.97
def GAPower(indi,m,g):
    '''
    GA寻找最低功率
    '''
    Population = []
    best = copy.deepcopy(indi)
    for i in range(m): #复制某个个体
        newIndi = copy.deepcopy(indi)
        Population.append(newIndi)
    
    GeneLength = len(Population[0].genec[0])
    for i in range(g):
        NewPopulation = []
        for j in range(int(len(Population)/2)):
            k = random.randint(0, len(Population)-1)
            X = Population.pop(k)
            k = random.randint(0, len(Population)-1)
            Y = Population.pop(k)   #随机选择XY交叉 
            
            n = random.randint(1, int(GeneLength/30))#交叉次数为基因长度除以10,取整
            Children = Crossover2SingleTarget(X, Y, n, TargetFlag= 1) #其中包括父代两个个体和子代两个个体一共四个[x,y,xc,yc]
            
            Children[0].MutationOnceSingleTarget(TargetFlag = 1)
            Children[0].Revise
            Children[1].MutationOnceSingleTarget(TargetFlag = 1)  #对父代变异
            Children[1].Revise
    
            NewPopulation += Children
        if len(NewPopulation)<m:
            print("错误，当前NewPopulation大小为%d，目标种群大小为%d"%(len(NewPopulation),m))
        Population = ChooseSingleTarget(NewPopulation, m)
        if Population[0].power<best.power:
            best = copy.deepcopy(Population[0])

    return best
        
def ga_zoutendijk(m,BSList,g, userNumber, LowestRate=1.4,  maxchannal = 3, Qmax = 1, Alpha = 10, B = 20, ConvergenceFlag = 0):
    '''
    单目标优化遗传算法，采用Zoutendijk优化
    '''
    if ConvergenceFlag == 1:
        res = []  #记录收敛性
    Population = []
    for i in range (m):  #产生m个个体形成种群

        indi = Individual(BSList, userNumber,maxchannal, Qmax, Alpha, B)
        #    def __init__(self, BS, userNum, maxc, Qmax, Alpha, B):
        indi.CalculateTotalRate  #计算一下速率，为gener赋值
        indi.Revise()
        indi.CalculateAll()
        indi.LowestRate = LowestRate*userNumber
        Population.append(indi)
    
    GeneLength = len(Population[0].genec[0])
    
    best = copy.deepcopy(Population[0])
    best.power = 26
    if ConvergenceFlag == 1:
        res.append(best.power)
    # print("开始执行遗传算法+Zoutendijk")
    for i in range(g):

        NewPopulation = []
        for j in range(int(len(Population)/2)):    
            k = random.randint(0, len(Population)-1)
            X = Population.pop(k)
            k = random.randint(0, len(Population)-1)
            Y = Population.pop(k)   #随机选择XY交叉 
            
            n = random.randint(1, int(GeneLength/30))#交叉次数为基因长度除以10,取整
            Children = Crossover2SingleTarget(X, Y, n) #其中包括父代两个个体和子代两个个体一共四个[x,y,xc,yc]
            
            #Children[0].Mutation(6)
            Children[0].MutationSingleTarget(6)
            Children[0].Revise()
            #Children[1].Mutation(6)  #对父代变异
            Children[1].MutationSingleTarget(6)
            Children[1].Revise()
    
            NewPopulation += Children
        for j in range(len(NewPopulation)):
            NewPopulation[j].Revise()  #修复个体
            NewPopulation[j].CalculateAll()
            
        if len(NewPopulation)<m:
            print("错误，当前NewPopulation大小为%d，目标种群大小为%d"%(len(NewPopulation),m))
        Population = ChooseSingleTarget(NewPopulation, m)
        NewPopulation = []
        for j in range(int(len(Population) / 2)):
            k = random.randint(0, len(Population) - 1)
            X = Population.pop(k)
            k = random.randint(0, len(Population) - 1)
            Y = Population.pop(k)  # 随机选择XY交叉

            n = random.randint(1, int(GeneLength / 10))  # 交叉次数为基因长度除以10,取整
            Children = Crossover2SingleTarget(X, Y, n, TargetFlag=1)  # 其中包括父代两个个体和子代两个个体一共四个[x,y,xc,yc]

            # Children[0].Mutation(6)
            Children[0].MutationSingleTarget(6, TargetFlag=1)
            Children[0].Revise()
            # Children[1].Mutation(6)  #对父代变异
            Children[1].MutationSingleTarget(6, TargetFlag=1)
            Children[1].Revise()

            NewPopulation += Children

        for j in range(len(NewPopulation)):
            NewPopulation[j].Revise()  #修复个体
            #Population[j].CheckChannel()
            NewPopulation[j].CalculateAll()

        PopulationZoutendijk(NewPopulation)

        
        for j in range(len(NewPopulation)):
            NewPopulation[j].Revise()  #修复个体
            #Population[j].CheckChannel()
            NewPopulation[j].CalculateAll()

        Population = ChooseSingleTarget(NewPopulation, m)
        for j in range(len(Population)):
            
            if Population[j].rate>Population[j].LowestRate and Population[j].power > pow and Population[j].power<best.power:
                best = copy.deepcopy(Population[j])
        if ConvergenceFlag == 1 and i%5==0:
            res.append(best.power)
      
    if ConvergenceFlag == 1:
        return res             
    else:
        return best.power

def muti_ga_zoutendijk(m,BSList,g, userNumber, begin = 1, end = 2, step = 0.1,  maxchannal = 3, Qmax = 1, Alpha = 10, B = 20):
    '''
    多次进行单目标GA，最低速率从begin到end，步长为step
    '''
    result = []
    for LowestRate in range(begin, end+step, step):
        #res = GASingleTarget(m=m, BSList = BSList, g =g, userNumber = userNumber, LowestRate = LowestRate, maxchannal = maxchannal, Qmax = Qmax, Alpha = Alpha, B = B)
        res = ga_zoutendijk(m=m, BSList = BSList, g =g, userNumber = userNumber, LowestRate = LowestRate, maxchannal = maxchannal, Qmax = Qmax, Alpha = Alpha, B = B)
        result.append(res)
    
    return result
    
def ga_convergence(m,BSList,g, userNumber, begin = 1, end = 2, step = 0.1,  maxchannal = 3, Qmax = 1, Alpha = 10, B = 20):
    '''
    单目标优化收敛性验证
    '''
    result = []
    LowestRate = 0.5
    
    res1 = ga(m=m, BSList = BSList, g =g, userNumber = userNumber, LowestRate = LowestRate, maxchannal = maxchannal, Qmax = Qmax, Alpha = Alpha, B = B, ConvergenceFlag=1)
    # res1 = ga_zoutendijk(m=m, BSList = BSList, g =g, userNumber = userNumber, LowestRate = LowestRate, maxchannal = maxchannal, Qmax = Qmax, Alpha = Alpha, B = B, ConvergenceFlag=1)
    # res1 = ga_pso(m=m, BSList = BSList, g =g, userNumber = userNumber, LowestRate = LowestRate, maxchannal = maxchannal, Qmax = Qmax, Alpha = Alpha, B = B, ConvergenceFlag=1)

    result.append(res1)
    return result
