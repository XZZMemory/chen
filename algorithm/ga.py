from classes.Individual import Individual
import copy
import random
from util import RnWLists


def cross_two_point(x, y, n, p=0.95, TargetFlag=0):
    '''
    重写的交叉函数，可以多点交叉

    交叉函数,输入两个个体X,Y，执行交叉，返回clist是两个子代个体
    对每个子基因执行单点交叉，全部执行完毕之后修复个体

    单目标，如果TargetFlag为0交叉信道，否则交叉功率
    '''
    xc = copy.deepcopy(x)  # 深拷贝父代个体x
    yc = copy.deepcopy(y)  # 深拷贝父代个体y
    if random.random() <= 0.95:
        t = x.BSNum  # 获取子基因个数
        l = len(x.genec[0])  # 获取当前子基因长度
        for i in range(t):  # 循环子基因个数次，每个基因进行交叉
            List = [k for k in range(l)]  # 基因位置编号
            slice = random.sample(List, n)  # 随机选取n个
            flag = 0
            for j in range(l):  # 循环子基因个数次
                if j in slice:  # 如果i在选择的序列中
                    flag = 0 if flag == 1 else 1  # 改变交叉标志
                if flag == 1:  # 交叉标志是1的话
                    '''
                    temp = xc.genec[i][j]
                    xc.genec[i][j] = yc.genec[i][j]
                    yc.genec[i][j] = temp

                    temp = xc.genep[i][j]
                    xc.genep[i][j] = yc.genep[i][j]
                    yc.genep[i][j] = temp
                    '''
                    if TargetFlag == 0:  # 单目标优化，如果flag为0，交叉信道，否则交叉功率
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
        if TargetFlag != 0:  # 如果优化的是功率，则修复
            xc.Revise()
            yc.Revise()
    clist = [x, y, xc, yc]  # 两个新生个体
    return clist  # 返回
def choose(Population, m):
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

def ga(m, BSList, g, userNumber, LowestRate=1.4, maxchannal=3, Qmax=1, Alpha=10, B=20, ConvergenceFlag=0):
    '''
    单目标优化遗传算法
    '''
    if ConvergenceFlag == 1:
        res = []  # 记录收敛性
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
    if ConvergenceFlag == 1:
        res.append(best.rp)
    print("开始执行遗传算法")
    for i in range(g):
        # print("开始第%d次迭代"%(i))
        NewPopulation = []
        for j in range(int(len(Population) / 2)):
            k = random.randint(0, len(Population) - 1)
            X = Population.pop(k)
            k = random.randint(0, len(Population) - 1)
            Y = Population.pop(k)  # 随机选择XY交叉

            n = random.randint(1, int(GeneLength / 10))  # 交叉次数为基因长度除以10,取整
            Children = cross_two_point(X, Y, n)  # 其中包括父代两个个体和子代两个个体一共四个[x,y,xc,yc]

            Children[0].MutationSingleTarget(6)
            Children[0].Revise()
            Children[1].MutationSingleTarget(6)  # 对父代变异
            Children[1].Revise()

            NewPopulation += Children

        if len(NewPopulation) < m:
            print("错误，当前NewPopulation大小为%d，目标种群大小为%d" % (len(NewPopulation), m))

        for j in range(len(Population)):
            Population[j].Revise()  # 修复个体
            Population[j].CalculateAll()

        Population = choose(NewPopulation, m)
        for j in range(len(Population)):
            if Population[j].rate > Population[j].LowestRate and Population[j].rp > best.rp:
                best = copy.deepcopy(Population[j])
        if ConvergenceFlag == 1 and i % 5 == 0:
            res.append(best.rp)

    indiSavePath = './data/indi/result_ga.txt'
    print("开始写入个体")
    RnWLists.writeTxt(indiSavePath, best.genec, flag=0)
    RnWLists.writeTxt(indiSavePath, best.genep, flag=1)

    if ConvergenceFlag == 1:
        return res
    else:
        return best.power
