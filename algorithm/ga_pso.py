from classes.Individual import Individual
import copy
import random
from util import RnWLists
from algorithm.ga import cross_two_point, choose

def ga_pso(m, BSList, g, userNumber, LowestRate=1.4, maxchannal=3, Qmax=1, Alpha=10, B=20, ConvergenceFlag=0):
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
            Children = cross_two_point(X, Y, n)  # 其中包括父代两个个体和子代两个个体一共四个[x,y,xc,yc]

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

        Population = choose(NewPopulation, m)

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

        Population = choose(NewPopulation, m)
        # PopulationZoutendijk(Population)


        for j in range(len(Population)):

            if Population[j].rate > Population[j].LowestRate and Population[j].power < best.power:
                best = copy.deepcopy(Population[j])
        if ConvergenceFlag == 1 and i%5==0:
            res.append(best.rp)

    indiSavePath = './data/indi/result_ga_pso.txt'
    print("开始写入个体")
    RnWLists.writeTxt(indiSavePath, best.genec, flag=0)
    RnWLists.writeTxt(indiSavePath, best.genep, flag=1)

    if ConvergenceFlag == 1:
        return res
    else:
        return best