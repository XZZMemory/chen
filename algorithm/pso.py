# _*_ coding:utf-8 _*_
__author__ = 'seal'
__data__ = '9/29/18'
import random
import copy

# 一个数与矩阵相乘
def get_multi_matrix(num, m):
    h = len(m)
    w = len(m[0])

    for i in range(h):
        for j in range(w):
            m[i][j] *= num

# 两个矩阵相减
def get_diff_matrix(m1, m2):

    h = len(m1)
    w = len(m1[0])

    res = [[0 for i in range(w)] for j in range(h)]

    for i in range(h):
        for j in range(w):
            res[i][j] = m1[i][j] - m2[i][j]
    return res

# 两个矩阵相加
def get_sum_matrix(m1, m2):
    h = len(m1)
    w = len(m1[0])

    for i in range(h):
        for j in range(w):
            m1[i][j] = m1[i][j] + m2[i][j]


# 种群速度更新函数
def update_v_population(population, v_population, w, c1, c2, g_best, z_best):

    for i in range(len(v_population)):
        get_multi_matrix(w, v_population[i])
        sum1 = get_diff_matrix(g_best[i].genep, population[i].genep)
        get_multi_matrix(c1, sum1)
        sum2 = get_diff_matrix(z_best.genep, population[i].genep)
        get_multi_matrix(c2, sum2)
        get_sum_matrix(v_population[i], sum1)
        get_sum_matrix(v_population[i], sum2)


# 种群位置更新函数
def update_population(population, v):
    for i in range(len(population)):
        get_sum_matrix(population[i].genep, v[i])
        population[i].Revise()

def pso(population1, power, time=5, w=0.5, c1=2, c2=2):
    population = copy.deepcopy(population1)
    n = len(population) # 种群中的个体数
    v_init = []
    count = 0
    temp = (1 / power) / 100
    p1 = population[0].genep
    h = len(p1)
    w = len(p1[0])

    # 初始化速度
    for i in range(n):
        v_indi = []
        for j in range(h):
            v_temp = []
            for k in range(w):
                v_temp.append(random.random()*temp)
            v_indi.append(v_temp)
        v_init.append(v_indi)

    # 初始化个体极值位置和群体极值位置
    g_best = copy.deepcopy(population)
    z_best = copy.deepcopy(sorted(population, key=lambda item: item.rp, reverse=True)[0])

    # 更新个体极值、群体极值和速度
    while count < time:
        update_population(population, v_init)

        for i in range(n):
            population[i].CalculateTotalRate()  # 计算一下速率，为gener赋值
            population[i].Revise()
            if(population[i].power == 0 or population[i].power < 0):
                population[i].power = 1
        update_v_population(population, v_init, 0.5, c1, c2, g_best, z_best)

        # 更新个体极值和群体极值
        for i in range(n):
            if g_best[i].rp < population[i].rp:
                g_best[i] = copy.deepcopy(population[i])
            if z_best.rp < population[i].rp:
                z_best = copy.deepcopy(population[i])
        count += 1
