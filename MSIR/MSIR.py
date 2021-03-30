#M.S.I.R Model  黎锦灏   518021910771
import scipy.integrate as spi
import numpy as np
import matplotlib.pyplot as plt

# 传播时间
T = 25

#参数初值
p = 0.9
b = 1
beta = 1.5
d = 0.3
alpha = 0.3
gamma = 0.3
sigma = 0.4
delta = 0.7

S0 = 3.6
I0 = 1.5
R0 = 1.2

# Initial Condition
INI = (S0,I0,R0)


def funcSIR(inivalue,_):
    Y = np.zeros(3)
    X = inivalue
    # Susceptible 易感个体变化
    Y[0] = (1 - p) * b - beta * X[0] * X[1] / (1 + sigma * X[0]) - d * X[0] +delta * X[2]
    # Infection 感染个体变化
    Y[1] = beta * X[0] * X[1] / (1 + sigma * X[0]) - (d + alpha + gamma) *X[1]
    # Recovery 治愈个体变化
    Y[2] = gamma * X[1] + p * b -(d + delta) * X[2]
    return Y

T_range = np.arange(0,T + 1)

RES = spi.odeint(funcSIR,INI,T_range)

plt.plot(RES[:,0],color = 'blue',label = 'Susceptible',marker = '')
plt.plot(RES[:,1],color = 'red',label = 'Infection',marker = '')
plt.plot(RES[:,2],color = 'green',label = 'Recovery',marker = '')
plt.title('ljh2000\'s M.S.I.R Model     2021.03.16')
plt.legend()
plt.xlabel('Time')
plt.ylabel('M.S.I.R')
plt.show()