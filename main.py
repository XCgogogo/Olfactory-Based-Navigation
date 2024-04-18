import numpy as np
import math

class map_2D:
    def __init__(self, m, n, initial_belief, sx, sy, sigma_x, sigma_y, beta, l_x, l_y, alpha_ij):
        self.m = m
        self.n = n
        self.M = m * n
        self.initial_belief = 1/self.M
        self.sx, self.sy = sx, sy
        self.sigma_x, self.sigma_y = sigma_x, sigma_y
        self.beta = 0.9
        self.Value = np.zeros(self.M)
    def calculate_f(self, i):
        return (i - 1) % self.m + 1

    def calculate_g(self, i):
        return (i - 1) // self.m + 1

    def calculate_i(self, f, g):
        return f + (g - 1) * self.m

class Agent:

    def __init__(self, env, gamma = 0.9, theta=1e-6):
        self.env = env
        self.gamma = 0.9
        self.theta = 1e-6
        self.action = np.array(
            [[-1, -1], [0, -1], [1, -1], [-1, 0], [1, 0], [-1, 1], [0, 1], [1, 1]]
        )
        self.MAXITER = 1000

    def calculate_Pij(self, i, j):
        # 通过map_2D类的方法获取i和j的f,g坐标
        f_i = self.env.calculate_f(i)
        g_i = self.env.calculate_g(i)
        f_j = self.env.calculate_f(j)
        g_j = self.env.calculate_g(j)

        # 计算x和y方向上的距离
        dx = f_j - f_i
        dy = g_j - g_i
        Pij = (np.exp(-(() ** 2 / (2 * σx_k ** 2) + (dy) ** 2 / (2 * σy_k ** 2))) /
               (2 * np.pi * σx_k * σy_k))


    def value_iter(self):
        for _ in range(self.MAXITER):
            delta = 0
            for i in range(1, self.env.M + 1):
                v = self.env.Value[i - 1]
                for action in self.action:
                    f, g = self.env.calculate_f_g(i)
                    new_f = f + action[0]
                    new_g = g + action[1]
                    if 1 <= new_f <= self.env.m and 1 <= new_g <= self.env.n:
                        new_i = self.env.calculate_i(new_f, new_g)



