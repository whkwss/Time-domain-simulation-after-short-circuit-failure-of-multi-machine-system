import numpy as np
import matplotlib.pyplot as plt

class Euler_Slove:
    # 初始化状态
    def __init__(self,x0,t_range,n,function_str):
        # 范围
        self.t_range = t_range
        # 步长
        self.n = n
        # 初值
        self.x0 = x0
        # 函数值2
        self.function_str = function_str
        # 用于储存微分方程的解的列表
        self.dx_dt_list = []
        self.x_list = []
        # 求解微分方程
        self.slove()

    # 可以通过这里修改函数来修改
    def function(self,x):
        return eval(self.function_str)
    def slove(self):

        step = 1 / self.n
        t = np.arange(self.t_range[0],self.t_range[1],step)
        dx_dt_list=[]
        x_list=[]
        # 带入初值
        x = self.x0
        dx_dt = self.function(x)

        # 开始循环
        for i in t:
            # 记录迭代前的值
            x_list.append(x)
            dx_dt_list.append(dx_dt)
            # 更新斜率
            dx_dt = self.function(x)
            # 计算横坐标变换值
            x_delta = dx_dt*step
            # 曲线上下一个点
            x = x+x_delta
        self.dx_dt_list = dx_dt_list
        self.x_list = x_list




class Modified_Euler_Slove:
    # 初始化状态
    def __init__(self,x0,t_range,n,function_str):
        # 范围
        self.t_range = t_range

        # 步长
        self.n = n

        # 初值
        self.x0 = x0

        # 函数值2
        self.function_str = function_str

        # 微分方程的解
        self.dx_dt_list = []
        self.x_list = []

        self.slove()

    # 可以通过这里修改函数来修改
    def function(self,x):
        return eval(self.function_str)
    def slove(self):

        step = 1 / self.n
        t = np.arange(self.t_range[0],self.t_range[1],step)
        dx_dt_list=[]
        x_list=[]
        # 带入初值
        x = self.x0
        dx_dt = self.function(x)

        # 开始循环
        for i in t:
            # 记录迭代前的值值
            x_list.append(x)
            dx_dt_list.append(dx_dt)
            # 更新斜率
            dx_dt_0 = self.function(x)
            # 按普通欧拉法计算横坐标变换值
            x_0= x+dx_dt*step
            # 计算欧拉法得到修正后的函数值
            dx_dt = self.function(x_0)
            # 修正后的x
            x = x + 1/2*(dx_dt_0+dx_dt)*step
        self.dx_dt_list = dx_dt_list
        self.x_list = x_list


# 只针对 dx/dt = λx的情况，因此传入的参数不一样
class Implicit_Trapezoid_Method:
    # 初始化状态
    def __init__(self,x0,t_range,n,lamuda):
        self.lamuda = lamuda
        # 范围
        self.t_range = t_range
        # 步长
        self.n = n
        # 初值
        self.x0 = x0
        # 函数值
        # 微分方程的解
        self.dx_dt_list = []
        self.x_list = []
        self.slove()

    # 可以通过这里修改函数来修改
    def function(self,x):
        return self.lamuda*x
    def slove(self):

        step = 1 / self.n
        t = np.arange(self.t_range[0],self.t_range[1],step)
        dx_dt_list=[]
        x_list=[]
        # 带入初值
        x = self.x0
        dx_dt = self.function(x)

        # 开始循环
        for i in t:
            # 记录迭代前的值值
            x_list.append(x)
            dx_dt_list.append(dx_dt)
            # 修正并求解方程后得到的x
            x = x*(1+self.lamuda*step/2)/(1-self.lamuda*step/2)
            # 使用修正后的x计算
            dx_dt = self.function(x)
        self.dx_dt_list = dx_dt_list
        self.x_list = x_list


if __name__ == '__main__':
    x0 =1
    t_range=[0,1]
    n=10
    function_str = '-20*x'
    t = np.arange(t_range[0],t_range[1],1/n)

    x1 = Euler_Slove(x0,t_range,n,function_str).x_list
    x2 = Modified_Euler_Slove(x0,t_range,n,function_str).x_list
    x3 = Implicit_Trapezoid_Method(x0,t_range,n,-20).x_list



    # 设置绘图参数
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False

    font1 = {'family': 'SimHei',
             'weight': 'normal',
             'size': 15,
             }

    plt.figure(figsize=(12, 4))
    plt.subplot(1, 3, 1)
    plt.plot(t, x1, linewidth=3, color='r', label='欧拉法求解')
    plt.xlabel('t', fontsize=24)
    plt.ylabel('x', fontsize=24)
    plt.legend(prop=font1)
    plt.grid()

    plt.subplot(1, 3, 2)
    plt.plot(t, x2, linewidth=3, color='b', label='改进欧拉法求解')
    plt.xlabel('t', fontsize=24)
    plt.ylabel('x', fontsize=24)
    plt.legend(prop=font1)
    plt.grid()

    plt.subplot(1, 3, 3)
    plt.plot(t, x3, linewidth=3, color='g', label='隐式梯形法求解')
    plt.xlabel('t', fontsize=24)
    plt.ylabel('x', fontsize=24)
    plt.legend(prop=font1)
    plt.grid()

    plt.show()