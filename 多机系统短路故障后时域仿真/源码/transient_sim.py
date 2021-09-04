import math
import pandas as pd
import numpy as np
from read_data import MyData
from power_flow import PowerFlow
from itertools import combinations
# 绘图部分

# import matplotlib
# matplotlib.use("Qt5Agg")  # 声明使用pyqt5
# from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg  # pyqt5的画布
# import matplotlib.pyplot as plt
# # matplotlib.figure 模块提供了顶层的Artist(图中的所有可见元素都是Artist的子类)，它包含了所有的plot元素
# from matplotlib.figure import Figure


def tand(x):
    return math.tan(math.radians(x))

def sind(x):
    return math.sin(math.radians(x))

def cosd(x):
    return math.cos(math.radians(x))


class TransientSim:
    def __init__(self,data_path):
        # 读入原始数据
        self.read_data(data_path)

    # 读入系统原始数据
    def read_data(self,data_path):
        # 网络结构及参数，包括发电机、支路数据、潮流数据
        self.base_data  = MyData(data_path)
        # 根据基础数据进行潮流计算
        self.power_flow = PowerFlow(self.base_data)
        # 这些数据都不应该被修改
        u_df = pd.read_csv('data_log/general_solution.csv', index_col=0)
        # 各节点电压
        self.u_powerflow = u_df['V'].values
        # 各节点功率
        self.s_powerflow = u_df['S_node'].values
        # 潮流计算用节点导纳矩阵
        self.Y_powerflow = pd.read_csv('data_log/Y_None_None.csv', index_col=0).to_numpy()

        # 读取的字符串转为复数
        for index, item in enumerate(self.u_powerflow):
            self.u_powerflow[index] = complex(item)
        for index, item in enumerate(self.s_powerflow):
            self.s_powerflow[index] = complex(item)
        for index1, item1 in enumerate(self.Y_powerflow):
            for index2, item2 in enumerate(item1):
                self.Y_powerflow[index1][index2] = complex(item2)
        self.node_list =[]
        self.branch_list = []
        for bus_data in self.base_data.bus_data_dict_list:
            node = bus_data['BusNumber']
            self.node_list.append(str(node))
        for branch_data in self.base_data.branch_data_dict_list:
            branch_from = branch_data['TNumber']
            branch_to = branch_data['ZNumber']
            branch_name = str(branch_from)+'_'+str(branch_to)
            self.branch_list.append(branch_name)

    # 设置仿真过程的数据
    def set_simu(self,fault_bus,fault_type,cut_branch,step,cut_time,total_time):
        # 切除时间
        self.cut_time = float(cut_time)
        # 故障类型
        self.fault_type = fault_type
        # 故障节点
        self.fault_bus = int(fault_bus)
        # 仿真总时间
        self.total_time = float(total_time)
        # 仿真步长
        self.step = float(step)
        # 切除支路

        self.cut_branch_head = int(cut_branch.split('_')[0])
        self.cut_branch_end = int(cut_branch.split('_')[1])

        # 额定频率
        self.f_B = 50

        t_in_fault = np.arange(0,self.cut_time+self.step,self.step)
        t_after_fault = np.arange(self.cut_time, self.total_time + self.step, self.step)


        self.t = [t_in_fault,t_after_fault]

    # 启动仿真
    def start_sim(self):
        gen_num = len(self.base_data.gen_data_dict_list)
        omega_B = 2*180*self.f_B

        self.init_sim_data()

        omgea_1 =np.array([omega_B for i in range(gen_num)])
        theta_1 = np.array(self.E0_theta)
        omega_2,theta_2 = self.run_sim(omgea_1,theta_1,self.t[0])

        # 获取故障切除后的转移阻抗
        self.Z, self.alpha = self.get_Z_alpha_cut_fault()
        omega3, theta3 = self.run_sim(omega_2, theta_2, self.t[1])



    # 根据输入的潮流结果及同步电机稳态向量关系计算系统状态量以及代数量在t=0刻的初值
    # 修改系统稳态工况下的导纳矩阵
    def init_sim_data(self):
        bus_num = len(self.base_data.bus_data_dict_list)
        gen_num = len(self.base_data.gen_data_dict_list)



        # 应该是一个n*t的矩阵，n是发电机数目，t是仿真的全部步长
        self.P = []
        self.omega =[]
        self.theta =[]
        # 获取故障下的转移阻抗
        self.Z ,self.alpha = self.get_Z_alpha_fault()
        # 获取发电机初态
        self.E0_amp,self.E0_theta = self.get_E0()
        self.P0 = self.get_P(self.E0_theta)
        self.PT = self.get_PT()

        self.Tj = np.zeros([gen_num])
        for gen_index,gen_data in enumerate(self.base_data.gen_data_dict_list):
            self.Tj[gen_index] =  gen_data['Tj']
        # 开始仿真，时间指针归零
        self.current_index = 0


    # 根据当前的相角获得当前功率
    def get_P(self,E_theta):
        gen_num = len(self.base_data.gen_data_dict_list)
        # 功率初值是浮点数
        P = np.zeros([gen_num],dtype=float)
        for i in range(gen_num):
            P_temp =  0
            for j in range(gen_num):
                theta_i_j = E_theta[i]-E_theta[j]
                aplha_i_j = self.alpha[i][j]
                # 幅值是不变的
                Egi = self.E0_amp[i]
                Egj = self.E0_amp[j]
                if i ==j:
                    aplha_i_j = - aplha_i_j
                P_temp += Egi * Egj * sind(theta_i_j - aplha_i_j) / self.Z[i][j]
            P[i] = P_temp
        return P


    # 得到原动机功率，默认就是稳态下发电机的输出功率
    def get_PT(self):
        gen_num = len(self.base_data.gen_data_dict_list)
        PT = np.zeros([gen_num],dtype=float)

        for index,gen_data in enumerate(self.base_data.gen_data_dict_list):
            # 电机连接的节点
            bus_index = gen_data['BUS']-1
            bus_data=None

            for bus_data_dict in self.base_data.bus_data_dict_list:
                # 获取当前节点数据
                if bus_data_dict['BusNumber']==bus_index+1:
                    bus_data = bus_data_dict

            loadP= bus_data['LoadP']
            loadQ= bus_data['LoadQ']


            # 获取发电机注入节点的功率，也就是原动机输入的功率
            S = (self.s_powerflow[bus_index]/100+complex(loadP,loadQ))
            PT[index]=S.real
        return PT

    def get_E0(self):
        bus_num = len(self.base_data.bus_data_dict_list)
        gen_num = len(self.base_data.gen_data_dict_list)
        E0 = np.zeros([gen_num],dtype=complex)
        E0_amp = np.zeros([gen_num], dtype=float)
        E0_theta = np.zeros([gen_num], dtype=float)
        for index,gen_data in enumerate(self.base_data.gen_data_dict_list):
            # 电机连接的节点
            bus_index = gen_data['BUS']-1
            bus_data=None

            for bus_data_dict in self.base_data.bus_data_dict_list:
                # 获取当前节点数据
                if bus_data_dict['BusNumber']==bus_index+1:
                    bus_data = bus_data_dict

            loadP= bus_data['LoadP']
            loadQ= bus_data['LoadQ']

            # 获取发电机连接的节点电压，反推发电机电势
            U = self.u_powerflow[bus_index]
            # 获取发电机注入节点的功率
            S = (self.s_powerflow[bus_index]/100+complex(loadP,loadQ))
            Xd1 = gen_data['Xd1']
            Xt = gen_data['Xt']
            X1 = Xd1+Xt

            # 计算该节点注入电流
            I = (S/U).conjugate()
            # 用末端电压和功率计算首端电压
            E0[index] = U+I*complex(0,X1)
            E0_amp[index] = abs(E0[index])
            E0_theta[index] = math.atan(E0[index].imag/E0[index].real)*180/3.14
        # 测试用数据，和题目里算出来的一样
        E0_amp =[1.52532,1.26143,1.12089]
        E0_theta =[38.64149,29.44408,15.52411]
        return E0_amp,E0_theta


    # 获取故障下的只含电机节点的转移阻抗矩阵以及角度，用于初态和后续暂态仿真的计算
    def get_Z_alpha_fault(self):
        bus_num = len(self.base_data.bus_data_dict_list)
        gen_num = len(self.base_data.gen_data_dict_list)

        # 初始化为潮流计算用的节点导纳矩阵

        Y1 = self.Y_powerflow

        zero_row = np.zeros([bus_num], dtype=complex)
        # 增加矩阵的维度,维度数为电机数
        for i in range(0, gen_num):
            Y1 = np.row_stack((Y1, zero_row))
        zero_cloumn = np.zeros([gen_num + bus_num], dtype=complex)
        for i in range(0, gen_num):
            Y1 = np.column_stack((Y1, zero_cloumn))

        # 生成端口-节点关联矩阵
        Ml = []
        for index, gen_data in enumerate(self.base_data.gen_data_dict_list):
            Ml_temp = np.zeros(bus_num + gen_num)
            # 发电机
            Ml_temp[index + bus_num] = 1
            Ml.append(Ml_temp)
        Ml = np.array(Ml)

        # 修改导纳矩阵
        # 修改矩阵,将暂态电抗并入导纳方程，发电机节点电压即为暂态电压，并且幅值不变，相角改变
        for index, gen_data in enumerate(self.base_data.gen_data_dict_list):
            bus_index = gen_data['BUS'] - 1

            Xd1 = gen_data['Xd1']
            Xt = gen_data['Xt']
            X1 = complex(0, Xd1 + Xt)

            # 修改自阻抗
            Y1[bus_num + index][bus_num + index] += 1 / X1
            # 修改连接母线节点的自阻抗
            Y1[bus_index][bus_index] += 1 / X1
            # 修改所连接的母线节点的自阻抗
            Y1[bus_index][bus_num + index] -= 1 / X1
            Y1[bus_num + index][bus_index] -= 1 / X1

        # 修改由负荷影响的矩阵
        for index, bus_data in enumerate(self.base_data.bus_data_dict_list):
            bus_index = bus_data['BusNumber'] - 1
            loadP = bus_data['LoadP']
            loadQ = bus_data['LoadQ']
            loadS = complex(loadP, loadQ)
            U = self.u_powerflow[bus_index]

            # 由于负荷接地，修改自导纳即可
            if loadP != 0 or loadQ != 0:
                Z_ld = loadS * (abs(U) ** 2) / (abs(loadS) ** 2)
                Y1[bus_index][bus_index] += 1 / Z_ld

        # 修改由附加电抗影响的矩阵部分

        # 计算附加阻抗
        self.zff2 = self.get_Zff2()
        self.zff0 = self.get_Zff0()

        self.z_delta = complex(0, 0)
        if self.fault_type == '单相接地':
            self.z_delta = self.zff0 + self.zff2
        elif self.fault_type  == '两相短路':
            self.z_delta = self.zff2
        elif self.fault_type  == '两相接地':
            self.z_delta = (self.zff0 * self.zff2) / (self.zff0 + self.zff2)
        elif self.fault_type == '三相短路':
            self.z_delta = 0

        # 只需修改故障节点的自导纳
        if self.z_delta != 0:
            y_delta = 1 / self.z_delta
            Y1[self.fault_bus - 1][self.fault_bus - 1] += y_delta

        Ml = Ml.T
        Y1 = Y1.astype(np.complex)

        Z1 = np.linalg.inv(Y1)


        Zeq = np.matmul(Ml.T, np.matmul(Z1, Ml))
        # 以电机连接节点为端口的戴维南等效电阻
        Zeq = Zeq.astype(np.complex)
        # 以电机连接节点为端口的诺顿等效电路
        Yeq = np.linalg.inv(Zeq)

        # 转移阻抗矩阵
        Z = np.zeros([gen_num, gen_num], dtype=float)
        alpha = np.zeros([gen_num, gen_num], dtype=float)
        for index1, item1 in enumerate(Yeq):
            for index2, item2 in enumerate(item1):
                # 转移阻抗是节点导纳矩阵的倒数，而不是求逆
                z = 1 / item2

                Z[index1][index2] = abs(z)
                # 角度注意是弧度角！！！
                alpha_temp = math.atan(z.imag / z.real) * 180 / 3.14
                # 转化到1，2象限
                if alpha_temp < 0:
                    alpha_temp = alpha_temp + 180
                alpha[index1][index2] = 90-alpha_temp

        return Z, alpha

    # 获取故障切除后的只含电机节点的转移阻抗矩阵以及角度，用于后续暂态仿真的计算
    def get_Z_alpha_cut_fault(self):
        bus_num = len(self.base_data.bus_data_dict_list)
        gen_num = len(self.base_data.gen_data_dict_list)

        # 初始化为潮流计算用的节点导纳矩阵

        Y1 = self.Y_powerflow

        zero_row = np.zeros([bus_num], dtype=complex)
        # 增加矩阵的维度,维度数为电机数
        for i in range(0, gen_num):
            Y1 = np.row_stack((Y1, zero_row))
        zero_cloumn = np.zeros([gen_num + bus_num], dtype=complex)
        for i in range(0, gen_num):
            Y1 = np.column_stack((Y1, zero_cloumn))

        # 生成端口-节点关联矩阵
        Ml = []
        for index, gen_data in enumerate(self.base_data.gen_data_dict_list):
            # bus_index = gen_data['BUS']-1
            Ml_temp = np.zeros(bus_num + gen_num)
            # 发电机
            Ml_temp[index + bus_num] = 1
            Ml.append(Ml_temp)
        Ml = np.array(Ml)

        # 修改导纳矩阵
        # 修改矩阵,将暂态电抗并入导纳方程，发电机节点电压即为暂态电压，并且幅值不变，相角改变
        for index, gen_data in enumerate(self.base_data.gen_data_dict_list):
            bus_index = gen_data['BUS'] - 1

            Xd1 = gen_data['Xd1']
            Xt = gen_data['Xt']
            X1 = complex(0, Xd1 + Xt)

            # 修改自阻抗
            Y1[bus_num + index][bus_num + index] += 1 / X1
            # 修改连接母线节点的自阻抗
            Y1[bus_index][bus_index] += 1 / X1
            # 修改所连接的母线节点的自阻抗
            Y1[bus_index][bus_num + index] -= 1 / X1
            Y1[bus_num + index][bus_index] -= 1 / X1

        # 修改由负荷影响的矩阵
        for index, bus_data in enumerate(self.base_data.bus_data_dict_list):
            bus_index = bus_data['BusNumber'] - 1
            loadP = bus_data['LoadP']
            loadQ = bus_data['LoadQ']
            loadS = complex(loadP, loadQ)
            U = self.u_powerflow[bus_index]

            # 由于负荷接地，修改自导纳即可
            if loadP != 0 or loadQ != 0:
                Z_ld = loadS * (abs(U) ** 2) / (abs(loadS) ** 2)
                Y1[bus_index][bus_index] += 1 / Z_ld


        # 修改线路切除导致阻抗增加的部分
        self.cut_branch_head = 1
        self.cut_branch_end = 2
        z_branch = None
        for branch_data in (self.base_data.branch_data_dict_list):
            if branch_data['TNumber']==self.cut_branch_head and branch_data['ZNumber']==self.cut_branch_end:
                r = branch_data['R1']
                x = branch_data['X1']
                z_branch = complex(r,x)

        # 移除原支路
        Y1[self.cut_branch_head-1][self.cut_branch_head-1] -= 1 / z_branch
        Y1[self.cut_branch_end - 1][self.cut_branch_end - 1] -= 1 / z_branch
        Y1[self.cut_branch_head - 1][self.cut_branch_end - 1] += 1 / z_branch
        Y1[self.cut_branch_end - 1][self.cut_branch_head - 1] += 1 / z_branch

        # 增加现支路
        Y1[self.cut_branch_head-1][self.cut_branch_head-1] += 1 / (z_branch*2)
        Y1[self.cut_branch_end - 1][self.cut_branch_end - 1] += 1 / (z_branch*2)
        Y1[self.cut_branch_head - 1][self.cut_branch_end - 1] -= 1 / (z_branch*2)
        Y1[self.cut_branch_end - 1][self.cut_branch_head - 1] -= 1 / (z_branch*2)


        Ml = Ml.T
        Y1 = Y1.astype(np.complex)

        Z1 = np.linalg.inv(Y1)


        Zeq = np.matmul(Ml.T, np.matmul(Z1, Ml))
        # 以电机连接节点为端口的戴维南等效电阻
        Zeq = Zeq.astype(np.complex)
        # 以电机连接节点为端口的诺顿等效电路
        Yeq = np.linalg.inv(Zeq)

        # 转移阻抗矩阵
        Z = np.zeros([gen_num, gen_num], dtype=float)
        alpha = np.zeros([gen_num, gen_num], dtype=float)
        for index1, item1 in enumerate(Yeq):
            for index2, item2 in enumerate(item1):
                # 转移阻抗是节点导纳矩阵的倒数，而不是求逆
                z = 1 / item2

                Z[index1][index2] = abs(z)
                # 角度注意是弧度角！！！
                alpha_temp = math.atan(z.imag / z.real) * 180 / 3.14
                # 转化到1，2象限
                if alpha_temp < 0:
                    alpha_temp = alpha_temp + 180
                alpha[index1][index2] = 90-alpha_temp

        return Z, alpha
    # 故障节点为1，需要获得负序看进去的戴维南，从而获得附加电抗
    def get_Zff2(self):
        # 母线数量
        bus_num = len(self.base_data.bus_data_dict_list)

        # 关联向量,长度为节点数
        Ml = np.zeros(bus_num)
        # 故障节点为1
        Ml[self.fault_bus-1] = 1

        # 形成负序节点导纳矩阵,节点数目仍为母线数目，因为发电机接地了
        Y_2 =np.zeros([bus_num,bus_num],dtype=complex)
        # 可以作为一个单独的程序,参数为电路的序，待完善
        for bus_data_dict in self.base_data.bus_data_dict_list:
            # 获取当前节点
            node_num =bus_data_dict['BusNumber']

            # 遍历该节点连接的母线
            for branch_data_dict in self.base_data.branch_data_dict_list:
                if  branch_data_dict['TNumber']==(node_num):
                    to_node_num = branch_data_dict['ZNumber']

                    # 考虑线路的阻抗对矩阵的影响
                    node_data_z = complex(branch_data_dict['R2'],
                                          branch_data_dict['X2'])
                    node_data_y = 1 / node_data_z

                    # 修改节点自导纳
                    Y_2[node_num-1, node_num-1] +=node_data_y
                    # 修改连接节点自导纳
                    Y_2[to_node_num-1, to_node_num-1] += node_data_y
                    # 修改互导纳
                    Y_2[node_num-1, to_node_num-1] -= node_data_y
                    Y_2[to_node_num-1, node_num-1] -= node_data_y




        # 修改导纳矩阵
        # 修改发电机阻抗对节点的影响
        for index,gen_data in enumerate(self.base_data.gen_data_dict_list):
            bus_index = gen_data['BUS']-1
            # Xd1 = gen_data['Xd1']
            X2 = gen_data['X2']
            Xt = gen_data['Xt']
            Xg2 = complex(0,X2+Xt)
            # 修改连接母线节点的自阻抗
            Y_2[bus_index][bus_index] += 1 / Xg2


        # 修改由负荷影响的矩阵
        for index,bus_data in enumerate(self.base_data.bus_data_dict_list):
            bus_index = bus_data['BusNumber']-1
            loadP = bus_data['LoadP']
            loadQ = bus_data['LoadQ']
            loadS = complex(loadP,loadQ)
            U = self.u_powerflow[bus_index]

            # 由于负荷接地，修改自导纳即可
            if loadP!=0 or loadQ!=0:
                Z_ld = loadS*(abs(U)**2)/(abs(loadS)**2)
                z = complex(0,0.35*abs(Z_ld))
                Y_2[bus_index][bus_index] += 1/z

        Y_2 = Y_2.astype(np.complex)
        Z_2 = np.linalg.inv(Y_2)
        Zeq = np.matmul(Ml.T, np.matmul(Z_2, Ml))
        return Zeq
    # 故障节点为1，需要获得零序看进去的戴维南，从而获得附加电抗
    def get_Zff0(self):
        # 母线数量
        bus_num = len(self.base_data.bus_data_dict_list)

        # 关联向量,长度为节点数
        Ml = np.zeros(bus_num)
        # 故障节点支为1
        Ml[self.fault_bus-1] = 1

        # 形成负序节点导纳矩阵,节点数目仍为母线数目，因为发电机接地了
        Y_0 =np.zeros([bus_num,bus_num],dtype=complex)
        # 可以作为一个单独的程序,参数为电路的序，晚点考虑
        for bus_data_dict in self.base_data.bus_data_dict_list:
            # 获取当前节点
            node_num =bus_data_dict['BusNumber']

            # 遍历该节点连接的母线
            for branch_data_dict in self.base_data.branch_data_dict_list:
                if  branch_data_dict['TNumber']==(node_num):
                    to_node_num = branch_data_dict['ZNumber']

                    # 考虑线路的阻抗对矩阵的影响
                    node_data_z = complex(branch_data_dict['R0'],
                                          branch_data_dict['X0'])
                    node_data_y = 1 / node_data_z

                    # 修改节点自导纳
                    Y_0[node_num-1, node_num-1] +=node_data_y
                    # 修改连接节点自导纳
                    Y_0[to_node_num-1, to_node_num-1] += node_data_y
                    # 修改互导纳
                    Y_0[node_num-1, to_node_num-1] -= node_data_y
                    Y_0[to_node_num-1, node_num-1] -= node_data_y

        # 修改导纳矩阵
        # 修改发电机阻抗对节点的影响
        for index,gen_data in enumerate(self.base_data.gen_data_dict_list):
            bus_index = gen_data['BUS']-1
            # 零序电路发电机三角形接法里不流通，只有变压器阻抗
            Xt = gen_data['Xt']
            Xg0 = complex(0,Xt)
            # 修改连接母线节点的自阻抗
            Y_0[bus_index][bus_index] += 1 / Xg0

        # 异步电动机以及多数负荷常常接成三角形，或者接成不接地星形
        # 零序电流不能流通，故不需要建立零序等值
        Y_0 = Y_0.astype(np.complex)
        Z_0 = np.linalg.inv(Y_0)
        Zeq = np.matmul(Ml.T, np.matmul(Z_0, Ml))
        return Zeq

    # 开始仿真,仿真只和初值和仿真时段有关
    def run_sim(self,omega_init,theta_init,t):

        gen_num = len(self.base_data.gen_data_dict_list)
        omega_B = 2*180*self.f_B

        # 应该是一个n*t的矩阵，n是发电机数目，t是仿真的全部步长
        P = [np.zeros(gen_num) for i in t]
        omega = [np.zeros(gen_num) for i in t]
        theta =[np.zeros(gen_num) for i in t]

        # 状态量微分值
        domega_dt = [np.zeros(gen_num) for i in t]
        dtheta_dt = [np.zeros(gen_num) for i in t]


        # 微分方程初值
        # 初始角度
        theta[0] = theta_init
        # 初始转速都为额定转速
        omega[0] = omega_init

        # 开始计算,从t=1开始
        for t_index,t_step in enumerate(t[:-1]) :
            # 用第i步状态量的初值计算代数量
            P[t_index] = self.get_P(theta[t_index])

            # dw/dt_n = Pm - Pen - D(omega -1)
            # D为常阻尼系数，需要计及机械阻尼时才有
            for gen_index in range(gen_num):
                # 前半步,计算初始时间段变化速度
                domega_dt[t_index][gen_index] = ((self.PT[gen_index] - P[t_index][gen_index])/self.Tj[gen_index])*180/3.14
                dtheta_dt[t_index][gen_index]=omega[t_index][gen_index]-omega_B

                # 求得时间末段近似值
                omega[t_index+1][gen_index] = omega[t_index][gen_index]+ domega_dt[t_index][gen_index]*self.step
                theta[t_index+1][gen_index] = theta[t_index][gen_index]+ dtheta_dt[t_index][gen_index]*self.step

            # 根据i+1的预测状态量预估i+1时刻的代数量
            P[t_index+1]= self.get_P(theta[t_index+1])

            # 计算时间末端变化速度
            for gen_index in range(gen_num):

                domega_dt[t_index+1][gen_index] = ((self.PT[gen_index] - P[t_index+1][gen_index])/self.Tj[gen_index])*180/3.14
                dtheta_dt[t_index+1][gen_index]= omega[t_index+1][gen_index]-omega_B

            # 校正数据
            for gen_index in range(gen_num):
                omega[t_index+1][gen_index] = omega[t_index][gen_index]+(domega_dt[t_index+1][gen_index]+domega_dt[t_index][gen_index])*self.step/2
                theta[t_index+1][gen_index] = theta[t_index][gen_index]+(dtheta_dt[t_index+1][gen_index]+dtheta_dt[t_index][gen_index])*self.step/2

        self.omega.append(omega)
        self.theta.append(theta)
        self.P.append(P)
        return omega[-1],theta[-1]
    def process_result(self,para_name):
        gen_num = len(self.base_data.gen_data_dict_list)
        t =[]
        for para_list in self.t:
            t.extend(para_list)
        t = np.array(t)
        t = t.reshape((-1,1))
        result_list = []
        index_list =[]
        if para_name == 'omega':
            omega = []
            for para_list in self.omega:
                omega.extend(para_list)
            theta = np.array(omega)
            for i in range(gen_num):
                temp_result = theta[:,i]
                temp_result = temp_result.reshape((-1,1))
                result_list.append(temp_result)
                index_list.append(i+1)

        if para_name == 'theta_abs':
            theta = []
            for para_list in self.theta:
                theta.extend(para_list)
            theta = np.array(theta)
            for i in range(gen_num):
                temp_result = theta[:,i]
                temp_result = temp_result.reshape((-1,1))
                result_list.append(temp_result)
                index_list.append(i + 1)

        if para_name == 'theta_opp':
            theta = []
            for para_list in self.theta:
                theta.extend(para_list)
            theta = np.array(theta)

            # 用组合的方式算出有哪些相对角
            opp_list =  list(combinations(range(gen_num), 2))
            for i_j in opp_list:
                i,j =i_j
                temp_result_i = theta[:, i]
                temp_result_j = theta[:, j]
                temp_result = temp_result_i-temp_result_j
                result_list.append(temp_result)
                index_list.append('{}{}'.format(str(i+1),str(j+1)))

        if para_name == 'P':
            P = []
            for para_list in self.P:
                P.extend(para_list)
            P = np.array(P)
            for i in range(gen_num):
                temp_result = P[:,i]
                temp_result = temp_result.reshape((-1,1))
                result_list.append(temp_result)
                index_list.append(i + 1)

        return t,result_list,index_list


    def run_sim_2(self):
        gen_num = len(self.base_data.gen_data_dict_list)
        # 要修改阻抗矩阵
        self.Z, self.alpha = self.get_Z_alpha_cut_fault()
        # 额定转速，单位为角度
        omega_B = 2*180*self.f_B
        # 时间置零
        i = 0
        step = 0.01
        # 故障发生时段
        t_staus1_range =[0.0,10]
        #均分，因为包括端点，所以数目应该是仿真步长+1
        self.t = np.arange(t_staus1_range[0],t_staus1_range[1]+step,step)

        t_len = len(self.t)


        # 应该是一个n*t的矩阵，n是发电机数目，t是仿真的全部步长
        P = [np.zeros(gen_num) for i in self.t]
        omega = [np.zeros(gen_num) for i in self.t]
        theta =[np.zeros(gen_num) for i in self.t]

        # 状态量微分值
        domega_dt = [np.zeros(gen_num) for i in self.t]
        dtheta_dt = [np.zeros(gen_num) for i in self.t]

        # 微分方程初值
        # 初始角度
        theta[0] =  self.theta
        # 初始转速都为额定转速
        omega[0] =  self.omega

        # 开始计算,从t=1开始
        for i in range(t_len)[:-1]:
            # 用第i步状态量的初值计算代数量
            P[i] = self.get_P(theta[i])

            # dw/dt_n = Pm - Pen - D(omega -1)
            # D为常阻尼系数，需要计及机械阻尼时才有
            for gen_index in range(gen_num):

                # 前半步,计算初始时间段变化速度
                domega_dt[i][gen_index] = ((self.PT[gen_index] - P[i][gen_index])/self.Tj[gen_index])*180/3.14
                dtheta_dt[i][gen_index]= omega[i][gen_index]-omega_B

                # 求得时间末段近似值
                omega[i+1][gen_index] = omega[i][gen_index]+ domega_dt[i][gen_index]*step
                theta[i+1][gen_index] = theta[i][gen_index]+ dtheta_dt[i][gen_index]*step

            # 根据i+1的预测状态量预估i+1时刻的代数量
            P[i+1]= self.get_P(theta[i+1])

            # 计算时间末端变化速度
            for gen_index in range(gen_num):
                domega_dt[i+1][gen_index] = ((self.PT[gen_index] - P[i+1][gen_index])/self.Tj[gen_index])*180/3.14
                dtheta_dt[i+1][gen_index]= omega[i+1][gen_index]-omega_B

            # 校正数据
            for gen_index in range(gen_num):
                omega[i+1][gen_index] = omega[i][gen_index]+(domega_dt[i+1][gen_index]+domega_dt[i][gen_index])*step/2
                theta[i+1][gen_index] = theta[i][gen_index]+(dtheta_dt[i+1][gen_index]+dtheta_dt[i][gen_index])*step/2

        theta_1 =[]
        theta_2 = []
        theta_3 = []

        for para_list in theta:
            theta_1.append(para_list[0])
            theta_2.append(para_list[1])
            theta_3.append(para_list[2])

        # 设置绘图参数

        plt.rcParams['font.sans-serif'] = ['SimHei']
        plt.rcParams['axes.unicode_minus'] = False

        font1 = {'family': 'SimHei',
                 'weight': 'normal',
                 'size': 15,
                 }

        plt.plot(self.t,theta_1,label='theta_1')
        plt.plot(self.t,theta_2,label='theta_2')
        plt.plot(self.t,theta_3,label='theta_3')
        plt.legend()
        plt.show()

        theta_12 = []
        theta_13 = []
        theta_23 = []
        for para_list in theta:
            theta_12.append(para_list[0]-para_list[1])
            theta_13.append(para_list[0]-para_list[2])
            theta_23.append(para_list[1]-para_list[2])

        plt.plot(self.t,theta_12,label='theta_12')
        plt.plot(self.t,theta_13,label='theta_13')
        plt.plot(self.t,theta_23,label='theta_23')
        plt.legend()
        plt.show()





if __name__ == '__main__':
     sim =TransientSim('data.txt')