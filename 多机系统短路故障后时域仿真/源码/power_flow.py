# sss=1
# 求正负零序的诺顿等值，包括故障前，故障后
# P134
# 线路数据

import math
import numpy as np
import pandas as pd
from read_data import MyData
class PowerFlow:
    def __init__(self,base_data):
        # 初始化
        self.status_init(base_data)
        # 计算节点导纳矩阵
        self.generate_node_admat()
        self.resort_my_admat()
        # 潮流解算
        self.cal_power_flow()
        # 输出结果
        self.output_result('general')
    #  状态初始化
    def status_init(self,base_data):
        # 误差精度
        self.error_tol =0.00001
        # 默认循环状态开启
        self.circle_status=True
        # 记录循环次数
        self.circle_count=0
        # 默认收敛情况
        self.coverge_status = False

        # n-1切除情况
        self.cut_node_1 =None
        self.cut_node_2 =None

        # 读取数据
        self.bus_data_dict_list = base_data.bus_data_dict_list
        self.branch_data_dict_list = base_data.branch_data_dict_list
        self.gen_data_dict_list = base_data.gen_data_dict_list

        # 雅可比矩阵计算部分
        # 生成nxn大小的矩阵，n为节点数
        self.Y=np.zeros([len(self.bus_data_dict_list), len(self.bus_data_dict_list)],dtype=complex)
        # 生成雅可比矩阵,默认只有1个平衡节点,其余全是PQ+PV
        self.Jacobbi=np.zeros([(len(self.bus_data_dict_list)-1)*2, (len(self.bus_data_dict_list)-1)*2])
        # 电压修正量
        self.e=np.zeros([len(self.bus_data_dict_list)])
        self.f = np.zeros([len(self.bus_data_dict_list)] )
        # 功率修正量
        self.PQUs = np.zeros([(len(self.bus_data_dict_list)-1)*2])
        self.PQU =  np.zeros([(len(self.bus_data_dict_list)-1) * 2])
        self.dPQU=  np.zeros([(len(self.bus_data_dict_list)-1) * 2])

        # 节点功率
        self.S_node = np.zeros([len(self.bus_data_dict_list)], dtype=complex)
        # 节点电流
        self.I = np.zeros([len(self.bus_data_dict_list), len(self.bus_data_dict_list)], dtype=complex)
        # 电压最终值直角坐标表达
        self.u=[]
        # 电压最终值极坐标表达
        self.u_amp = []
        self.u_angle = []
        # 节点功率最终值
        self.S =[]




    # 生成节点导纳矩阵
    def generate_node_admat(self):

        for bus_data_dict in self.bus_data_dict_list:
            # 获取当前节点
            node_num =bus_data_dict['BusNumber']
            # 遍历该节点连接的母线
            for branch_data_dict in self.branch_data_dict_list:
                if  branch_data_dict['TNumber']==(node_num):
                    to_node_num = branch_data_dict['ZNumber']

                    # 考虑线路的阻抗对矩阵的影响

                    node_data_z = complex(branch_data_dict['R1'],
                                          branch_data_dict['X1'])
                    node_data_y = 1 / node_data_z

                    # 修改节点自导纳
                    self.Y[node_num-1, node_num-1] +=node_data_y
                    # 修改连接节点自导纳
                    self.Y[to_node_num-1, to_node_num-1] += node_data_y
                    # 修改互导纳
                    self.Y[node_num-1, to_node_num-1] -= node_data_y
                    self.Y[to_node_num-1, node_num-1] -= node_data_y

    # 节点导纳矩阵核对无误
    def resort_my_admat(self):
        # 按节点类型重排
        self.bus_data_dict_list = sorted(self.bus_data_dict_list, key=lambda x: x["Type"])
        self.my_sort_list = [data_dict['BusNumber'] for data_dict in self.bus_data_dict_list]
        my_sort_list= [i-1 for i in self.my_sort_list]

        # 重排Y
        mat = pd.DataFrame(self.Y)
        mat.to_csv('data_log/Y_{}_{}.csv'.format(self.cut_node_1,self.cut_node_2),encoding='utf-8')
        mat=mat.reindex(index=my_sort_list,columns=my_sort_list)
        self.Y=mat.to_numpy()



    # 根据节点类型生成对应的雅可比矩阵
    def cal_power_flow(self):

        # 矩阵
        self.G = self.Y.real
        self.B = self.Y.imag

        # 生成雅可比矩阵,默认只有1个平衡节点,其余全是PQ+PV
        self.Jacobbi=np.zeros([(len(self.bus_data_dict_list)-1)*2, (len(self.bus_data_dict_list)-1)*2])
        # 电压修正量
        self.e=np.zeros([len(self.bus_data_dict_list)])
        self.f = np.zeros([len(self.bus_data_dict_list)] )
        # 功率修正量
        self.PQUs = np.zeros([(len(self.bus_data_dict_list)-1)*2])
        self.PQU =  np.zeros([(len(self.bus_data_dict_list)-1) * 2])
        self.dPQU=  np.zeros([(len(self.bus_data_dict_list)-1) * 2])


       # 先读取数据形成节点电压以及电压初始状态
        for node_num,bus_data_dict in enumerate(self.bus_data_dict_list):
            node_type = (bus_data_dict['Type'])


            #  如果是平衡节点
            if node_type==3:
                voltage = bus_data_dict['FinV']
                angle = bus_data_dict['FinAngle']
                voltage_comp = voltage * complex(math.cos(angle), math.sin(angle))
                self.e[node_num] = voltage_comp.real
                self.f[node_num] = voltage_comp.imag

            # 如果是PV节点
            elif node_type ==2:
                voltage = bus_data_dict['FinV']
                self.PQUs[2*node_num] = (bus_data_dict['GenP']-bus_data_dict['LoadP'])
                self.PQUs[2*node_num+1] = voltage*voltage
                self.e[node_num] = voltage
                self.f[node_num] = 0


            else:
                #PQ节点设置初始电压为1即可
                self.PQUs[2*node_num] = (bus_data_dict['GenP']-bus_data_dict['LoadP'])
                self.PQUs[2*node_num+1] = (bus_data_dict['GenQ']-bus_data_dict['LoadQ'])
                self.e[node_num] = 1
                self.f[node_num] = 0



        while self.circle_status==True:
            # 根据当前节点电压计算先计算PQU,再得出dPQU
            for node_num,bus_data_dict in enumerate(self.bus_data_dict_list):
                node_type = bus_data_dict['Type']

                # 如果是PQ节点
                sum1 = 0
                sum2 = 0
                for j in range(len(self.bus_data_dict_list)):
                    sum1 += self.G[node_num][j] * self.e[j] - self.B[node_num][j] * self.f[j]
                    sum2 += self.G[node_num][j] * self.f[j] + self.B[node_num][j] * self.e[j]


                if node_type==3:
                    pass
                elif node_type==2:
                    self.PQU[2 * node_num] = self.e[node_num] * sum1 + self.f[node_num] * sum2
                    self.PQU[2*node_num+1] = self.e[node_num]**2+self.f[node_num]**2
                else:
                    self.PQU[2 * node_num] = self.e[node_num] * sum1 + self.f[node_num] * sum2
                    self.PQU[2*node_num+1] = self.f[node_num] * sum1 - self.e[node_num] * sum2

            # 计算偏差量
            self.dPQU =self.PQUs-self.PQU

            # 如果不满足精度要求
            if max(abs(self.dPQU)) > self.error_tol:
                if self.circle_count>100:
                    print("潮流不收敛！")
                    # 停止迭代
                    self.circle_status=False
                    self.coverge_status=False
                else:
                    self.circle_count += 1
                    # 开始遍历所有节点,生成雅可比矩阵
                    for i in range(len(self.bus_data_dict_list)-1):
                        node_type = self.bus_data_dict_list[i]['Type']
                        for j in range(len(self.bus_data_dict_list)-1):
                            if i!=j:
                                # Pidej
                                self.Jacobbi[2*i][2*j] = -(self.G[i][j] * self.e[i] + self.B[i][j] * self.f[i])
                                # Pidfj
                                self.Jacobbi[2*i][2*j+1] = self.B[i][j] * self.e[i] - self.G[i][j] * self.f[i]

                                if node_type ==3:
                                    pass
                                elif node_type == 2:
                                    # U2dej
                                    self.Jacobbi[2*i+1][2*j]=0
                                    # U2dfj
                                    self.Jacobbi[2*i+1][2*j+1]=0
                                else :
                                    # Qidej=Pidfj
                                    self.Jacobbi[2 * i + 1][2 * j] = self.Jacobbi[2 * i][2 * j + 1]
                                    # Qidfj=-Pidej
                                    self.Jacobbi[2 * i + 1][2 * j + 1] = -self.Jacobbi[2 * i][2 * j]
                            else:
                                sum1 = 0
                                for k in range(len(self.bus_data_dict_list)):
                                    sum1 += (self.G[i][k] * self.e[k] - self.B[i][k] * self.f[k])
                                sum2 = 0
                                for k in range(len(self.bus_data_dict_list) ):
                                    sum2 += (self.G[i][k] * self.f[k] + self.B[i][k] * self.e[k])
                                    # Pidej
                                    self.Jacobbi[2*i][2*j] = -sum1-(self.G[i][j] * self.e[i] + self.B[i][j] * self.f[i])
                                    # Pidfj
                                    self.Jacobbi[2*i][2*j+1] =-sum2+ self.B[i][j] * self.e[i] - self.G[i][j] * self.f[i]

                                if node_type ==3:
                                    pass
                                elif node_type == 2:
                                    # U2dej
                                    self.Jacobbi[2*i+1][2 *j] = -2*self.e[i]
                                    # U2dfj
                                    self.Jacobbi[2*i+1][2*j + 1] = -2*self.f[i]
                                else:
                                    # Qidej
                                    self.Jacobbi[2 * i + 1][2 * j] = sum2 + self.B[i][j] * self.e[i] - self.G[i][j] * \
                                                                     self.f[i]
                                    # Qidfj
                                    self.Jacobbi[2 * i + 1][2 * j + 1] = -sum1 + (
                                                self.G[i][j] * self.e[i] + self.B[i][j] * self.f[i])


                    try:
                        self.dU=np.linalg.solve(self.Jacobbi, -self.dPQU)
                        # 叠加修正量
                        for i in range(len(self.bus_data_dict_list) - 1):
                            self.e[i] += self.dU[2 * i]
                            self.f[i] += self.dU[2 * i + 1]
                    except:
                        print('该方程无解！')
                        self.circle_status = False
                        self.coverge_status = False


            else:
                # 满足条件则终止当前循环
                self.circle_status = False
                self.coverge_status = True
                print('循环结束！循环次数：',self.circle_count)

                # 计算最终各节点电压
                for i in range(len(self.bus_data_dict_list)):
                    temp_amplitude = math.sqrt(self.e[i] ** 2 + self.f[i] ** 2)
                    temp_phase_angle = math.atan(self.f[i] / self.e[i]) * 180 / math.pi

                    self.u.append(complex(self.e[i], self.f[i]))
                    self.u_amp.append(temp_amplitude)
                    self.u_angle.append(temp_phase_angle)

                # 计算最终各节点功率
                for i in range(len(self.bus_data_dict_list)):
                    I=0
                    Ui = complex(self.e[i], self.f[i])
                    # 计算支路的首端功率
                    for j in range(len(self.bus_data_dict_list)):
                        # 计算节点注入的共轭值
                        Uj = complex(self.e[j], self.f[j])
                        I+=(self.Y[i][j])*(Uj)
                    # 计算各节点的功率 S = 电压 X 注入电流的共轭值
                    self.S_node[i]=Ui*I.conjugate()*100

                # 遍历被连接的节点,先计算线路之间的功率，再计算线损从而得到线路电流
                for branch_data in self.branch_data_dict_list:

                    node_num = branch_data['TNumber']
                    to_node_num = branch_data['ZNumber']


                    r=branch_data['R1']
                    x = branch_data['X1']

                    z=complex(r,x)
                    y=1/z
                    # 找到支路在节点导纳矩阵中对应的位置
                    i = self.my_sort_list.index(node_num)
                    j = self.my_sort_list.index(to_node_num)
                    Ui = complex(self.e[i], self.f[i])
                    Uj = complex(self.e[j], self.f[j])

                    # 计算首端功率，这里不能直接用节点导纳矩阵的互导纳算，因为实际流过电线的充电功率被算在自导纳里了
                    # 换算到基准变比侧
                    # 充电电流


                    S_loss = (Ui/k-Uj)*(((Ui/k-Uj)*y).conjugate())

                    i = math.sqrt(abs(S_loss*y)/3)
                    self.S.append({'branch':(node_num,to_node_num),'loss':S_loss*100,'current':i})



    def output_result(self,number):
        if self.coverge_status==True:
            # 保存读取的参数
            # 保存格式化后的母线数据
            df_bus = pd.DataFrame(self.bus_data_dict_list)
            df_bus.to_csv('data_log/bus_data' + '.csv', encoding='utf-8')

            # 保存格式化后的支路数据
            df_branch = pd.DataFrame(self.branch_data_dict_list)
            df_branch=df_branch.set_index(["TNumber"])
            df_branch.to_csv('data_log/branch_data' + '.csv', encoding='utf-8')

            # 保存潮流计算的解，包括电压以及节点功率
            solution_dict = {'V': self.u, 'V_amp': self.u_amp, 'V_angle': self.u_angle, 'S_node': self.S_node}
            solution_df = pd.DataFrame(solution_dict, index=self.my_sort_list)
            # 按实际节点顺序重排
            solution_df = solution_df.reindex(range(1, len(self.bus_data_dict_list) + 1))
            solution_df.to_csv('data_log/{}_solution.csv'.format(number))
            print('切去{}_{}支路,成功'.format(self.cut_node_1,self.cut_node_2))

            # 保存功率损耗以及电流计算
            my_S = pd.DataFrame(self.S)
            my_S.to_csv('{}_loss_current.csv'.format(number), encoding='utf-8')
        else:
            if self.cut_node_1!=None:
                print('切去{}_{}支路后方程无解'.format(self.cut_node_1,self.cut_node_2))




if __name__ == '__main__':
    my_data = MyData('data.txt')
    my_pf = PowerFlow(my_data)