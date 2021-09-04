class MyData:
    def __init__(self,data_path):
        self.data_block_list = ['BUS','BRANCH','GEN_T']
        self.bus_data_name_list=['BusNumber','Type','FinV','FinAngle','LoadP','LoadQ','GenP','GenQ']
        self.branch_data_name_list =['BranchName','TNumber','ZNumber','R1','X1','R2','X2','R0','X0']
        self.gen_data_name_list=['Name','BUS','Xd1','X2','Tj','Xt']

        # 初始化母线和支路数据
        self.bus_data_dict_list=[]
        self.branch_data_dict_list = []
        self.gen_data_dict_list=[]
        self.read_data(data_path)


    def read_data(self,data_path):
        with open('data/'+data_path, 'r', encoding='utf-8') as f:
            # data_process = 'data_process.txt'
            data_lines =f.readlines()
            nums =len(data_lines)
            rows_get=[]
            for i in range(nums):
                line =data_lines[i]
                if '//' in line:
                    pass
                else:
                    rows_get.append(data_lines[i])

            # 删除注释的文本
            data_processed = ''
            for data_line in rows_get:
                data_processed += data_line
            with open('data_log/data_processed.txt','w',encoding='utf-8') as f:
                f.write(data_processed)
            data_block=data_processed.split('\n\n')

            data_block_dict = dict(zip(self.data_block_list,data_block))
            #
            bus_data = data_block_dict['BUS'].split('\n')[2:]
            branch_data = data_block_dict['BRANCH'].split('\n')[2:]
            gen_data =data_block_dict['GEN_T'].split('\n')[2:]

            # 读取总线数据
            for data_line in bus_data:
                data_list=[]
                # 数据列表
                data_list.append(int(data_line[1-1:12-1-1]))
                data_list.append(int(data_line[12-1:17-1-1]))
                data_list.append(float(data_line[17-1:23-1-1]))
                data_list.append(float(data_line[23-1:32-1-1]))
                data_list.append(float(data_line[32-1:39-1-1]))
                data_list.append(float(data_line[39-1:46-1-1]))
                data_list.append(float(data_line[46-1:51-1-1]))
                data_list.append(float(data_line[51-1:60-1-1]))

                data_dict = dict(zip(self.bus_data_name_list,data_list))
                self.bus_data_dict_list.append(data_dict)

            for data_line in branch_data:
                data_list = []
                # 数据列表
                data_list.append(int(data_line[1 - 1:9 - 1 - 1]))
                data_list.append(int(data_line[9 - 1:17 - 1 - 1]))
                data_list.append(int(data_line[17 - 1:22 - 1 - 1]))
                data_list.append(float(data_line[22 - 1:26 - 1 - 1]))
                data_list.append(float(data_line[26 - 1:32 - 1 - 1]))
                data_list.append(float(data_line[32 - 1:36 - 1 - 1]))
                data_list.append(float(data_line[36 - 1:44 - 1 - 1]))
                data_list.append(float(data_line[44 - 1:49 - 1 - 1]))
                data_list.append(float(data_line[49 - 1:55 - 1 - 1]))

                data_dict = dict(zip(self.branch_data_name_list,data_list))
                self.branch_data_dict_list.append(data_dict)

            for data_line in gen_data:
                data_list = []
                # 数据列表
                data_list.append(int(data_line[1 - 1:9 - 1 - 1]))
                data_list.append(int(data_line[9 - 1:14 - 1 - 1]))
                data_list.append(float(data_line[14 - 1:20 - 1 - 1]))
                data_list.append(float(data_line[20 - 1:26 - 1 - 1]))
                data_list.append(float(data_line[26 - 1:31 - 1 - 1]))
                data_list.append(float(data_line[31 - 1:40 - 1 - 1]))

                data_dict = dict(zip(self.gen_data_name_list, data_list))
                self.gen_data_dict_list.append(data_dict)



if __name__ == '__main__':
    data = MyData()