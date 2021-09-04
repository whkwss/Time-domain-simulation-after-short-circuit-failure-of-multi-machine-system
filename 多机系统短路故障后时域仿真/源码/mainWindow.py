# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'qtui.py'
#
# Created by: PyQt5 UI code generator 5.13.2
#
# WARNING! All changes made in this file will be lost!


import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QMessageBox

from transient_sim import TransientSim
from ui_mainwindow import Ui_MainWindow

import matplotlib
from matplotlib.figure import Figure
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

class MyWindow(QMainWindow, Ui_MainWindow):

    def __init__(self, parent=None):
        super(MyWindow, self).__init__(parent)
        self.setupUi(self)


        self.fig = Figure()

        self.canvas = FigureCanvas(self.fig)

        self.widget_pic.hboxlayout= QtWidgets.QHBoxLayout(self.widget_pic)
        self.widget_pic.hboxlayout.addWidget(self.canvas)

        # 加载变量列表,绝对速度和相对速度
        var_list = ['omega','theta_opp','theta_abs','P']
        self.com_var.clear()
        for var in var_list:
            self.com_var.addItem(var)  # 将串口设备名添加到组合框中

        self.mySignalSet()
        # 设置视频监控页面
    def signal_import_data(self):
        directory = QtWidgets.QFileDialog.getOpenFileName(None, "选取文件夹", "data/","All Files (*);;Text Files (*.txt)")
        directory= directory[0].split('/')
        file_name = directory[-1]
        self.lab_data.setText(file_name)
        self.transient_sim = TransientSim(file_name)
        self.set_fault_node(self.transient_sim.node_list)
        self.set_fault_type()
        self.set_cut_branch(self.transient_sim.branch_list)


    def mySignalSet(self):
        self.btn_import_data.clicked.connect(self.signal_import_data)
        self.btn_start_simu.clicked.connect(self.run_sim)

        # 修改需要显示的图像
        self.com_var.currentTextChanged.connect(self.draw_pic)

    def draw_pic(self):
        self.fig.clf()
        self.axe = self.fig.add_subplot(111)
        para_name = self.com_var.currentText()
        t,result_list,index_list= self.transient_sim.process_result(para_name)

        for index,result in enumerate (result_list):
            self.axe.plot(t, result, label='{}_{}'.format(para_name,index_list[index]))
        self.axe.set_xlabel('t/s')
        self.axe.set_ylabel("δ/deg")
        self.axe.legend()
        # 必须重新调整子图，不然不会显示坐标
        self.fig.tight_layout()
        self.canvas.draw()

    def set_fault_node(self,node_list):
        # 检测所有存在的节点
        self.com_fault_node.clear()
        for node in node_list:
            self.com_fault_node.addItem(node)
    def set_cut_branch(self,branch_list):
        # 检测所有存在的节点
        self.com_cut_branch.clear()
        for node in branch_list:
            self.com_cut_branch.addItem(node)

    def set_fault_type(self):
        # 检测所有存在的节点
        fault_type_list = ['单相接地', '两相短路', '两相接地', '三相短路']
        self.com_fault_type.clear()
        for fault_type in fault_type_list:
            self.com_fault_type.addItem(fault_type)  # 将串口设备名添加到组合框中

    def run_sim(self):
        # try:
        fault_bus = self.com_fault_node.currentText()
        fault_type = self.com_fault_type.currentText()
        cut_branch = self.com_cut_branch.currentText()
        cut_time = self.line_cut_time.text()
        step = self.line_step.text()
        total_time = self.line_total_time.text()
        self.transient_sim.set_simu(fault_bus,fault_type,cut_branch,step,cut_time,total_time)
        self.transient_sim.start_sim()


        self.draw_pic()



if __name__ == '__main__':

    app = QApplication(sys.argv)
    myWin = MyWindow()
    myWin.show()
    sys.exit(app.exec_())

