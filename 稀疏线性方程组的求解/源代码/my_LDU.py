import numpy as np
from  scipy.sparse import coo_matrix,dok_matrix
from itertools import combinations_with_replacement
import pandas as pd

def myLDU(A):
    # 首先获取A的上三角三元组格式的
    A = np.triu(A)

    A = coo_matrix(A)
    rowInd = A.row
    # 获得矩阵的行数
    node_list = list(set(rowInd))
    for p in node_list[:-1]:
        A = coo_matrix(A,dtype=np.float)
        rowInd = A.row
        colInd = A.col
        A = dok_matrix(A,dtype=np.float)
        # 获取该行中所有的非零元列坐标
        col_index_list = []
        for row_index,i in enumerate(rowInd) :
            if i == p:
                # 根据行的索引找到对应的列
                col_index_list.append(colInd[row_index])
        # 进行规格化计算,由于已经取了对角元，对角元需要保留
        # 对节点p发出的所有互边的边权进行修正
        for j in (col_index_list[1:]):
            # 根据列索引找到对应的元素
            A[p,j] =  A[p,j]/ A[p,p]
        # 消去运算
        # 该节点发出的所有边的排列组合
        col_index_list.sort()
        i_j_list = list(combinations_with_replacement(col_index_list[1:], 2))
        for i_j in (i_j_list):
             i = i_j[0]
             j = i_j[1]
             A[i,j] = A[i,j] - A[p,i]*A[p,j]*A[p,p]
    D =np.diag(A.diagonal())
    dim = (A.shape[0])
    U = np.eye(dim,dtype=np.float)
    U += A-D
    L = U.T

    return L,D,U
# 引入中间矢量y和z，Lz =b,Dy = z,Ux = y
def my_forward_substitution(b,L):

    L = coo_matrix(L,dtype=np.float)
    rowInd = L.row
    colInd = L.col
    # 获取需要的数据后把A变为可以修改的矩阵
    L = dok_matrix(L,dtype=float)
    # 获得矩阵的行数
    node_list = list(set(rowInd))
    node_list.sort()
    z = b
    z = z.astype(float)

    for j in  node_list[:-1]:
        if z[j]!=0:
            # 获取要被消去的行的非零元的列坐标
            row_index_list = []
            for col_index,p in enumerate(colInd) :
                if p == j:
                    # 根据列的索引找到对应的行
                    row_index_list.append(rowInd[col_index])
            row_index_list.sort()
            for  i in (row_index_list[1:]):
                # 根据列索引找到对应的元素
                    z[i] = z[i] - L[i,j]*z[j]
    return z

# 规则化运算
def my_norm_substitution(z,D):
    y = z/D.diagonal()
    y = y.astype(float)
    return y
def my_backward_substitution(y,U):
    U = coo_matrix(U,dtype=np.float)
    rowInd = U.row
    # 获取需要的数据后把A变为可以修改的矩阵
    U = dok_matrix(U,dtype=np.float)
    # 获得矩阵的行数
    node_list = list(set(rowInd))
    node_list.sort(reverse=True)
    x = y
    x = x.astype(float)
    for j in node_list[:-1]:
        if x[j] != 0:
            U = coo_matrix(U)
            rowInd = U.row
            colInd = U.col
            U = dok_matrix(U)
            # 扫描该列非零元的行号
            row_index_list = []
            for col_index, p in enumerate(colInd):
                if p == j:
                    # 根据行的索引找到对应的行
                    row_index_list.append(rowInd[col_index])
            row_index_list.sort(reverse=True)

            for i in (row_index_list[1:]):
                # 根据列索引找到对应的元素
                x[i] = x[i] - U[i, j] * x[j]

    return x

def my_slove(A,b):
    L,D,U = myLDU(A)
    z = my_forward_substitution(b, L)
    y = my_norm_substitution(z,D)
    x = my_backward_substitution(y, U)
    return x
if __name__ == '__main__':
    A = np.array([[2,0,0,0,0,0,1,0,0,0,0,1],
             [0,2,0,0,1,0,0,0,1,0,0,0],
             [0,0,2,0,0,1,0,0,1,0,0,0],
             [0,0,0,2,0,0,1,0,0,1,0,0],
             [0,1,0,0,2,0,0,0,0,0,0,1],
             [0,0,1,0,0,2,0,0,0,1,0,0],
             [1,0,0,1,0,0,2,0,0,0,0,0],
             [0,0,0,0,0,0,0,2,1,0,1,0],
             [0,1,1,0,0,0,0,1,2,0,0,0],
             [0,0,0,1,0,1,0,0,0,2,1,0],
             [0,0,0,0,0,0,0,1,0,1,2,1],
             [1,0,0,0,1,0,0,0,0,0,1,2]],dtype=np.float)
    b = np.array([4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5],dtype=float)
    x1 = my_slove(A,b)
    print('x1:',x1)
    print('x2:',np.linalg.solve(A,b))
    print('x3:', np.matmul(np.linalg.inv(A),b))
