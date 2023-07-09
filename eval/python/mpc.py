from ast import Global
import pandas as pd
import numpy as np
import sys
import warnings
import time
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn import metrics
import time
from math import ceil
from math import log2
import multiprocessing
warnings.filterwarnings('ignore')

# time 
time_coordinator = 0

round_shr = 0 # 1
round_rec = 0 # 1
round_mult = 0 # 2
round_div = 0 # 8
round_argmax = 0
round_feature = 0
round_sigmoid = 0 # 20+8*2+20
round_split = 0
round_predict = 0

online_communication = multiprocessing.Value('d', 0)
offline_communication = multiprocessing.Value('d', 0)

# np.set_printoptions(threshold=np.inf)

clientNum = 4 # for wine, we assume each party has 3 features (use j%3 to get feature_id on specific party)
divided_feature_num = 2 # we set mannually
# for credit, divided_feature_num=6, FEATURE_NUM = 23
# for BostonHousing, divided_feature_num=4, FEATURE_NUM = 13
# for audit_risk, divided_feature_num=4, FEATURE_NUM = 13
# for diabetes, divided_feature_num=2, FEATURE_NUM = 8
# for iris, divided_feature_num=1, FEATURE_NUM = 4
# for cal_housing, divided_feature_num=2, FEATURE_NUM = 8
B = 32
# T_f is the length of fraction part
# T_t is the truncation param
T_f,T_t,M = np.uint64(2**20),np.uint64(2**44),np.uint64(2**64-1)


def gen_uint64():
    ran = np.random.randint(0, 2**64 , size=1, dtype=np.uint64)
    return ran[0]

def double_uint64(db_val): # mult by T_f
    unsigned_long_x = np.dtype(np.uint64) 
    long_long_x = np.int64(db_val*T_f)
    unsigned_long_x = np.uint64(long_long_x)
    return unsigned_long_x

def double_uint64_array(db_array):
    uint64_data = np.zeros_like(db_array, dtype = np.uint64)
    for i in range(0,db_array.shape[0]):
        for j in range(0,db_array.shape[1]):
            uint64_data[i][j] = double_uint64(db_array[i][j])
    return uint64_data

def uint64_double(uint64_val): # div by T_f
    float_res = np.float64(np.int64(uint64_val))
    restore = float_res / T_f
    return restore

def uint64_double_array(uint64_array):
    restore = np.zeros_like(uint64_array, dtype = np.float64)
    for i in range(0,uint64_array.shape[0]):
        for j in range(0,uint64_array.shape[1]):
            restore[i][j] = uint64_double(uint64_array[i][j])
    return restore

def truncate_array(uint64_array):# uint64_shared_array has been expanded by *T_f
    trunc = np.zeros_like(uint64_array, dtype = np.uint64)
    for i in range(0,uint64_array.shape[0]):
        for j in range(0,uint64_array.shape[1]):
            trunc[i][j] = truncate(uint64_array[i][j])
    return trunc

def truncate(num):
    """
    input-
        uint64 number
    output-
        uint64 truncated number
    """
    res = num // T_f # in Python,We use '//' instead of '/'  to conduct integer division when we want to get int() result.
    
    flag = np.int64(num)
    if flag < 0:
        #print("num < 0")
        res= res - T_t
    return res



shared_A_Train_X = [0,0,0,0] # ClientA's shared data in client 1,2,3,4
shared_B_Train_X = [0,0,0,0] # ClientB's shared data in client 1,2,3,4
shared_C_Train_X = [0,0,0,0] # ClientC's shared data in client 1,2,3,4
shared_D_Train_X = [0,0,0,0] # ClientD's shared data in client 1,2,3,4



def getData():
    dataset = pd.read_csv("../data/cal_housing.csv")
    X, y = dataset.iloc[:, :-1], dataset.iloc[:, -1:]
    train_X, test_X, train_y, test_y = train_test_split(X, y,
					test_size = 0.2, random_state = 9595)
    feature_num = train_X.shape[1]
    print(feature_num)
    # Client A
    A_Train_X = train_X.iloc[:,0:divided_feature_num*1].values
    #Client B
    B_Train_X = train_X.iloc[:,divided_feature_num*1:divided_feature_num*2].values
    #Client C
    C_Train_X = train_X.iloc[:,divided_feature_num*2:divided_feature_num*3].values    
    #Client D
    D_Train_X = train_X.iloc[:,divided_feature_num*3:].values
    D_Train_y = train_y.values
    return A_Train_X,B_Train_X,C_Train_X,D_Train_X,D_Train_y,test_X.values,test_y.values # no return value, train data are gloabal variables.

# 1 round $$$
def SHR_array(u64_array):
    """
    Share raw data in secret sharing
     input
        - uint64
     output
        - uint64
    """
    tmp_array = np.copy(u64_array)
    ran = np.array([np.random.randint(0, 2**64 , size=(tmp_array.shape[0],tmp_array.shape[1]), dtype=np.uint64) for i in range(clientNum - 1)])
    tmp_array -= np.sum(ran, axis=0).astype('uint64') # x0 = x - x1 - x2, self_value's shape is the same as data
    tmp_array = np.expand_dims(tmp_array, axis=0) # Turn u64_data into array_element
    dataList = np.concatenate([ran, tmp_array], axis=0) # 4 elements in datalist
    return dataList

def SHR_array_test(u64_array):
    """
    Share raw data in secret sharing
     input
        - uint64
     output
        - uint64
    """
    tmp_array = np.copy(u64_array)
    ran = np.array([np.random.randint(0, 2**64 , size=(tmp_array.shape[0],tmp_array.shape[1]), dtype=np.uint64) for i in range(clientNum - 1)])
    tmp_array -= np.sum(ran, axis=0).astype('uint64') # x0 = x - x1 - x2, self_value's shape is the same as data
    tmp_array = np.expand_dims(tmp_array, axis=0) # Turn u64_data into array_element
    dataList = np.concatenate([ran, tmp_array], axis=0) # 4 elements in datalist
    return dataList

def SHR_array_online(u64_array):
    global online_communication
    online_communication.value += 64*(clientNum-1)*u64_array.shape[0]*u64_array.shape[1]
    """
    Share raw data in secret sharing
     input
        - uint64
     output
        - uint64
    """
    tmp_array = np.copy(u64_array)
    ran = np.array([np.random.randint(0, 2**64 , size=(tmp_array.shape[0],tmp_array.shape[1]), dtype=np.uint64) for i in range(clientNum - 1)])
    tmp_array -= np.sum(ran, axis=0).astype('uint64') # x0 = x - x1 - x2, self_value's shape is the same as data
    tmp_array = np.expand_dims(tmp_array, axis=0) # Turn u64_data into array_element
    dataList = np.concatenate([ran, tmp_array], axis=0) # 4 elements in datalist
    return dataList

def SHR_array_online(u64_array):
    global online_communication
    online_communication.value += 64*(clientNum-1)*u64_array.shape[0]*u64_array.shape[1]
    """
    Share raw data in secret sharing
     input
        - uint64
     output
        - uint64
    """
    tmp_array = np.copy(u64_array)
    ran = np.array([np.random.randint(0, 2**64 , size=(tmp_array.shape[0],tmp_array.shape[1]), dtype=np.uint64) for i in range(clientNum - 1)])
    tmp_array -= np.sum(ran, axis=0).astype('uint64') # x0 = x - x1 - x2, self_value's shape is the same as data
    tmp_array = np.expand_dims(tmp_array, axis=0) # Turn u64_data into array_element
    dataList = np.concatenate([ran, tmp_array], axis=0) # 4 elements in datalist
    return dataList

def SHR_array_offline(u64_array):
    """
    Share raw data in secret sharing
     input
        - uint64
     output
        - uint64
    """
    global offline_communication
    offline_communication.value += 64*(clientNum)*u64_array.shape[0]*u64_array.shape[1]

    tmp_array = np.copy(u64_array)
    ran = np.array([np.random.randint(0, 2**64 , size=(tmp_array.shape[0],tmp_array.shape[1]), dtype=np.uint64) for i in range(clientNum - 1)])
    tmp_array -= np.sum(ran, axis=0).astype('uint64') # x0 = x - x1 - x2, self_value's shape is the same as data
    tmp_array = np.expand_dims(tmp_array, axis=0) # Turn u64_data into array_element
    dataList = np.concatenate([ran, tmp_array], axis=0) # 4 elements in datalist
    return dataList
# 1 round $$$
def Rec(shared_list):
    restore_uint64 = np.sum(shared_list).astype('uint64')
    return uint64_double(restore_uint64)

# 1 round $$$
def Rec_array(shared_list):
    restore_uint64 = np.sum(shared_list,axis=0).astype('uint64')
    return uint64_double_array(restore_uint64)

def Rec_array_test(shared_list):
    restore_uint64 = np.sum(shared_list,axis=0).astype('uint64')
    return uint64_double_array(restore_uint64)

def Rec_array_online(shared_list):
    global online_communication
    online_communication.value += 64*(clientNum-1)*shared_list.shape[0]*shared_list.shape[1]
    restore_uint64 = np.sum(shared_list,axis=0).astype('uint64')
    return uint64_double_array(restore_uint64)

# 1 round $$$
# Offline 1 round
# Online 2 ROUND
# Round1: reconstruct(e,f) # parallel #不需要算进去
# Round2: customed Truncation $$$
def SSMult_array(shared_X,shared_Y): # Vectorized
    """
        Input -- 2 lists of 4 matrix elements: shared_X = n*m shared_Y=m*k (shape = 4,5000,3)
        Output-- 1 list of 4 matrix elements
    """

    A = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],shared_X[0].shape[1]), dtype=np.uint64)
    A2 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],shared_X[0].shape[1]), dtype=np.uint64)
    A3 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],shared_X[0].shape[1]), dtype=np.uint64)
    A4 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],shared_X[0].shape[1]), dtype=np.uint64)
    A1 = A - A2 - A3 - A4

    B = np.random.randint(0, 2**64 , size=(shared_Y[0].shape[0],shared_Y[0].shape[1]), dtype=np.uint64)
    B2 = np.random.randint(0, 2**64 , size=(shared_Y[0].shape[0],shared_Y[0].shape[1]), dtype=np.uint64)
    B3 = np.random.randint(0, 2**64 , size=(shared_Y[0].shape[0],shared_Y[0].shape[1]), dtype=np.uint64)
    B4 = np.random.randint(0, 2**64 , size=(shared_Y[0].shape[0],shared_Y[0].shape[1]), dtype=np.uint64)
    B1 = B - B2 - B3 - B4

    C = A.dot(B)
    C2 = np.random.randint(0, 2**64 , size=(C.shape[0],C.shape[1]), dtype=np.uint64)
    C3 = np.random.randint(0, 2**64 , size=(C.shape[0],C.shape[1]), dtype=np.uint64)
    C4 = np.random.randint(0, 2**64 , size=(C.shape[0],C.shape[1]), dtype=np.uint64)
    C1 = C - C2 - C3 - C4

    # Clients' part
    # shared_E = shared_X - A # list[clientNum]
    # shared_F = shared_Y - B # list[clientNum]
    shared_E1 = shared_X[0] - A1
    shared_E2 = shared_X[1] - A2
    shared_E3 = shared_X[2] - A3
    shared_E4 = shared_X[3] - A4

    shared_F1 = shared_Y[0] - B1
    shared_F2 = shared_Y[1] - B2
    shared_F3 = shared_Y[2] - B3
    shared_F4 = shared_Y[3] - B4



    # aggregate E,F
    E = shared_E1 + shared_E2 + shared_E3 + shared_E4
    F = shared_F1 + shared_F2 + shared_F3 + shared_F4
    #print(E,F)

    res0 = E.dot(F) + E.dot(B1) + A1.dot(F) + C1
    res1 = E.dot(B2) + A2.dot(F) + C2
    res2 = E.dot(B3) + A3.dot(F) + C3
    res3 = E.dot(B4) + A4.dot(F) + C4

    # truncate
    res0 = truncate_array(res0)
    # truncate: client2 gather shared from 3,4
    ran2 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],shared_Y[0].shape[1]), dtype=np.uint64)
    ran3 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],shared_Y[0].shape[1]), dtype=np.uint64)
    res1 = truncate_array(res1 + res2 + res3) - ran2 - ran3
    res2 = ran2
    res3 = ran3
    res = np.array([res0,res1,res2,res3])
    return res

# 总计 1 round $$$
# Offline 1 round
# Online 2 ROUND
# Round1: reconstruct(e,f) # parallel #不需要算进去
# Round2: customed Truncation $$$
# 减去triple生成的时间
def SSHadamardProduct(shared_X,shared_indicator):

    """
        Input: shared_X is (4,5000,1) (4,5000,1)
        Output: shared = (4,5000,1)
    """
    global online_communication, offline_communication
    

    time_start=time.time()

    A = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    A2 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    A3 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    A4 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    A1 = A - A2 - A3 - A4
    offline_communication.value += 64*(clientNum)*shared_X.shape[0]*shared_X.shape[1]

    B = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    B2 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    B3 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    B4 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    B1 = B - B2 - B3 - B4
    offline_communication.value += 64*(clientNum)*shared_X.shape[0]*shared_X.shape[1]

    C = A*B
    C2 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    C3 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    C4 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    C1 = C - C2 - C3 - C4
    offline_communication.value += 64*(clientNum)*shared_X.shape[0]*shared_X.shape[1]

    time_end=time.time()
    time_mult_beaver = time_end - time_start
    global time_coordinator 
    time_coordinator += time_mult_beaver

    # Clients' part
    # shared_E = shared_X - A # list[clientNum]
    # shared_F = shared_Y - B # list[clientNum]
    shared_E1 = shared_X[0] - A1
    shared_E2 = shared_X[1] - A2
    shared_E3 = shared_X[2] - A3
    shared_E4 = shared_X[3] - A4

    shared_F1 = shared_indicator[0] - B1
    shared_F2 = shared_indicator[1] - B2
    shared_F3 = shared_indicator[2] - B3
    shared_F4 = shared_indicator[3] - B4



    # Rec E,F
    E = shared_E1 + shared_E2 + shared_E3 + shared_E4
    online_communication.value += 64*(clientNum-1)*shared_X.shape[0]*shared_X.shape[1]
    F = shared_F1 + shared_F2 + shared_F3 + shared_F4
    online_communication.value += 64*(clientNum-1)*shared_X.shape[0]*shared_X.shape[1]

    res0 = E*F + E*B1 + A1*F + C1
    res1 = E*B2 + A2*F + C2
    res2 = E*B3 + A3*F + C3
    res3 = E*B4 + A4*F + C4

    # truncate
    res0 = truncate_array(res0)
    # truncate: client2 gather shared from 3,4
    ran2 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    ran3 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    res1 = truncate_array(res1 + res2 + res3) - ran2 - ran3
    online_communication.value += 64*(2)*shared_X.shape[0]*shared_X.shape[1]
    res2 = ran2
    res3 = ran3
    res = np.array([res0,res1,res2,res3])
    return res

def SSHadamardProduct_test(shared_X,shared_indicator):

    """
        Input: shared_X is (4,5000,1) (4,5000,1)
        Output: shared = (4,5000,1)
    """
    
    time_start=time.time()

    A = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    A2 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    A3 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    A4 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    A1 = A - A2 - A3 - A4

    B = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    B2 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    B3 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    B4 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    B1 = B - B2 - B3 - B4

    C = A*B
    C2 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    C3 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    C4 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    C1 = C - C2 - C3 - C4

    time_end=time.time()
    time_mult_beaver = time_end - time_start
    global time_coordinator 
    time_coordinator += time_mult_beaver

    # Clients' part
    # shared_E = shared_X - A # list[clientNum]
    # shared_F = shared_Y - B # list[clientNum]
    shared_E1 = shared_X[0] - A1
    shared_E2 = shared_X[1] - A2
    shared_E3 = shared_X[2] - A3
    shared_E4 = shared_X[3] - A4

    shared_F1 = shared_indicator[0] - B1
    shared_F2 = shared_indicator[1] - B2
    shared_F3 = shared_indicator[2] - B3
    shared_F4 = shared_indicator[3] - B4



    # Rec E,F
    E = shared_E1 + shared_E2 + shared_E3 + shared_E4
    F = shared_F1 + shared_F2 + shared_F3 + shared_F4

    res0 = E*F + E*B1 + A1*F + C1
    res1 = E*B2 + A2*F + C2
    res2 = E*B3 + A3*F + C3
    res3 = E*B4 + A4*F + C4

    # truncate
    res0 = truncate_array(res0)
    # truncate: client2 gather shared from 3,4
    ran2 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    ran3 = np.random.randint(0, 2**64 , size=(shared_X[0].shape[0],1), dtype=np.uint64)
    res1 = truncate_array(res1 + res2 + res3) - ran2 - ran3
    res2 = ran2
    res3 = ran3
    res = np.array([res0,res1,res2,res3])
    return res

# 20 次乘法
def SSDiv(shared_x,shared_y):
    """
    Newton-Raphson algorithm
    Y = 2**6, 
    z0 = 1/Y = 1/2**6
    ..........................................
    Input: uint64 shared_x,shared_y = (4,1,1)
    Output: uint64 res = (4,)
    ..........................................
    """
    shared_z0 = SHR_array_test(double_uint64_array(np.array([1.0/2**20],ndmin = 2))) # 3*logY > 20 会爆掉， 我们用100放缩分母
    for i in range(20):
        shared_z0_square = SSHadamardProduct(shared_z0,shared_z0) # return (4,1,1)
        shared_y_mult_z0_square = SSHadamardProduct(shared_y,shared_z0_square)
        shared_z1 = 2*shared_z0 - shared_y_mult_z0_square # 可以直接用2乘吗？要转uint64再share吗？
        shared_z0 = shared_z1
    return SSHadamardProduct(shared_z0,shared_x)

def SSDiv_test(shared_x,shared_y):
    """
    Newton-Raphson algorithm
    Y = 2**6, 
    z0 = 1/Y = 1/2**6
    ..........................................
    Input: uint64 shared_x,shared_y = (4,1,1)
    Output: uint64 res = (4,)
    ..........................................
    """
    shared_z0 = SHR_array_test(double_uint64_array(np.array([1.0/2**20],ndmin = 2))) # 3*logY > 20 会爆掉， 我们用100放缩分母
    for i in range(20):
        shared_z0_square = SSHadamardProduct_test(shared_z0,shared_z0) # return (4,1,1)
        shared_y_mult_z0_square = SSHadamardProduct_test(shared_y,shared_z0_square)
        shared_z1 = 2*shared_z0 - shared_y_mult_z0_square # 可以直接用2乘吗？要转uint64再share吗？
        shared_z0 = shared_z1
    return SSHadamardProduct_test(shared_z0,shared_x)

def SSDiv100(shared_x,shared_y):
    shared_100 = SHR_array_test(double_uint64_array(np.array([100],ndmin = 2)))
    shared_y_div_100 = SSDiv_test(shared_y,shared_100)
    tmp = SSDiv_test(shared_x,shared_y_div_100)
    res = SSDiv(tmp,shared_100)
    return res



def getParty(feature_j):
    # for wine.csv
    return (feature_j // divided_feature_num)+1

def shift(xs, n):
    if n >= 0:
        return np.r_[np.full(n, xs[0]), xs[:-n]]
    else:
        return np.r_[xs[-n:], np.full(-n, False)]

def shift_array(G, n):
    C0 = shift(G[0],n)
    C1 = shift(G[1],n)
    C2 = shift(G[2],n)
    C3 = shift(G[3],n)
    C = np.array([C0,C1,C2,C3])
    return C

def num2BoolCom(alpha):
    a = bin(alpha& 0xffffffffffffffff).replace('0b', '')[::-1]
    alpha = np.zeros(64, dtype=bool)
    for i in range(len(a)):
        if a[i]=='1':
            alpha[63-i]=True
    #print(alpha)
    return alpha

def bool2NumCom(b_seq):
    res = 0
    for i in range(64):
        if b_res[0][63-i] == True:
            res += 2**i
    return bin(res)

def b_rec(X):
    b_res = X[0] ^ X[1] ^ X[2] ^ X[3]
    res = 0
    for i in range(64):
        if b_res[63-i] == True:
            res += 2**i
    return res

def bit_decomposition(X):
    """
    Input: (4,N,1)
    Output:(4,N,64)
    """
    X_boolean = np.ones((clientNum, 64), dtype=bool)
    for i in range(clientNum):
            X_boolean[i] = num2BoolCom(X[i])
    return X_boolean

def bit_share(shared_X):
    shared_X0 = np.array([shared_X[0],0,0,0],dtype=np.uint64)
    shared_X1 = np.array([0,shared_X[1],0,0],dtype=np.uint64)
    shared_X2 = np.array([0,0,shared_X[2],0],dtype=np.uint64)
    shared_X3 = np.array([0,0,0,shared_X[3]],dtype=np.uint64)
    b_shared_X0 = bit_decomposition(shared_X0) #(4,1,64)
    b_shared_X1 = bit_decomposition(shared_X1) #(4,1,64)
    b_shared_X2 = bit_decomposition(shared_X2) #(4,1,64)
    b_shared_X3 = bit_decomposition(shared_X3) #(4,1,64)
    return b_shared_X0,b_shared_X1,b_shared_X2,b_shared_X3


def ppa4(x, y):
    """
    Parallel prefix adder (PPA), using the Kogge-Stone adder topology.
    """

    k = 64
    keep_masks = []
    for i in range(ceil(log2(k))):
        keep_masks.append((1 << (2 ** i)) - 1)
    """
    For example, if prot.nbits = 64, then keep_masks is:
    keep_masks = [0x0000000000000001, 0x0000000000000003, 0x000000000000000f,
                0x00000000000000ff, 0x000000000000ffff, 0x00000000ffffffff]
    """

    G = B_AND(x,y)
    P = x ^ y
    for i in range(ceil(log2(k))):
        # k_mask = prot.define_constant(
        #     np.ones(x.shape, dtype=np.object) * keep_masks[i],
        #     apply_scaling=False,
        #     share_type=BOOLEAN,
        # )

        # G1 = G << (2 ** i)
        # P1 = P << (2 ** i)
        G1 = shift_array(G,-(2 ** i))
        P1 = shift_array(P,-(2 ** i))
        
        b_mask = num2BoolCom(keep_masks[i])
        P1 = P1 ^ b_mask


        # G = G ^ (P & G1)
        # P = P & P1
        G = G ^ B_AND(P,G1)
        P = B_AND(P,P1)
    # G stores the carry-in to the next position
    #C = G << 1
    C = shift_array(G,-1)
    P = x ^ y
    res = C ^ P
    return res

def B_AND(x,y):
    """
    Input  = 4,64
    Output = 4,64
    """
    # x=(x0,x1), y=(y0,y1)
    z = np.ones((4,64), dtype=bool)

    for i in range(64):
        z[0][i],z[1][i],z[2][i],z[3][i] = mpcMultiBool(x[0][i],x[1][i],x[2][i], x[3][i], y[0][i], y[1][i], y[2][i], y[3][i])
    return z

def mpcMultiBool(x1, x2, x3, x4, y1, y2, y3, y4) :
    global online_communication, offline_communication

    # 生乘triple
    u = False  # np.random.randint(0, math.sqrt(MAX))
    v = True  ##np.random.randint(0, math.sqrt(MAX))
    w = u & v
    a1 = True  # np.random.randint(0, math.sqrt(MAX))
    a2 = True
    a3 = True
    a4 = u ^ a1 ^ a2 ^ a3
    offline_communication.value += clientNum

    b1 = True  # np.random.randint(0, math.sqrt(MAX))
    b2 = True
    b3 = True
    b4 = v ^ b1 ^ b2 ^ b3
    offline_communication.value += clientNum
    
    c1 = True  # np.random.randint(0, math.sqrt(MAX))
    c2 = True
    c3 = True
    c4 = w ^ a1 ^ a2 ^ a3
    offline_communication.value += clientNum

    e = x1 ^ x2 ^ x3 ^ x4 ^ a1 ^ a2 ^ a3 ^ a4 
    online_communication.value += clientNum-1
    f = y1 ^ y2 ^ y3 ^ y4 ^ b1 ^ b2 ^ b3 ^ b4 
    online_communication.value += clientNum-1

    res1 = e & f ^ f & a1 ^ e & b1 ^ c1
    res2 = f & a2 ^ e & b2 ^ c2
    res3 = f & a3 ^ e & b3 ^ c3
    res4 = f & a4 ^ e & b4 ^ c4

    # return res1^ res2^ res3^ res4
    return res1, res2, res3, res4

def msb_b(shared_test):
    b_shared_X0,b_shared_X1,b_shared_X2,b_shared_X3 = bit_share(shared_test)
    tmp = ppa4(b_shared_X0,b_shared_X1)
    # print(b_rec(tmp))

    tmp1 = ppa4(tmp,b_shared_X2)
    tmp2 = ppa4(tmp1,b_shared_X3)
    return np.array([[tmp2[0][0]],[tmp2[1][0]],[tmp2[2][0]],[tmp2[3][0]]])


def B2A(shared_bit):
    global online_communication, offline_communication

    r = np.random.choice([True,False])
    r0 = False
    r1 = True
    r2 = False
    r3 = r^r0^r1^r2
    offline_communication.value += clientNum

    arith_r = np.array([[r]], dtype = np.float64)
    shared_r = SHR_array_offline(double_uint64_array(arith_r))

    arith_1 = np.array([[1]], dtype = np.float64)
    shared_1 = SHR_array_test(double_uint64_array(arith_1))


    z0 = shared_bit[0][0] ^ r0
    z1 = shared_bit[1][0] ^ r1
    z2 = shared_bit[2][0] ^ r2
    z3 = shared_bit[3][0] ^ r3

    z = z0^z1^z2^z3
    offline_communication.value += clientNum

    if z == True:
        b_arith = shared_1 - shared_r
    else:
        b_arith = shared_r

    return b_arith

def SSArgMax(shared_array):
    size_array = shared_array[0].shape[0]

    # set the first element as the biggest one
    tmp = np.array([[0]], dtype = np.float64)
    shared_max = SHR_array_test(double_uint64_array(tmp)) #extra round
    

    shared_max[0] = shared_array[0][0]
    shared_max[1] = shared_array[1][0]
    shared_max[2] = shared_array[2][0]
    shared_max[3] = shared_array[3][0]

    shared_max_index = SHR_array_test(double_uint64_array(tmp)) 

    arith_1 = np.array([[1]], dtype = np.float64)
    shared_1 = SHR_array_test(double_uint64_array(arith_1)) #extra round

    for i in range(1,size_array):

        #print("max value=",Rec_array(shared_max))

        shared_tmp = SHR_array(double_uint64_array(tmp)) #extra round
        shared_tmp[0] = shared_array[0][i]
        shared_tmp[1] = shared_array[1][i]
        shared_tmp[2] = shared_array[2][i]
        shared_tmp[3] = shared_array[3][i]

        #print("cur value=",Rec_array(shared_tmp))

        shared_xMINUSy = shared_max - shared_tmp
        sign = B2A(msb_b(shared_xMINUSy)) # +1 
        shared_max = SSHadamardProduct(shared_xMINUSy,shared_1-sign) + shared_tmp # +1

        cur_index = np.array([[i]], dtype = np.float64)
        shared_cur_index = SHR_array_test(double_uint64_array(cur_index)) 
        shared_max_index = SSHadamardProduct(shared_max_index-shared_cur_index,shared_1-sign) + shared_cur_index # +1

        #print("shared_max_index=",Rec_array(shared_max_index))
    max_index = Rec_array_online(shared_max_index) # +1

    return ceil(max_index[0][0])







def SSSigmoid_fm(shared_y_predict):
    global round_sigmoid
    round_sigmoid += (20+8*2+20)
    shared_res = np.zeros((4,shared_y_predict[0].shape[0],1), dtype = np.uint64)

    m = 2**8
    shared_m = SHR_array_test(double_uint64_array(m*np.ones(shape=(shared_y_predict[0].shape[0],1),dtype = np.float64)))
    shared_1 = SHR_array_test(double_uint64_array(np.ones(shape=(shared_y_predict[0].shape[0],1),dtype = np.float64)))


    shared_x_div_m = SSDiv(shared_y_predict,shared_m)
    shared_term = shared_1 + shared_x_div_m
    shared_term_2 = SSHadamardProduct(shared_term,shared_term)
    shared_term_4 = SSHadamardProduct(shared_term_2,shared_term_2)
    shared_term_8 = SSHadamardProduct(shared_term_4,shared_term_4)
    shared_term_16 = SSHadamardProduct(shared_term_8,shared_term_8)
    shared_term_32 = SSHadamardProduct(shared_term_16,shared_term_16)
    shared_term_64 = SSHadamardProduct(shared_term_32,shared_term_32)
    shared_term_128 = SSHadamardProduct(shared_term_64,shared_term_64)
    shared_term_256 = SSHadamardProduct(shared_term_128,shared_term_128)


    shared_res = SSDiv(shared_term_256,shared_term_256 + shared_1)

    return shared_res




if __name__ == "__main__":
    test_1 = 1*np.ones((10,1),dtype=np.float64)
    test_4 = 4*np.ones((10,1),dtype=np.float64)
    shared_test_1 = SHR_array_online(double_uint64_array(test_1))
    shared_test_4 = SHR_array_online(double_uint64_array(test_4))
    print(Rec_array(SSDiv(shared_test_1,shared_test_4)))