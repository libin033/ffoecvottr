#coding:utf-8
import sys
import os
import numpy as np
import matplotlib.pyplot as plt





if __name__ == '__main__':
    ifp = open('rtcoef_fy4_1_agri.dat', 'r')
    lines1 = ifp.readlines()
    ifp.close()

    ifp = open('fy4a_agri.dat', 'r')
    lines2 = ifp.readlines()
    ifp.close()

    data1 = lines1[1278:2391]
    data2 = lines2[1262:2375]

    data = []
    for item in data1:
        data.append(item.split())
    data1 = np.array(data, dtype=np.float64)

    data = []
    for item in data2:
        data.append(item.split())
    data2 = np.array(data, dtype=np.float64)

    flag = (data2 > 1.0) | (data2 < -1) | (data1 > 1.0) | (data1 < -1) | (data1 == 0) | (data2 == 0)
    print(data2[flag])

    data1[flag] = np.nan
    data2[flag] = np.nan

    plt.plot(data1, data2, 'r.')
    plt.title('water')
    plt.show()

