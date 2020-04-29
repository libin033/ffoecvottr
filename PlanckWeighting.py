#coding:utf-8
import os
import sys
import numpy as np
import datetime
import struct

exepath = os.path.dirname(__file__)
sys.path.append(exepath)
from config import *


def Planck_t2r(temp, freq):
    '''
    Units: mW/(m2.cm-1.sr
    :param temp: bright temperature
    :param freq: wave number
    :return:
    '''
    # Radiation_C1 = 1.191042953E-5
    # Radiation_C2 = 1.4387774
    a = Radiation_C1 * freq * freq * freq
    b = (Radiation_C2 * freq) / temp
    c = np.exp(b) - 1.0
    radi = a / c
    return radi

def Planck_r2t(rad, freq):
    '''
    Units: mW/(m2.cm-1.sr
    :param temp: bright temperature
    :param freq: wave number
    :return:
    '''

    # Radiation_C1 = 1.191042953E-5
    # Radiation_C2 = 1.4387774
    vs = freq
    bt = (Radiation_C2*vs/np.log(Radiation_C1*vs*vs*vs/(rad)+1.0))

    return bt

def GetProfData(strfilename):
#
    if not os.path.isfile(strfilename):
        raise Exception(strfilename + ' is not exist, will exit...')

    data_prof = np.full(shape=(54, 8), fill_value=-9999.0, dtype=np.float32)

    fp = open(strfilename, 'rb')
    for i in range(54):
        for j in range(8):
            data = fp.read(4)
            data_prof[i, j] = struct.unpack("f", data)[0]
    fp.close()

    return data_prof

def GetSurfData(strfilename):
# lat,lon,psfc,landseamask,t2m,u10m,v10m,q2m,ts
    if not os.path.isfile(strfilename):
        raise Exception(strfilename + ' is not exist, will exit...')

    data_surf = np.full(shape=(9, ), fill_value=-9999.0, dtype=np.float32)

    fp = open(strfilename, 'rb')
    for i in range(9):
        data = fp.read(4)
        data_surf[i] = struct.unpack("f", data)[0]
    fp.close()

    return data_surf

def InterpLevel2Layer(data_in, data_press=None):

    data_out = np.full(shape=(53, ), fill_value= FILL_VALUE, dtype= np.float32)

    for i in range(53) :
        x1 = np.log(Dict_Chan_Info['Pressure_REF'][i])
        x2 = np.log(Dict_Chan_Info['Pressure_REF'][i+1])

        y1 = data_in[i]
        y2 = data_in[i+1]

        x = (x1 + x2) * 0.5

        data_out[i] = ((x2 - x) * y1 + (x - x1) * y2) / (x2 - x1)

    return data_out

def PlanckWeighting():
    '''
    普朗克权重计算
    :return: 53层的温度
    '''
    # P, T, Wv, CO2, O3, N2O, CO, CH4
    data_temp = GetProfData(strProfFileName)
    # data_rad = np.full(shape=(53, ), fill_value=FILL_VALUE, dtype=np.float32)
    # for i in range(8):
    #     data_prof[:, i] = InterpLevel2Layer(data_temp[:, i])


    return InterpLevel2Layer(data_temp[:, 1])

if __name__ == '__main__':

    PlanckWeighting()

    # data_temp = GetProfData(strProfFileName)
    # # data_surf = GetSurfData(strSurfFileName)
    #
    # data_prof = np.full(shape=(53, 8), fill_value=FILL_VALUE, dtype=np.float32)
    # for i in range(8) :
    #     data_prof[:, i] = InterpLevel2Layer(data_temp[:, i])
    #     # print(data_prof[:, i])








