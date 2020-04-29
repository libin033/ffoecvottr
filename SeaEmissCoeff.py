#coding:utf-8
'''
计算海洋发射率系数

'''



import sys
import os
import numpy as np
from scipy import optimize
from matplotlib import pyplot as plt
from scipy import interpolate


exepath = os.path.dirname(__file__)
sys.path.append(exepath)
sys.path.append(os.path.join(exepath, '..'))

from config import *


def func(p, x):
    a = p[0]
    b = p[1]
    c = p[2]

    return a + b * x + c * x**2


def least_square_w(p, x, y):
    '''
    y = a - b * x - c * x**2
    最小二乘法残差计算（自变量为风速 w）
    :param p:
    :param x:
    :param y:
    :return:
    '''
    a = p[0]
    b = p[1]
    c = p[2]

    return (a + b * x + c * x**2) - y

def least_squareC(p, x1, x2, y):
    '''
    最小二乘法残差计算c7、c8、c9
    :param p:
    :param x1: 卫星天顶角
    :param x2: 风速
    :param y:
    :return:
    '''
    c7 = p[0]
    c8 = p[1]
    c9 = p[2]

    c1 = CoeffA[0]
    c2 = CoeffA[1]
    c3 = CoeffA[2]
    c4 = CoeffB[0]
    c5 = CoeffB[1]
    c6 = CoeffB[2]

    A = c1 + c2 * x2 + c3 * x2**2
    B = c4 + c5 * x2 + c6 * x2**2

    # return (A + (B - A) * np.exp(((c9 - 75.0)**2 + (x1 - c9)**2)/ (c7 + c8 * x2))) - y
    return ((c1 + c2 * x2 + c3 * x2**2) + \
            ((c4 + c5 * x2 + c6 * x2**2) - (c1 + c2 * x2 + c3 * x2**2))  \
            * np.exp(((c9 - 75.0)**2 - (x1 - c9)**2)/ (c7 + c8 * x2))) - y

def least_squareD(p, x1, x2, y):
    '''
    最小二乘法残差计算c7、c8、c9
    :param p:
    :param x1: 卫星天顶角
    :param x2: 风速
    :param y:
    :return:
    '''
    d1 = p[0]
    d2 = p[1]


    c1 = CoeffA[0]
    c2 = CoeffA[1]
    c3 = CoeffA[2]
    c4 = CoeffB[0]
    c5 = CoeffB[1]
    c6 = CoeffB[2]
    c7 = CoeffC[0]
    c8 = CoeffC[1]
    c9 = CoeffC[2]


    # A = c1 + c2 * x2 + c3 * x2**2
    # B = c4 + c5 * x2 + c6 * x2**2

    return ((c1 + c2 * x2 + c3 * x2**2) + \
            ((c4 + c5 * x2 + c6 * x2**2) - (c1 + c2 * x2 + c3 * x2**2))  \
            * np.exp(((c9 - 75.0)**2 - (x1 - c9)**2)/ (c7 + c8 * x2)))  \
            + (279.0 - 301.2)*(d1 + np.exp((d2 * x1**2) / 75.0**2 )) - y


def least_square1(x, c7, c8, c9 ):
    '''
    最小二乘法残差计算c7、c8、c9
    :param p:
    :param x1: 卫星天顶角
    :param x2: 风速
    :param y:
    :return:
    '''
    # c7 = p[0]
    # c8 = p[1]
    # c9 = p[2]

    c1 = CoeffA[0]
    c2 = CoeffA[1]
    c3 = CoeffA[2]
    c4 = CoeffB[0]
    c5 = CoeffB[1]
    c6 = CoeffB[2]

    x1 = x[0, :]
    x2 = x[1, :]

    # A = c1 + c2 * x2 + c3 * x2 ** 2
    # B = c4 + c5 * x2 + c6 * x2 ** 2

    # return (A + (B - A) * np.exp(((c9 - 75.0)**2 + (x1 - c9)**2)/ (c7 + c8 * x2))) - y
    return (c1 + c2 * x2 + c3 * x2 ** 2) + \
           ((c4 + c5 * x2 + c6 * x2 ** 2) - (c1 + c2 * x2 + c3 * x2 ** 2)) \
           * np.exp(((c9 - 75.0) ** 2 - (x1 - c9) ** 2) / (c7 + c8 * x2))
    # return (c1 + c2 * x2 + c3 * x2**2) + \
    #         ((c4 + c5 * x2 + c6 * x2**2) - (c1 + c2 * x2 + c3 * x2**2)) \
    #         * np.exp(((75.0 - x1)*(75.0 + x1 - 2 * c9)) / (c7 + c8 * x2))


def least_square2(x, coeff ):
    '''
    最小二乘法残差计算c7、c8、c9
    :param p:
    :param x1: 卫星天顶角
    :param x2: 风速
    :param y:
    :return:
    '''


    c1 = coeff[0]
    c2 = coeff[1]
    c3 = coeff[2]
    c4 = coeff[3]
    c5 = coeff[4]
    c6 = coeff[5]
    c7 = coeff[6]
    c8 = coeff[7]
    c9 = coeff[8]

    d1 = coeff[9]
    d2 = coeff[10]

    x1 = x[0, :]
    x2 = x[1, :]

    # A = c1 + c2 * x2 + c3 * x2 ** 2
    # B = c4 + c5 * x2 + c6 * x2 ** 2

    # return (A + (B - A) * np.exp(((c9 - 75.0)**2 + (x1 - c9)**2)/ (c7 + c8 * x2))) - y
    return (c1 + c2 * x2 + c3 * x2 ** 2) + \
           ((c4 + c5 * x2 + c6 * x2 ** 2) - (c1 + c2 * x2 + c3 * x2 ** 2)) \
           * np.exp(((c9 - 75.0) ** 2 - (x1 - c9) ** 2) / (c7 + c8 * x2)) \
           + (279.0 - 301.2) * (d1 + np.exp((d2 * x1 ** 2) / 75.0 ** 2))


def Emiss_Integrate_Simpson(wavenumber, response, emissivity, count):
    '''
    辛普森复式卷积法
    :param wavenumber:
    :param emissivity:
    :param response:
    :param count:
    :return:
    '''

    sum1 = 0
    sum2 = 0
    for i in range(int((count-1) / 2.0)):
        x1 = i * 2
        x2 = i * 2 + 1
        x3 = (i + 1) * 2
        deltx = wavenumber[x3] - wavenumber[x1]

        sum1 += (deltx / 6) * (emissivity[x1]*response[x1] \
                                + 4 * emissivity[x2]*response[x2]
                                + emissivity[x3]*response[x3])
        sum2 += (deltx / 6) * (response[x1] + 4*response[x2] + response[x3])

    if sum2 == 0:
        return np.nan
    else:
        return sum1 / sum2

def Emiss_Integrate_Trapezoid(wavenumber, response, emissivity, count):
    '''
    梯形卷积法
    :param wavenumber:
    :param emissivity:
    :param response:
    :param count:
    :return:
    '''

    sum1 = 0
    sum2 = 0
    for i in range(count-1):
        x1 = i
        x2 = i + 1
        deltx = wavenumber[x2] - wavenumber[x1]
        sum1 += (deltx / 2) * (emissivity[x1]*response[x1] \
                                + emissivity[x2]*response[x2])
        sum2 += (deltx / 2) * (response[x1] + response[x2])

    if sum2 == 0:
        return np.nan
    else:
        return  sum1 / sum2



def AnalysisTxt(filename, dtype=np.float64):
    '''
    读取TXT文本文件
    :param filename:
    :return:
    '''
    with open(filename, 'r') as fp :
        tmp = fp.readlines()
        fp.close()
    data = []
    for item in tmp[1:] :   # 剔除标题行
        tmpdata = item.split()
        # tmpdata = np.array(tmpdata, dtype= np.float64)
        data.append(tmpdata)

    return np.array(data, dtype=dtype)

def ZenithEmissiCoeff(path_SeaEmiss, wavenumber, response, zenith, Emissi_Flag = False, Coeff_Flag = False):
    '''
    计算不同卫星天顶角海表发射率拟合系数
    :param wavenumber:
    :param response:
    :param zenith:
    :param Emissi_Flag:
    :param Coeff_Flag:
    :return:
    '''

    WindSpeed = np.arange(RTTOV_Param['Wind_Range'][0], RTTOV_Param['Wind_Range'][1]+0.1, RTTOV_Param['Wind_Step'])

    # 按不同角度获取发射率
    Data_Emissity = []
    for iwind in WindSpeed:

        EmissFile = os.path.join(path_SeaEmiss, 'IREMIS_Emissivity_ws_%d' %(int(iwind)),
                                 'emissvity_zenith_%.2f.txt' % (zenith))
        if not os.path.isfile(EmissFile):
               raise Exception(EmissFile + ' is not exist, please check this file exist...')
        data = AnalysisTxt(EmissFile)
        wave = data[:, 0]
        emiss = data[:, 1]

        index = np.where((wave >= wavenumber[0] - WaveNumStep) & (wave <= wavenumber[-1] + WaveNumStep))
        wavetmp = wave[index]
        emisstmp = emiss[index]

        # 对发射率进行线性插值
        fit2 = interpolate.interp1d(wavetmp, emisstmp, kind='slinear')
        emission = fit2(wavenumber)

        # 发射率与光谱响应之间卷积
        tmpdata = Emiss_Integrate_Simpson(wavenumber, response, emission, wavenumber.shape[0])
        # data = Emiss_Integrate_Trapezoid(wave, resp, emiss, wave.shape[0])

        Data_Emissity.append(tmpdata)
    Data_Emissity = np.array(Data_Emissity, dtype=np.float64)

    # 非线性拟合
    # a = optimize.curve_fit(func, SatZ/60.0, Data_Emissity)[0]

    if Coeff_Flag :
        # 最小二乘法拟合
        coeff0 = [1.0, 0, 0]  # 初始值
        coeff = optimize.leastsq(least_square_w, coeff0, args=(WindSpeed, Data_Emissity))[0]
        # print('Emissity coeff:', coeff)

    if Emissi_Flag and Coeff_Flag :
        return Data_Emissity, coeff
    elif Emissi_Flag :
        return Data_Emissity
    elif Coeff_Flag :
        return coeff
    else:
        raise Exception("Calculate ")


def SeaEmissCoeff(wave, resp, CentreWaveNumber):
    '''
    计算地表发射率系数
    :param GPFuncTxt: 光谱系数文件
    :param CentreWaveNumber:通道中心波数
    :return:
    '''

    SatZenith = np.arange(RTTOV_Param['Zenith_Range'][0], RTTOV_Param['Zenith_Range'][1] + 0.1,
                          RTTOV_Param['Zenith_Step'])
    WindSpeed = np.arange(RTTOV_Param['Wind_Range'][0], RTTOV_Param['Wind_Range'][1] + 0.1, RTTOV_Param['Wind_Step'])

    # CentreWaveNumber = 10000.0 / 7.1
    print('CentreWaveNumber:', CentreWaveNumber)

    # 统一波数，对光谱响应进行线性插值
    wavenumber = np.arange(np.ceil(np.nanmin(wave)), np.floor(np.nanmax(wave)), WaveNumStep)
    # wavenumber = np.arange(np.nanmin(wave), np.nanmax(wave), WaveNumStep)
    fit1 = interpolate.interp1d(wave, resp, kind='slinear')
    response = fit1(wavenumber)

    global CoeffA, CoeffB, CoeffC, CoeffD

    # 计算最大、最小天顶角、所有风速范围、参考海表温度条件下的发射率
    EmissA, CoeffA = ZenithEmissiCoeff(Path_Sea_Emissivity_T0, wavenumber, response, 0.0, Emissi_Flag=True, Coeff_Flag=True)
    EmissB, CoeffB = ZenithEmissiCoeff(Path_Sea_Emissivity_T0, wavenumber, response, 75.0, Emissi_Flag=True, Coeff_Flag=True)
    print('Emissity CoeffA:', CoeffA)
    print('Emissity CoeffB:', CoeffB)


    Emiss_T0 = []
    Emiss_T1 = []
    for isatz in SatZenith :
        # 计算T0=301.2K的海表发射率
        m_Emiss = ZenithEmissiCoeff(Path_Sea_Emissivity_T0, wavenumber, response, isatz, Emissi_Flag=True)
        Emiss_T0.append(m_Emiss)

        # 计算T1=279
        m_Emiss = ZenithEmissiCoeff(Path_Sea_Emissivity_T1, wavenumber, response, isatz, Emissi_Flag=True)
        Emiss_T1.append(m_Emiss)

    Emiss_T0 = np.array(Emiss_T0, dtype=np.float64)
    Emiss_T1 = np.array(Emiss_T1, dtype=np.float64)
    WindSpeed = np.array(WindSpeed, dtype=np.float64)
    SatZenith = np.array(SatZenith, dtype=np.float64)
    # print(Emiss_All.shape)
    # print(Emiss_All)

    # SatZenith = SatZenith * np.pi / 180.0
    WS, SatZ = np.meshgrid(WindSpeed, SatZenith)

    xdata = []
    xdata.append(SatZ.reshape(-1))
    xdata.append(WS.reshape(-1))
    xdata = np.array(xdata)
    ydata = Emiss_T0.reshape(-1)

    # 初始值
    coeff0 = [1900.0, 30.0, 150.0]

    # CoeffC = optimize.curve_fit(least_square1, xdata, ydata, p0=coeff0)[0]
    # print('Emissity CoeffC:', CoeffC)
    CoeffC = optimize.leastsq(least_squareC, coeff0, args=(SatZ.reshape(-1), WS.reshape(-1), Emiss_T0.reshape(-1)))[0]
    print('Emissity CoeffC:', CoeffC)


    # 对于光谱在770 − 1230cm-1区间的通道，需要考虑温度的影响
    if np.nanmin(wave) >= 770 and np.nanmin(wave) <= 1230 :
        # 初始值
        coeff1 = [0.0, 0.0]
        CoeffD = optimize.leastsq(least_squareD, coeff1, args=(SatZ.reshape(-1), WS.reshape(-1), Emiss_T1.reshape(-1)))[0]
    elif np.nanmax(wave) >= 770 and np.nanmax(wave) <= 1230 :
        # 初始值
        coeff1 = [0.0, 0.0]
        CoeffD = optimize.leastsq(least_squareD, coeff1, args=(SatZ.reshape(-1), WS.reshape(-1), Emiss_T1.reshape(-1)))[0]
    else:
        CoeffD = [0.0, 0.0]
    print('Emissity CoeffD:', CoeffD)


    if DEBUG_TEST :
        # 此模块做验证系数的精度
        y = least_square1(xdata, CoeffC[0], CoeffC[1], CoeffC[2])

        coef = [0.98099368E+00, -0.49583451E-05, -0.24431857E-06,
                0.81908622E+00,  0.57013733E-02, -0.11514579E-03,
                0.19814249E+04,  0.32061268E+02, 0.15502974E+03,
                -0.10000010E+01,  0.60428910E-05]
        y1 = least_square2(xdata, coef)
        data = np.array(y - y1).reshape((31, 21))

        print(data)

        print(np.where(np.abs(data) > 0.005))

    Coeff = np.full(shape=(11,), fill_value=FILL_VALUE, dtype=np.float64)

    Coeff[0:3] = CoeffA[:]
    Coeff[3:6] = CoeffB[:]
    Coeff[6:9] = CoeffC[:]
    Coeff[9:11] = CoeffD[:]

    return Coeff





if __name__ == '__main__':

    for i in range(Dict_Chan_Info['BandTotalNum']):
        GPFuncTxt = os.path.join(Path_Response, 'FY4A_AGRI_B%02d_SRF.txt' % (i + 8))
        CentreWaveNumber = 10000.0 / Dict_Chan_Info['WaveLength'][i]


        if DEBUG_TEST:
            GPFuncTxt = os.path.join(Path_Response, 'rtcoef_fy4_1_agri_srf_ch%02d.txt' % (9))
            data = np.loadtxt(GPFuncTxt, dtype=np.float64)
            wave = data[:, 0]
            resp = data[:, 1]
        else:
            data = np.loadtxt(GPFuncTxt, dtype=np.float64)
            wave = 10000.0 / data[:, 0]
            resp = data[:, 1]

        coeff = SeaEmissCoeff(wave, resp, CentreWaveNumber)

        print(coeff)


