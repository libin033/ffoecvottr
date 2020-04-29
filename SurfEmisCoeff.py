#coding:utf-8
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


def func(x, a, b, c):

    return a - b * x**N1 - c * x**N2


def least_square(p, x, y):
    '''
    最小二乘法残差计算
    :param p:
    :param x:
    :param y:
    :return:
    '''
    a = p[0]
    b = p[1]
    c = p[2]

    return (a - b * x**N1 - c * x**N2) - y


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


def SurfEmissCoeff(wave, resp, CentreWaveNumber):
    '''
    计算地表发射率系数
    :param ChanID:通道ID
    :param CentreWaveNumber:通道中心波数
    :return:
    '''
    SatZ = np.arange(0.0, 60.01, AnglStep)

    # CentreWaveNumber = 10000.0 / 7.1
    print('CentreWaveNumber:', CentreWaveNumber)

    # GPFuncTxt = os.path.join(Path_Response, 'rtcoef_fy4_1_agri_srf_ch%02d.txt' %(ChanID))
    # data = np.loadtxt(GPFuncTxt, dtype=np.float64)
    # wave = data[:, 0]
    # resp = data[:, 1]

    # GPFuncTxt = os.path.join(r'D:\DZZ\Code\data\fy4a_agri_response', 'FY4A_AGRI_B10_SRF.txt')
    # data = np.loadtxt(GPFuncTxt, dtype=np.float64)
    # wave = 10000.0/data[:, 0]
    # resp = data[:, 1]

    # 统一波数，对光谱响应进行线性插值
    wavenumber = np.arange(np.ceil(np.nanmin(wave)), np.floor(np.nanmax(wave)), WaveNumStep)
    # wavenumber = np.arange(np.nanmin(wave), np.nanmax(wave), WaveNumStep)
    fit1 = interpolate.interp1d(wave, resp, kind='slinear')
    response = fit1(wavenumber)

    # 中心波数小于750cm-1,N1、N2分别取3、6，否则，取4、8
    global N1, N2
    if CentreWaveNumber <= 750:
        N1 = 3
        N2 = 6
    else:
        N1 = 4
        N2 = 8

    # 按不同角度获取发射率
    Data_Emissity = []
    for item in SatZ:
        EmissFile = os.path.join(Path_Emissity,
                                 'emissvity_zenith_%.2f.txt' % (item))
        if not os.path.isfile(EmissFile):
            continue
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

    # 最小二乘法拟合
    coeff0 = [1.0, 0, 0]  # 初始值
    coeff = optimize.leastsq(least_square, coeff0, args=(SatZ / 60.0, Data_Emissity))[0]
    print('Emissity coeff:', coeff)

    # 对处理结果进行分析
    # x = np.arange(0, 61, 5.0) / 60.0
    #
    # y = func(x, coeff[0], coeff[1], coeff[2])
    # p = [0.9814065, 0.0161783, 0.0222909]
    # y1 = func(x, p[0], p[1], p[2])
    #
    # plt.plot(x * 60, y, 'r.')
    # plt.plot(x * 60, y1, 'b.')
    # plt.plot(SatZ, Data_Emissity, 'g.')
    #
    # plt.show()
    # print(y1 -y)

    return coeff




if __name__ == '__main__':


    SurfEmissCoeff(9, 10000.0 / 7.1)


