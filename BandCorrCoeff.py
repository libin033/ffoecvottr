#coding:utf-8

import os
import numpy as np
from scipy import optimize
from matplotlib import pyplot as plt
from scipy import interpolate
from config import *

def fx(temp, freq, resp):
    Radiation_C1 = 1.1910439e-5
    Radiation_C2 = 1.438769
    a = Radiation_C1 * freq * freq * freq
    b = (Radiation_C2 * freq) / temp
    c = np.exp(b) - 1.0
    radi = a / c
    return radi * resp

def integrate(temp, wavenumber, response, nstep):
    '''
    光谱积分
    :param temp:
    :param wavenumber:
    :param response:
    :param deltx:
    :param nstep:
    :return:
    '''

    a = wavenumber[0]
    b = wavenumber[-1]
    deltx = wavenumber[2] - wavenumber[0]
    # print(nstep, deltx)
    totals = 0.0
    sum1 = 0
    sum2 = 0
    for i in range(int((nstep-1) / 2.0)):
        x1 = i * 2
        x2 = i * 2 + 1
        x3 = (i + 1) * 2
        # print(x1, x2, x3)
        sum1 += (deltx / 6) * (fx(temp, wavenumber[x1], response[x1]) \
                                + 4 * fx(temp, wavenumber[x2],response[x2])
                                + fx(temp, wavenumber[x3], response[3]))

        sum2 += (deltx / 6) * (response[x1] + 4*response[x2] + response[3])

    totals = sum1 / sum2
    return totals


def Planck_t2r(temp, freq):
    '''
    Units: mW/(m2.cm-1.sr
    :param temp: bright temperature
    :param freq: wave number
    :return:
    '''
    Radiation_C1 = 1.191042953E-5
    Radiation_C2 = 1.4387774
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

    Radiation_C1 = 1.191042953E-5
    Radiation_C2 = 1.4387774
    vs = freq
    bt = (Radiation_C2*vs/np.log(Radiation_C1*vs*vs*vs/(rad)+1.0))

    return bt

def func(x, a, b):
    return (x - b) / a

def Process(GPFuncTxt):
    data = np.loadtxt(GPFuncTxt, dtype=np.float)
    wavelength = data[:, 0]
    response = data[:, 1]
    wavenumber = 10000.0 / wavelength
    print(wavenumber.shape)

    index = np.nanargmax(response)
    centerwavenumber = 0.1604836734E+04  # 10000/6.25
    print('centerwavenumber:', centerwavenumber)

    Te = []
    Num = wavenumber.shape[0]
    for BT in np.arange(180.0, 340.0001, 1.0):
        Rad = integrate(BT, wavenumber, response, STEP, Num)
        bt = Planck_r2t(Rad, centerwavenumber)

        Te.append(bt)
    Te = np.array(Te)
    BT = np.arange(180.0, 340.0001, 1.0)
    a = optimize.curve_fit(func, Te, BT)
    print(a)
    # print(a, b)
    #
    # plt.plot(Te, BT, 'b.')
    # y = (Te - b) / a
    #
    # plt.plot(Te, BT, 'r.')
    # # print(BT - y)
    #
    # plt.show()

def Cal_Band_Correction_Coefficients(wave, resp, centerwavenumber):
    '''

    :param wave:
    :param resp:
    :return:
    '''

    fit = interpolate.interp1d(wave, resp, kind= 'slinear')
    wavenumber = np.arange(np.ceil(np.nanmin(wave)), np.floor(np.nanmax(wave)), WaveNumStep / 2.0) #0.125)

    response = fit(wavenumber)
    # centerwavenumber = 0.1404847737E+04
    print('centerwavenumber:', centerwavenumber)

    Te = []
    Num = wavenumber.shape[0]
    for BT in np.arange(180.0, 340.0001, 1.0):
        Rad = integrate(BT, wavenumber, response, Num)
        bt = Planck_r2t(Rad, centerwavenumber)
        # print(bt)
        Te.append(bt)
    Te = np.array(Te)
    BT = np.arange(180.0, 340.0001, 1.0)
    coeff = optimize.curve_fit(func, Te, BT)[0]

    return coeff
    # print(a, b)
    # plt.plot(Te, BT, 'b.')
    # y = (Te - b) / a
    # plt.plot(Te, y, 'r.')
    # # print(BT - y)
    # plt.xlim(180, 340)
    # plt.ylim(180, 340)
    # plt.show()

    # z1 = np.polyfit(Te, BT, 1)
    # p1 = np.poly1d(z1)
    #
    # print(z1)
    # print(p1)

if __name__ == '__main__':
    # GPFuncTxt = os.path.join(r'D:\DZZ\仿真', 'FY4A_AGRI_B09_SRF.txt')
    # Process(GPFuncTxt)

    for i in range(Dict_Chan_Info['BandTotalNum']):
        GPFuncTxt = os.path.join(Path_Response, 'FY4A_AGRI_B%02d_SRF.txt' % (i + 8))
        data = np.loadtxt(GPFuncTxt, dtype=np.float64)
        wave = 10000 / data[:, 0]  # 波长转波数（cm-1）
        resp = data[:, 1]
        coeff = Cal_Band_Correction_Coefficients(wave, resp, 10000.0 / Dict_Chan_Info['WaveLength'][i])
        print(coeff)





