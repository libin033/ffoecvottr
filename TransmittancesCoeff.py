#coding:utf-8
import os
import sys
import numpy as np
import datetime
from scipy import interpolate
from matplotlib import pyplot as plt
from scipy import optimize
exepath = os.path.dirname(__file__)
sys.path.append(exepath)
from config import *
import PlanckWeighting

def func_line(x, a, b, c):

    return a - b * x - c * x


def func_least_square_N(x, y, PredNum):
    # here, create lambda functions for Line fit
    # tpl is a tuple that contains the parameters of the fit
    funcLine = lambda tpl, x: np.dot(x, tpl)
    # func is going to be a placeholder for funcLine,funcQuad or whatever
    # function we would like to fit
    func = funcLine
    # ErrorFunc is the diference between the func and the y "experimental" data
    ErrorFunc = lambda tpl, x, y: func(tpl, x) - y
    # tplInitial contains the "first guess" of the parameters
    tplInitial = np.zeros(shape=(PredNum), dtype=np.float64)
    # leastsq finds the set of parameters in the tuple tpl that minimizes
    # ErrorFunc=yfit-yExperimental
    coeff, success = optimize.leastsq(ErrorFunc, tplInitial, args=(x, y))
    # print('coeff:', coeff)
    y1 = funcLine(coeff, x)

    # plt.plot(y1, y, 'r.')
    # # plt.xlim((0, 7))
    # # plt.ylim((0, 7))
    # # plt.plot(x, y, 'g.')
    # # plt.plot(x, y1, 'b.')
    # plt.show()
    # exit(0)
    return coeff

def Trans_Integrate_Simpson(wavenumber, response, transmission, count, prof_temp=1.0):
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

        # 普朗克权重：通过分层温度廓线计算对应层的辐亮度作为各层权重
        if prof_temp != 1.0 :
            planck_weight1 = PlanckWeighting.Planck_t2r(prof_temp, wavenumber[x1])
            planck_weight2 = PlanckWeighting.Planck_t2r(prof_temp, wavenumber[x2])
            planck_weight3 = PlanckWeighting.Planck_t2r(prof_temp, wavenumber[x3])
        else:
            planck_weight1 = 1.0
            planck_weight2 = 1.0
            planck_weight3 = 1.0

        sum1 += (deltx / 6) * (transmission[x1] * response[x1] * planck_weight1 \
                                + 4 * transmission[x2] * response[x2] * planck_weight2
                                + transmission[x3]*response[x3] * planck_weight3)
        sum2 += (deltx / 6) * (response[x1]* planck_weight1 \
                               + 4*response[x2]* planck_weight2 \
                               + response[x3]* planck_weight3)

        # print(x1, sum1, transmission[x1], transmission[x2], transmission[x3])

    if sum2 == 0:
        return np.nan
    else:
        return sum1 / sum2

def GasDataIntegral(wavenumber, response, Data_wn, Data_src, prof_temp):

    # 进行线性插值
    fit2 = interpolate.interp1d(Data_wn, Data_src, kind='slinear')
    Data_dst = fit2(wavenumber)

    # 发射率与光谱响应之间卷积
    tmpdata = Trans_Integrate_Simpson(wavenumber, response, Data_dst, wavenumber.shape[0], prof_temp)

    return tmpdata

def Cal_Water_Ozone_600_2600(wavenumber, response, chan_id):
    '''
    分通道进行处理
    输入混合气体的透过率，分别计算水汽、臭氧、混合气体的透过率
    :param wavenumber: 通道波数
    :param response: 通道光谱相应函数
    :param chan_id: 通道ID, 从0开始
    :return: 水汽、臭氧、混合气体的透过率（6个天顶角、83条线、53层）
    '''
    Ozone = np.full(shape=(6, 83, 53), fill_value=FILL_VALUE)
    Water = np.full(shape=(6, 83, 53), fill_value=FILL_VALUE)
    MixGas = np.full(shape=(6, 83, 53), fill_value=FILL_VALUE)


    wave = np.arange(600.0, 2600.0, WaveNumStep)
    index = np.where((wave > wavenumber[0] - WaveNumStep) & (wave < wavenumber[-1] + WaveNumStep))
    wave = wave[index]

    if Dict_Chan_Info['Planck_Weighted_Flag'][chan_id] == 1:
        # 获取53层的温度廓线
        Prof_T = PlanckWeighting.PlanckWeighting()
    else:
        Prof_T = np.ones(shape=(53,), dtype=np.float32)

    # 6个天顶角、83条线、53层
    for i in range(6) :
        izenith = Dict_Chan_Info['Zenith'][i]
        for ipoint in range(Dict_Chan_Info['Point']):
            for ilevel in range(Dict_Chan_Info['Level']):
                MiwFileName = os.path.join(Path_Transmit_Miw, '600_2600/Gas_miw',
                                           'zenith_%02d' %(izenith),
                                           'point_%02d' %(ipoint+1),
                                           'level_%03d.dat' %(ilevel+1))
                MixFileName = os.path.join(Path_Transmit_Mix, '600_2600/Gas_mix',
                                           'zenith_%02d' % (izenith),
                                           'point_%02d' % (ipoint+1),
                                           'level_%03d.dat' % (ilevel+1))
                MwoFileName = os.path.join(Path_Transmit_Mwo, '600_2600/Gas_mwo',
                                           'zenith_%02d' % (izenith),
                                           'point_%02d' % (ipoint+1),
                                           'level_%03d.dat' % (ilevel+1))
                print('start to do ',izenith, ipoint, ilevel)
                if not os.path.isfile(MiwFileName):
                    raise Exception(MiwFileName + ' is not exist, please check it')

                if not os.path.isfile(MixFileName):
                    raise Exception(MixFileName + ' is not exist, please check it')

                if not os.path.isfile(MwoFileName):
                    raise Exception(MwoFileName + ' is not exist, please check it')

                Data_Miw = np.loadtxt(MiwFileName, dtype=np.float64)
                Data_Mix = np.loadtxt(MixFileName, dtype=np.float64)
                Data_Mwo = np.loadtxt(MwoFileName, dtype=np.float64)

                Data_Miw = Data_Miw[index]
                Data_Mix = Data_Mix[index]
                Data_Mwo = Data_Mwo[index]

                # 剔除数据中0值
                ind = np.where( (Data_Miw != 0) & (Data_Mix != 0) & (Data_Mwo != 0))
                Data_Miw = Data_Miw[ind]
                Data_Mix = Data_Mix[ind]
                Data_Mwo = Data_Mwo[ind]
                wavetmp = wave[ind]

                # 计算水汽、臭氧的透过率
                # Data_Water = Data_Miw / Data_Mix
                # Data_Ozone = Data_Mwo / Data_Miw

                Data_Water = Data_Miw
                Data_Ozone = Data_Mwo



                # 对水汽、臭氧、混合气体进行卷积
                Data_Miw = GasDataIntegral(wavenumber, response, wavetmp, Data_Water, Prof_T[ilevel])
                Data_Mix = GasDataIntegral(wavenumber, response, wavetmp, Data_Mix, Prof_T[ilevel])
                Data_Ozone = GasDataIntegral(wavenumber, response, wavetmp, Data_Ozone, Prof_T[ilevel])

                # print(Data_Miw, Data_Mix, Data_Ozone)
                Water[i, ipoint, ilevel] = Data_Miw
                Ozone[i, ipoint, ilevel] = Data_Ozone
                MixGas[i, ipoint, ilevel] = Data_Mix


    return Water, Ozone, MixGas


def Cal_Water_Ozone_2000_4000(wavenumber, response, chan_id):
    '''
    进行处理8通道
    输入混合气体的透过率，分别计算水汽、臭氧、混合气体的透过率
    :param wavenumber: 通道波数
    :param response: 通道光谱相应函数
    :param chan_id: 通道ID, 从0开始
    :return: 水汽、臭氧、混合气体的透过率（6个天顶角、83条线、53层）
    '''
    Ozone = np.full(shape=(6, 83, 53), fill_value=FILL_VALUE)
    Water = np.full(shape=(6, 83, 53), fill_value=FILL_VALUE)
    MixGas = np.full(shape=(6, 83, 53), fill_value=FILL_VALUE)


    wave = np.arange(2000.0, 4000.0, WaveNumStep)
    index = np.where((wave > wavenumber[0] - WaveNumStep) & (wave < wavenumber[-1] + WaveNumStep))
    wave = wave[index]

    if Dict_Chan_Info['Planck_Weighted_Flag'][chan_id] == 1:
        # 获取53层的温度廓线
        Prof_T = PlanckWeighting.PlanckWeighting()
    else:
        Prof_T = np.ones(shape=(53,), dtype=np.float32)

    # 6个天顶角、83条线、53层
    for i in range(6) :
        izenith = Dict_Chan_Info['Zenith'][i]
        for ipoint in range(Dict_Chan_Info['Point']):
            for ilevel in range(Dict_Chan_Info['Level']):
                MiwFileName = os.path.join(Path_Transmit_Miw, '2000_4000/Gas_miw',
                                           'zenith_%02d' %(izenith),
                                           'point_%02d' %(ipoint+1),
                                           'level_%03d.dat' %(ilevel+1))
                MixFileName = os.path.join(Path_Transmit_Miw, '2000_4000/Gas_mix',
                                           'zenith_%02d' % (izenith),
                                           'point_%02d' % (ipoint+1),
                                           'level_%03d.dat' % (ilevel+1))
                MwoFileName = os.path.join(Path_Transmit_Miw, '2000_4000/Gas_mwo',
                                           'zenith_%02d' % (izenith),
                                           'point_%02d' % (ipoint+1),
                                           'level_%03d.dat' % (ilevel+1))
                print('start to do ',izenith, ipoint, ilevel)
                if not os.path.isfile(MiwFileName):
                    raise Exception(MiwFileName + ' is not exist, please check it')

                if not os.path.isfile(MixFileName):
                    raise Exception(MixFileName + ' is not exist, please check it')

                if not os.path.isfile(MwoFileName):
                    raise Exception(MwoFileName + ' is not exist, please check it')

                Data_Miw = np.loadtxt(MiwFileName, dtype=np.float64)
                Data_Mix = np.loadtxt(MixFileName, dtype=np.float64)
                Data_Mwo = np.loadtxt(MwoFileName, dtype=np.float64)

                Data_Miw = Data_Miw[index]
                Data_Mix = Data_Mix[index]
                Data_Mwo = Data_Mwo[index]

                # 剔除数据中0值
                ind = np.where( (Data_Miw != 0) & (Data_Mix != 0) & (Data_Mwo != 0))
                Data_Miw = Data_Miw[ind]
                Data_Mix = Data_Mix[ind]
                Data_Mwo = Data_Mwo[ind]
                wavetmp = wave[ind]

                # 计算水汽、臭氧的透过率
                # Data_Water = Data_Miw / Data_Mix
                # Data_Ozone = Data_Mwo / Data_Miw

                Data_Water = Data_Miw
                Data_Ozone = Data_Mwo

                # 对水汽、臭氧、混合气体进行卷积
                Data_Miw = GasDataIntegral(wavenumber, response, wavetmp, Data_Water, Prof_T[ilevel])
                Data_Mix = GasDataIntegral(wavenumber, response, wavetmp, Data_Mix, Prof_T[ilevel])
                Data_Ozone = GasDataIntegral(wavenumber, response, wavetmp, Data_Ozone, Prof_T[ilevel])

                # print(Data_Miw, Data_Mix, Data_Ozone)
                Water[i, ipoint, ilevel] = Data_Miw
                Ozone[i, ipoint, ilevel] = Data_Ozone
                MixGas[i, ipoint, ilevel] = Data_Mix


    return Water, Ozone, MixGas



def GetGasTrans():
    '''
    读取光谱响应函数和水汽、臭氧、混合气体的透过率
    :return: 输出水汽、臭氧、混合气体的透过率到临时npy
    '''
    for i in range(Dict_Chan_Info['BandTotalNum']):
        chan_id = i + 8
        GPFuncTxt = os.path.join(Path_Response, 'FY4A_AGRI_B%02d_SRF.txt' % (chan_id))
        CentreWaveNumber = 10000.0 / Dict_Chan_Info['WaveLength'][i]
        print(GPFuncTxt)
        data = np.loadtxt(GPFuncTxt, dtype=np.float64)
        wave = 10000.0 / data[:, 0]
        resp = data[:, 1]

        # 统一波数，对光谱响应进行线性插值
        wavenumber = np.arange(np.ceil(np.nanmin(wave)), np.floor(np.nanmax(wave)), WaveNumStep)
        # wavenumber = np.arange(np.nanmin(wave), np.nanmax(wave), WaveNumStep)
        fit1 = interpolate.interp1d(wave, resp, kind='slinear')
        response = fit1(wavenumber)

        # if chan_id == 8:
        #     Water, Ozone, MixGas = Cal_Water_Ozone_8(wavenumber, response)
        # else:
        #     Water, Ozone, MixGas = Cal_Water_Ozone(wavenumber, response)

        if wavenumber[0] >= 600 and wavenumber[-1] <= 2600:
            Water, Ozone, MixGas = Cal_Water_Ozone_600_2600(wavenumber, response, i)
        elif wavenumber[0] >= 2000 and wavenumber[-1] <= 4000:
            Water, Ozone, MixGas = Cal_Water_Ozone_2000_4000(wavenumber, response, i)
        else:
            raise Exception("Wavenumber [%d, %d] is not in [600, 2600] or [2000, 6000] " %(wavenumber[0], wavenumber[-1]))

        np.save(os.path.join(Path_Transmit_Gas, 'CH%02d_Water.npy' % (chan_id)), Water)
        np.save(os.path.join(Path_Transmit_Gas, 'CH%02d_Ozone.npy' % (chan_id)), Ozone)
        np.save(os.path.join(Path_Transmit_Gas, 'CH%02d_MixGas.npy' % (chan_id)), MixGas)



def LoadPredictorFile(filename):
    '''
    读取预报因子数据
    :param filename:预报因子文件名称
    :return:
    '''
    if not os.path.isfile(filename):
        raise Exception(filename + ' is not exist, please check it')
    # print(filename)
    fp = open(filename, 'r')

    lines = fp.readlines()
    point = Dict_Chan_Info['Point']
    level = Dict_Chan_Info['Level'] + 1
    data_r = []
    for i in range(point):
        starline = i * level+1
        endline = (i+1) * level
        tmplines = lines[starline:endline]
        data = []
        for item in tmplines :
            data.append(item.split())

        data_r.append(data)
    fp.close()

    return np.array(data_r, dtype=np.float64)

def GetPredictorData():
    Data_Pred_Water = []
    Data_Pred_Ozone = []
    Data_Pred_Mixgas = []

    for i in range(6):
        izenith = Dict_Chan_Info['Zenith'][i]
        WaterFileName = os.path.join(Path_Predict_Water, 'predictor_water_zenith_%02d.txt' % (izenith))
        OzoneFileName = os.path.join(Path_Predict_Ozone, 'predictor_ozone_zenith_%02d.txt' % (izenith))
        MixgasFileName = os.path.join(Path_Predict_Mixgas, 'predictor_mixgas_zenith_%02d.txt' % (izenith))

        predictor_water = LoadPredictorFile(WaterFileName)
        predictor_ozone = LoadPredictorFile(OzoneFileName)
        predictor_mixgas = LoadPredictorFile(MixgasFileName)

        Data_Pred_Water.append(predictor_water)
        Data_Pred_Ozone.append(predictor_ozone)
        Data_Pred_Mixgas.append(predictor_mixgas)

    Data_Pred_Water = np.array(Data_Pred_Water, dtype=np.float64)
    Data_Pred_Ozone = np.array(Data_Pred_Ozone, dtype=np.float64)
    Data_Pred_Mixgas = np.array(Data_Pred_Mixgas, dtype=np.float64)

    print(Data_Pred_Water.shape)
    print(Data_Pred_Ozone.shape)
    print(Data_Pred_Mixgas.shape)

    return Data_Pred_Water, Data_Pred_Ozone, Data_Pred_Mixgas

def Cal_Trans_Coeff(pretype, Data_Pred):
    '''
    计算透过率系数
    :return:
    '''

    PredNum = Data_Pred.shape[-1]
    ChanNum = Dict_Chan_Info['BandTotalNum']
    Level = Dict_Chan_Info['Level']

    Coeff_R = np.full(shape=(PredNum, ChanNum, Level), fill_value=-999.0, dtype=np.float64)

    for ichan in range(ChanNum): # 逐通道循环
        filename = os.path.join(Path_Transmit_Gas, 'CH%02d_%s.npy' % (ichan + 8, pretype))
        if not os.path.isfile(filename):
            print(filename + ' is not exist, will continue!!!')
            continue
        Water = np.load(filename)
        Data_Optical_Depth = np.log(Water) * -1.0

        deltax = np.full_like(Water, fill_value=-999.0, dtype=np.float64)
        deltax[:,:, 0] = Data_Optical_Depth[:, :, 0]
        deltax[:, :, 1:] = Data_Optical_Depth[:, :, 1:] - Data_Optical_Depth[:,:,:-1]

        for ilevel in range(Dict_Chan_Info['Level']) : # 逐层循环
            y = deltax[:, :, ilevel]
            x = Data_Pred[:, :, ilevel, :]

            x = x.reshape((-1, PredNum))
            y = y.reshape(-1)
            coef = func_least_square_N(x, y, PredNum)
            Coeff_R[:, ichan, ilevel] = coef



    print(Coeff_R.shape)
    return Coeff_R



if __name__ == '__main__':

    # 读取分成透过率文件并输出npy
    GetGasTrans()
    exit(0)
    Data_Pred_Water, Data_Pred_Ozone, Data_Pred_Mixgas = GetPredictorData()
    print('Data_Pred_Water.shape:', Data_Pred_Water.shape)
    print('Data_Pred_Ozone.shape:', Data_Pred_Ozone.shape)
    print('Data_Pred_Mixgas.shape:', Data_Pred_Mixgas.shape)
    Cal_Trans_Coeff('water', Data_Pred_Water)
    # Cal_Trans_Coeff('ozone', Data_Pred_Ozone)
    # Cal_Trans_Coeff('mixgas', Data_Pred_Mixgas)







