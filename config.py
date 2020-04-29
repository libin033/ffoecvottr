#coding:utf-8
import os
import sys


DEBUG_TEST = False

FILL_VALUE = -999.0

WaveNumStep = 0.25
AnglStep = 5.0

# 普朗克系数
Radiation_C1 = 1.191042953E-5
Radiation_C2 = 1.4387774

RTTOV_Param = {
    'RTTOV_Version': 12,        # rttov版本
    'Zenith_Range': [0.0, 75.0],# 卫星天顶角范围
    'Zenith_Step': 2.5,         # 卫星天顶角步长
    'Wind_Range': [0.0, 20.0],  # 风速范围
    'Wind_Step': 1.0,           # 风速步长
    'WaveNumStep': 0.25,        # 光谱分辨率
    'AnglStep': 5.0,            # 卫星天顶角步长
    'Ref_Zenith_Angle': 75.0,   # 参考角度
    'Ref_Tskin': 301.2,         # 参考温度
    'IR_SEA_EMIS_No': 11,       # 每个通道海表发射率系数个数
    'SatHeight': 35786.0,            # Nominal satellite height (km)
}


Dict_Chan_Info = {
    # 卫星载荷参数
    'Platform': 52,
    'Sat_ID' :  1,
    'Instrument' : 99,
    'BandTotalNum' : 7,
    'WaveLength' : [3.72, 6.25, 7.10, 8.50, 10.80, 12.00, 13.50],
    'Band_ID' : [8, 9, 10, 11, 12, 13, 14],
    'Planck_Weighted_Flag':[1, 0, 0, 0, 0, 0, 0],  # 标识哪些通道做普朗克权重校正


    'Dict_WaveNumber':{
        8:  10000.0/3.72,
        9:  10000.0/6.25,
        10: 10000.0/7.10,
        11: 10000.0/8.50,
        12: 10000.0/10.80,
        13: 10000.0/12.00,
        14: 10000.0/13.50
    },

    # 混合气体参数
    'Zenith':[0, 37, 48, 55, 60, 64],
    'Level': 53,
    'Point': 83,

    'profile_envelope':{
        'mixgas': 'Mixgas.txt',
        'ozone': 'Ozone.txt',
        'temperature': 'Temperature.txt',
        'water': 'Water.txt',
    },
    'reference_profile':{
        'mixgas': 'Mixgas.txt',
        'ozone': 'Ozone.txt',
        'water': 'Water.txt',
    },

    'Pressure_REF' : [0.005,   0.0131,   0.0304,    0.0644,    0.1263,    0.2324,    0.4052,    0.6749,
                 1.0801,   1.6691,   2.5011,    3.6462,    5.1864,     7.215,    9.8368,   13.1672,
                17.3308,  22.4601,  28.6937,   36.1735,    45.043,   55.4433,   67.5109,   81.3744,
                97.1505, 114.9415, 134.8318,  156.8846,  181.1394,  207.6092,  236.2784,  267.1012,
                   300., 334.8648, 371.5529,  409.8893,  449.6677,  490.6516,  532.5769,  575.1538,
               618.0706, 660.9965, 703.5863,  745.4841,  786.3278,  825.7546,  863.4047,  898.9275,
               931.9853, 962.2587, 989.4510, 1013.2923, 1033.5436,  1050. ],
}

Path_Home = r'D:\DZZ\data'

# 光谱响应函数
Path_Response = os.path.join(Path_Home, r'fy4a_agri_response')

# 发射率系数
Path_Emissity = os.path.join(Path_Home, r'ISEM_Emissivity_ws_0')

# 参考温度301.2K的海表发射率（20个风速）
Path_Sea_Emissivity_T0 = os.path.join(Path_Home, r'IREMIS_Emissivity_ws_301.2')
Path_Sea_Emissivity_T1 = os.path.join(Path_Home, r'IREMIS_Emissivity_ws_279')

# 透过率系数
# Path_Transmit_Miw = os.path.join(Path_Home, r'transmission/V1/Gas_miw')
# Path_Transmit_Mix = os.path.join(Path_Home, r'transmission/V1/Gas_mix')
# Path_Transmit_Mwo = os.path.join(Path_Home, r'transmission/V1/Gas_mwo')
# Path_Transmit_Gas = os.path.join(Path_Home, r'transmission/V1/Gas')
Path_Transmit_Miw = os.path.join(Path_Home, r'transmission/v2/')
Path_Transmit_Mix = os.path.join(Path_Home, r'transmission/v2/')
Path_Transmit_Mwo = os.path.join(Path_Home, r'transmission/v2/')

# Path_Transmit_Miw = os.path.join(Path_Home, r'transmission/v2_8/Gas_miw')
# Path_Transmit_Mix = os.path.join(Path_Home, r'transmission/v2_8/Gas_mix')
# Path_Transmit_Mwo = os.path.join(Path_Home, r'transmission/v2_8/Gas_mwo')

Path_Transmit_Gas = os.path.join(Path_Home, r'transmission/v2/Gas')

# 预报因子数据
Path_Predict_Water = os.path.join(Path_Home, r'predictors/water')
Path_Predict_Ozone = os.path.join(Path_Home, r'predictors/ozone')
Path_Predict_Mixgas = os.path.join(Path_Home, r'predictors/mixgas')

## 廓线和地面的83层数据
strProfFileName = os.path.join(Path_Home, 'profiles_83P.dat')
strSurfFileName = os.path.join(Path_Home, 'surface_83P.dat')




