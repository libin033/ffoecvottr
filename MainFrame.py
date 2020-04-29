#coding:utf-8
import os
import sys
import numpy as np
import datetime


exepath = os.path.dirname(__file__)
sys.path.append(exepath)

import PlanckWeighting
import BandCorrCoeff
import SurfEmisCoeff
import SeaEmissCoeff
from config import *
from TransmittancesCoeff import GetPredictorData, Cal_Trans_Coeff


def WriteRefFile(ofp, filename):
    if not os.path.isfile(filename) :
        raise Exception(filename + ' is not exist!!!')
    print(filename)
    fp = open(filename, 'r')
    lines = fp.readlines()
    fp.close()
    ofp.writelines(lines)

def CoeffProcess(strCoeffFile):
    NowTime = datetime.datetime.now()

    ofp = open(strCoeffFile, 'w')

    ofp.write(' ! RTTOV coefficient file fy4-1   agri-ir\n')
    ofp.write(' ! Automatic creation by subroutine rttov_write_ascii_coef\n')
    ofp.write(' ! RTTOV library version 12.1.0\n')
    ofp.write(' ! ------------------------------------------------------\n')
    ####################################################################################
    ################################# IDENTIFICATION ###################################
    ####################################################################################
    ofp.write('IDENTIFICATION\n')
    ofp.write(' !\n')
    ofp.write('%4d%4d%4d       ! Platform  sat_id  instrument\n'
              % (Dict_Chan_Info['Platform'], Dict_Chan_Info['Sat_ID'], Dict_Chan_Info['Instrument']))
    ofp.write(' fy4-1   agri-ir\n')
    ofp.write(' ir                ! Sensor type [ir,mw,hi,po]\n')
    ofp.write('%3d                ! RTTOV coefficient file version number\n' %(RTTOV_Param['RTTOV_Version']))
    ofp.write(' Created by rttov_lbl_make_coef.exe\n')
    ofp.write(' %s %s %s        ! Creation date\n'
              % (NowTime.strftime('%Y'), NowTime.strftime('%m'), NowTime.strftime('%d')))
    ofp.write(' ! ------------------------------------------------------\n')
    ####################################################################################
    ################################# LINE-BY-LINE #####################################
    ####################################################################################
    ofp.write('LINE-BY-LINE\n')
    ofp.write(' ! Line-by-line and other information\n')
    ofp.write(' !\n')
    ofp.write('LBLRTM_DB1/wv-ozone/l54/avg/fy4-1-agri-ir\n')
    ofp.write('LBLRTM\n')
    ofp.write('Created at SHK, September 2020\n')
    ofp.write('pascal.brunel@meteo.fr\n')
    ofp.write('Software = lblrtm_v12.2\n')
    ofp.write('Continuum = contnm.f(Revision: 16421) mt_ckd_2.5.2 (May 2020)\n')
    ofp.write('LBL data = lblrtm_v12.2/aer_v_3.2\n')
    ofp.write('Profiles = based on ECMWF_83P 2013 (marco.matricardi@ecmwf.int) revised for 1970-202x\n')
    ofp.write('  PROFILES_ECMWF_83_2016_CO2FIX 2 variable gases h2o, o3\n')
    ofp.write('  PROFILES_ECMWF_83_2016_CO2VAR 3 variable gases h2o, co2, o3\n')
    ofp.write('  PROFILES_ECMWF_83_2016_CH4VAR 6 variable gases h2o, co2, o3, n2o, co, ch4\n')
    ofp.write('  PROFILES_ECMWF_83_2016_SO2VAR 7 variable gases h2o, co2, o3, n2o, co, ch4, so2\n')
    ofp.write('  Minor gases from ECMWF_83P 2013 (marco.matricardi@ecmwf.int) and US76\n')
    ofp.write('Secants for LbL = 14\n')
    ofp.write('LbL spectral range 175-3300cm-1  LbL interpolated/averaged at 0.001cm-1\n')
    ofp.write('Spectral resolution for storage %4.2fcm-1\n' %(RTTOV_Param['WaveNumStep']))
    ofp.write('High resolution sounders, convolution at 0.001cm-1\n')
    ofp.write('Secant angles used for coefficient generation =\n')
    ofp.write('            0 37 48 55 60 64\n')
    ofp.write('Some channels were taken from rtcoef_fy4_1_agri-ir.H5--planck-weighted:\n')
    ofp.write(' 1\n')
    ofp.write(' ! ------------------------------------------------------\n')
    ####################################################################################
    ################################# FAST_MODEL_VARIABLES #############################
    ####################################################################################
    ofp.write('FAST_MODEL_VARIABLES\n')
    ofp.write(' !\n')
    ofp.write(' !\n')
    ofp.write(' RTTOV7            ! Fast model name\n')
    ofp.write('    7              ! Fast model version compatibility level\n')
    ofp.write('      7            ! Number of channels described in the coef file\n')
    ofp.write('    3              ! Number of gases described in the coef file\n')
    ofp.write('    0              ! PC compatibility level\n')
    ofp.write('    0              ! Zeeman flag\n')
    ofp.write(' Mixed_gases       ! Gas identification\n')
    ofp.write('   10  10  54      ! Variables/predictors  levels (pressure/absorber)\n')
    ofp.write(' Water_vapour      ! Gas identification\n')
    ofp.write('   15  15  54      ! Variables/predictors  levels (pressure/absorber)\n')
    ofp.write(' Ozone             ! Gas identification\n')
    ofp.write('   11  11  54      ! Variables/predictors  levels (pressure/absorber)\n')
    ofp.write(' ! ------------------------------------------------------\n')
    ####################################################################################
    #################### README_SPECTRAL_RESPONSE_FUNCTION #############################
    ####################################################################################
    ofp.write('README_SPECTRAL_RESPONSE_FUNCTION\n')
    ofp.write(' !\n')
    ofp.write('README file for FY4-1 AGRI filter functions  May 2020\n')
    ofp.write('Original file SRF_FY-4A_AGRI.zip provided by Qifeng Lu '
              'CMA luqf@cma.gov.cn the 16 April 2018 (email to Pascal Brunel)\n')
    ofp.write('Data are converted from micro-metres to wavenumber\n')
    ofp.write('Truncate to roughly +-3 half-width at half maximum to avoid very large wings\n')
    ofp.write(' ! ------------------------------------------------------\n')
    ####################################################################################
    ######################## FILTER_FUNCTIONS ##########################################
    ####################################################################################
    ofp.write('FILTER_FUNCTIONS\n')
    ofp.write(' !\n')
    ofp.write(' ! Channel number (from instrument original description)\n')
    ofp.write(' ! Channel status\n')
    ofp.write(' ! Central wavenumber\n')
    ofp.write(' ! Band correction coefficients (offset, slope)\n')
    ofp.write(' ! Gamma correction factor\n')
    # 计算波段校正系数Band correction coefficients
    for i in range(Dict_Chan_Info['BandTotalNum']):
        GPFuncTxt = os.path.join(Path_Response, 'FY4A_AGRI_B%02d_SRF.txt' % (i + 8))
        data = np.loadtxt(GPFuncTxt, dtype=np.float)
        wave = 10000 / data[:, 0]  # 波长转波数（cm-1）
        resp = data[:, 1]
        coeff = BandCorrCoeff.Cal_Band_Correction_Coefficients(wave, resp, 10000.0 / Dict_Chan_Info['WaveLength'][i])
        print('Band_Correction_Coefficients:', 10000.0 / Dict_Chan_Info['WaveLength'][i], coeff)
        # chanid    1   wavenumber  offset  slope   Gamma
        ofp.write('%6d%5d%19.10e%19.10e%19.10e%19.10e\n'
                  % (Dict_Chan_Info['Band_ID'][i], 1,
                     10000.0 / Dict_Chan_Info['WaveLength'][i], coeff[1], coeff[0], 1))
    ofp.write(' ! ------------------------------------------------------\n')
    ####################################################################################
    ######################### FUNDAMENTAL_CONSTANTS ####################################
    ####################################################################################
    ofp.write('FUNDAMENTAL_CONSTANTS\n')
    ofp.write(' !\n')
    ofp.write(' ! Units of constants for spectral radiance\n')
    ofp.write(' ! First radiation constant (mW/(m2.sr.cm-4))\n')
    ofp.write(' ! Second radiation constant (cm.K)\n')
    ofp.write('  1.191042953E-05  1.4387774 ! Planck constants\n')
    ofp.write('%11.1f                  ! Nominal satellite height (km)\n' %(RTTOV_Param['SatHeight']))
    ofp.write(' ! ------------------------------------------------------\n')
    ####################################################################################
    ################################# PLANCK_WEIGHTED ##################################
    ####################################################################################
    ofp.write('PLANCK_WEIGHTED\n')
    ofp.write(' !\n')
    ofp.write(' ! Channel number\n')
    ofp.write(' ! Planck-weighted flag (1 => yes; 0 => no)\n')
    for i in range(Dict_Chan_Info['BandTotalNum']):
        ofp.write('%5d%5d\n' % (i + 1, Dict_Chan_Info['Planck_Weighted_Flag'][i]))
    ofp.write(' ! ------------------------------------------------------\n')

    ####################################################################################
    ################################# SSIREM ###########################################
    ####################################################################################
    ofp.write('SSIREM\n')
    ofp.write(' !\n')
    ofp.write(' ! Channel number\n')
    ofp.write(' ! 5 coefficients for emissivity model SSIREM\n')
    ofp.write('  1   ! Version number\n')
    # 计算地表发射率系数
    for i in range(Dict_Chan_Info['BandTotalNum']):
        CentreWaveNumber = 10000.0 / Dict_Chan_Info['WaveLength'][i]

        GPFuncTxt = os.path.join(Path_Response, 'FY4A_AGRI_B%02d_SRF.txt' % (i + 8))
        data = np.loadtxt(GPFuncTxt, dtype=np.float64)
        wave = 10000.0 / data[:, 0]
        resp = data[:, 1]


        coeff = SurfEmisCoeff.SurfEmissCoeff(wave, resp, CentreWaveNumber)
        if CentreWaveNumber <= 750:
            N1 = 3
            N2 = 6
        else:
            N1 = 4
            N2 = 8
        ofp.write('%6d%12.7f%12.7f%12.7f%4.1f%4.1f\n'
                  % (i + 1, coeff[0], coeff[1], coeff[2], N1, N2))
    ofp.write(' ! ------------------------------------------------------\n')
    ####################################################################################
    ################################# IR_SEA_EMIS ######################################
    ####################################################################################
    ofp.write('IR_SEA_EMIS\n')
    ofp.write(' !\n')
    ofp.write('%3d   ! Version number\n' % (1))
    ofp.write('%3d   ! Number of coefficients per channel\n' % (RTTOV_Param['IR_SEA_EMIS_No']))
    ofp.write(' %8.3f%8.3f   ! Reference zenith angle and Tskin values\n'
              %(RTTOV_Param['Ref_Zenith_Angle'], RTTOV_Param['Ref_Tskin']))
    # 天顶角75度

    num = 0
    # 计算海表发射率系数
    for i in range(Dict_Chan_Info['BandTotalNum']):
        CentreWaveNumber = 10000.0 / Dict_Chan_Info['WaveLength'][i]

        GPFuncTxt = os.path.join(Path_Response, 'FY4A_AGRI_B%02d_SRF.txt' % (i + 8))
        data = np.loadtxt(GPFuncTxt, dtype=np.float64)
        wave = 10000.0 / data[:, 0]
        resp = data[:, 1]

        # 逐通道计算系数
        coeff = SeaEmissCoeff.SeaEmissCoeff(wave, resp, CentreWaveNumber)
        for i in range(11):
            ofp.write('%16.8E' % (coeff[i]))
            num += 1
            if num % 5 == 0:
                ofp.write('\n')
    if num % 5 != 0:
        ofp.write('\n')
    ofp.write(' ! ------------------------------------------------------\n')
    ####################################################################################
    ################################# REFERENCE_PROFILE 标准廓线##########################
    ####################################################################################
    ofp.write('REFERENCE_PROFILE\n')
    ofp.write(' !\n')
    ofp.write(' ! Reference pressure (hPa), reference temperature (K) and\n')
    ofp.write(' ! reference/background volume mixing ratios (ppmv) for each gas\n')
    ofp.write(' ! Note that mixing ratio is "missing" for mixed gases\n')
    ofp.write(' !     Mixed_gases\n')
    refprofile = os.path.join(Path_Home, 'reference_profile', Dict_Chan_Info['reference_profile']['mixgas'])
    WriteRefFile(ofp, refprofile)

    ofp.write(' !     Water_vapour\n')
    refprofile = os.path.join(Path_Home, 'reference_profile', Dict_Chan_Info['reference_profile']['water'])
    WriteRefFile(ofp, refprofile)

    ofp.write(' !     Ozone\n')
    refprofile = os.path.join(Path_Home, 'reference_profile', Dict_Chan_Info['reference_profile']['ozone'])
    WriteRefFile(ofp, refprofile)

    ofp.write(' ! ------------------------------------------------------\n')
    ####################################################################################
    ########################### PROFILE_ENVELOPE 包络线##################################
    ####################################################################################
    ofp.write('PROFILE_ENVELOPE\n')
    ofp.write(' !\n')
    ofp.write(' ! Reference pressure (hPa), temperature max and min (K) and\n')
    ofp.write(' ! volume mixing ratio max and min (ppmv) for each gas\n')
    ofp.write(' !      Temperature\n')
    refprofile = os.path.join(Path_Home, 'profile_envelope', Dict_Chan_Info['profile_envelope']['temperature'])
    WriteRefFile(ofp, refprofile)

    ofp.write(' !     Mixed_gases\n')
    refprofile = os.path.join(Path_Home, 'profile_envelope', Dict_Chan_Info['profile_envelope']['mixgas'])
    WriteRefFile(ofp, refprofile)

    ofp.write(' !     Water_vapour\n')
    refprofile = os.path.join(Path_Home, 'profile_envelope', Dict_Chan_Info['profile_envelope']['water'])
    WriteRefFile(ofp, refprofile)

    ofp.write(' !     Ozone\n')
    refprofile = os.path.join(Path_Home, 'profile_envelope', Dict_Chan_Info['profile_envelope']['ozone'])
    WriteRefFile(ofp, refprofile)

    ofp.write(' ! ------------------------------------------------------\n')
    ####################################################################################
    ################################# FAST_COEFFICIENTS#################################
    ####################################################################################
    Data_Pred_Water, Data_Pred_Ozone, Data_Pred_Mixgas = GetPredictorData()
    print('Data_Pred_Water.shape:', Data_Pred_Water.shape)
    print('Data_Pred_Ozone.shape:', Data_Pred_Ozone.shape)
    print('Data_Pred_Mixgas.shape:', Data_Pred_Mixgas.shape)
    ofp.write('FAST_COEFFICIENTS\n')
    ofp.write(' !\n')
    ofp.write(' ! Transmission coefficients\n')
    ofp.write(' ! Order of the gases:\n')
    ofp.write(' !     Mixed_gases\n')
    ofp.write(' !     Water_vapour\n')
    ofp.write(' !     Ozone\n')
    ofp.write('Mixed_gases\n')

    coeff_mixgas = Cal_Trans_Coeff('mixgas', Data_Pred_Mixgas)
    PreNum, ChanNum, Level = coeff_mixgas.shape
    coeff_mixgas = coeff_mixgas.reshape(-1)
    num = 0
    for i in range(PreNum*ChanNum*Level):
        ofp.write('%16.8E' %(coeff_mixgas[i]))
        num += 1
        if num % 5 == 0:
            ofp.write('\n')
    if num % 5 != 0:
        ofp.write('\n')

    ofp.write('Water_vapour\n')
    coeff_mixgas = Cal_Trans_Coeff('water', Data_Pred_Water)
    PreNum, ChanNum, Level = coeff_mixgas.shape
    coeff_mixgas = coeff_mixgas.reshape(-1)
    for i in range(PreNum * ChanNum * Level):
        ofp.write('%16.8E' % (coeff_mixgas[i]))
        if (i + 1) % 5 == 0:
            ofp.write('\n')

    ofp.write('Ozone\n')
    coeff_mixgas = Cal_Trans_Coeff('ozone', Data_Pred_Ozone)
    PreNum, ChanNum, Level = coeff_mixgas.shape
    coeff_mixgas = coeff_mixgas.reshape(-1)
    for i in range(PreNum * ChanNum * Level):
        ofp.write('%16.8E' % (coeff_mixgas[i]))
        if (i + 1) % 5 == 0:
            ofp.write('\n')
    if (i+1) % 5 != 0:
        ofp.write('\n')

    ofp.write(' ! ------------------------------------------------------\n')
    ofp.write('END\n')

    ofp.close()
    ####################################################################################
    #################################  END   ###########################################
    ####################################################################################


if __name__ == '__main__':

    # 传入输出系数文件名
    CoeffProcess('fy4a_agri.dat')