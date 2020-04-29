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
from config import *


if __name__ == '__main__':

    ifp =open('rtcoef_fy4_1_agri.dat', 'r')
    lines = ifp.readlines()
    ifp.close()

    num = 0

    for num in range(3604+1):
        item = lines[num]
        if 'REFERENCE_PROFILE\n' == item:
            refprofile = os.path.join(Path_Home, 'reference_profile', Dict_Chan_Info['reference_profile']['mixgas'])
            ofp = open(refprofile, 'w')
            starline = num + 6
            endline = starline + 54
            line = lines[starline:endline]
            ofp.writelines(line)
            ofp.close()
            num = endline
            refprofile = os.path.join(Path_Home, 'reference_profile', Dict_Chan_Info['reference_profile']['water'])
            ofp = open(refprofile, 'w')
            starline = num + 1
            endline = starline + 54
            line = lines[starline:endline]
            ofp.writelines(line)
            ofp.close()
            num = endline
            refprofile = os.path.join(Path_Home, 'reference_profile', Dict_Chan_Info['reference_profile']['ozone'])
            ofp = open(refprofile, 'w')
            starline = num + 1
            endline = starline + 54
            line = lines[starline:endline]
            ofp.writelines(line)
            ofp.close()
            num = endline

        print(num)
        if 'PROFILE_ENVELOPE\n' == item:
            print(num, item)
            refprofile = os.path.join(Path_Home, 'profile_envelope', Dict_Chan_Info['profile_envelope']['temperature'])
            ofp = open(refprofile, 'w')
            starline = num + 5
            endline = starline + 54
            line = lines[starline:endline]
            ofp.writelines(line)
            ofp.close()
            num = endline
            refprofile = os.path.join(Path_Home, 'profile_envelope', Dict_Chan_Info['profile_envelope']['mixgas'])
            ofp = open(refprofile, 'w')
            starline = num + 1
            endline = starline + 54
            line = lines[starline:endline]
            ofp.writelines(line)
            ofp.close()
            num = endline
            refprofile = os.path.join(Path_Home, 'profile_envelope', Dict_Chan_Info['profile_envelope']['water'])
            ofp = open(refprofile, 'w')
            starline = num + 1
            endline = starline + 54
            line = lines[starline:endline]
            ofp.writelines(line)
            ofp.close()
            num = endline
            refprofile = os.path.join(Path_Home, 'profile_envelope', Dict_Chan_Info['profile_envelope']['ozone'])
            ofp = open(refprofile, 'w')
            starline = num + 1
            endline = starline + 54
            line = lines[starline:endline]
            ofp.writelines(line)
            ofp.close()
            num = endline
        num += 1



