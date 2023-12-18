import numpy as np
import math
import pandas as pd
from decimal import *


def sampleCP(baseconf="coords.xyz", template="g16temp.gjf", scan="outfile.gjf", level='level', q2_min='q2min', q2_max='q2max', points='points'):

    getcontext().prec = 6
    ###########files########

    N = 5
    m = 2

    baseconf = pd.read_csv(conf_file, skiprows=0, engine='python', header=None, delim_whitespace=True)
    df2 = baseconf.rename({0:'atom',1:'X', 2:'Y', 3:'Z'}, axis=1)

    x_conv = df2['X'].values
    y_conv = df2['Y'].values
    z_conv = df2['Z'].values
    z_C = []
    for c in range(0, N):
        print(z_conv[c])
        z_C.append(z_conv[c])
    print(z_conv, type(z_conv))
    coords=[]
    CN_coords = {}
    r_i =[]
    sin_j=[]
    sin_jm = []
    cos_j=[]
    cos_jm=[]
    r_sys=[]
    a_m = []

    for i in range(0, N):
        jm = math.sin(2 * math.pi * i * m / N)
        km = math.cos(2 * math.pi * i * m / N)
        a = (2*math.pi*m*i)/N
        print(i, a)
        sin_jm.append(jm)
        cos_jm.append(km)
        a_m.append(a)

    ####Calculate CP coordinates
    q_m = (2/N)**(0.5)*(sum((z_C* np.asarray(cos_jm))**2) + sum((z_C* np.asarray(sin_jm))**2))
    theta_m_degs = math.degrees(math.atan(sum(z_C* np.asarray(sin_jm)))/(sum(z_C* np.asarray(cos_jm))))
    if math.cos(theta_m_degs)<0.0:
        theta_m_degs=math.degrees(theta_m_degs+180)


    print("Cremer-Pople coordinates at m=", m, "with q_m=", q_m, "and theta_m", theta_m_degs, "in degrees")

    dtheta_m = []
    dq_m = []
    for t in range(0, 378, 18):
        dtheta_m.append(math.radians(t))
    for q in np.linspace( float(q2_min), float(q2_max), int(points)):
        dq_m.append(q)

    with open(template, 'r') as temp:
        contents = temp.read()

    for q in dq_m:
        for theta in dtheta_m:
            scan_point = []
            print(q, theta)
            for i in range(0, N):
                z = (2/N)**(1/2)*q*math.cos(theta+a_m[i])
                scan_point.append(z)
                print("writing....", z, q, theta)
            print(scan_point, len(scan_point))
            f = open(scan % (q, theta), 'w')
            smpl = contents.replace('%level%', level).replace('%qm%', str(q)).replace('%thetam%', str(theta)).replace('%z1%', str(scan_point[0]))\
            .replace('%z2%', str(scan_point[1])).replace('%z3%', str(scan_point[2])).replace('%z4%', str(scan_point[3])) \
            .replace('%z5%', str(scan_point[4]))
            f.write(smpl)
            f.close()

#Inputs

conf_file = "C:/Users/charv/Documents/CP_sampling/axis_test/CP_package/C5_base.xyz"
template='C:/Users/charv/Documents/CP_sampling/axis_test/CP_package/R5_g16temp.gjf'
scan= 'C:/Users/charv/Documents/CP_sampling/axis_test/CP_package/C5-%f_%f.gjf'
dft = 'b3lyp/cc-pVTZ'

sampleCP(baseconf=conf_file, template=template, scan=scan, level=dft, q2_min='0.0', q2_max='0.6', points='10')

