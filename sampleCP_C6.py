import numpy as np
import math
import pandas as pd
from decimal import *

def scanCP(conf_file='base.xyz',template='g16temp.gjf', scan='outscan.gjf', level='level', Qmin='Qmin', Qmax='Qmax', points='points'):

    #####make initial structure
    getcontext().prec = 6
    N = 6
    m_1 = 2
    m_2 = 3
    ###Read file

    conf = pd.read_csv(conf_file, skiprows=7, engine='python', header=None, delim_whitespace=True)
    df2 = conf.rename({0:'atom',1:'X', 2:'Y', 3:'Z'}, axis=1)


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
    sin_jm1 = []
    cos_j=[]
    cos_jm1=[]
    r_sys=[]
    a_m1 = []
    cr_cos=[]

    for i in range(0, N):
        jm_1 = math.sin(2 * math.pi * i * m_1 / N)
        km_1 = math.cos(2 * math.pi * i * m_1 / N)
        cr_2 = math.cos(i*math.pi)
        a_1 = (2*math.pi*i*m_1)/N
        sin_jm1.append(jm_1)
        cos_jm1.append(km_1)
        cr_cos.append(float(cr_2))
        a_m1.append(a_1)
    print(cos_jm1)
    print(sin_jm1)
    print("crown", cr_cos, type(cr_cos))
    print(z_C, type(z_C))

    ####Calculate CP coordinates
    #pseudorotational subspace
    q_m1 = (2/N)**(0.5)*(sum((z_C* np.asarray(cos_jm1))**2) + sum((z_C* np.asarray(sin_jm1))**2))
    print(q_m1)
    theta_m1_degs = math.degrees(math.atan(sum(z_C* np.asarray(sin_jm1)))/(sum(z_C* np.asarray(cos_jm1))))
    print(theta_m1_degs)
    # if math.cos(theta_m1_degs)<0.0:
    #     theta_m_degs=math.degrees(theta_m1_degs+180)

    #crown form
    q_m2=(1/N)**(1/2)*sum(z_C*np.asarray(cr_cos))

    print("Initial Cremer-Pople coordinates at m=", m_1, "{q_2=", q_m1, ",theta_2=", theta_m1_degs, "} crown amplitude {q_3=", q_m2,"}")

    ###Spherical

    Q_ampl = ((q_m1)**2 + (q_m2)**2)**(1/2)
    print("Positve total puckering amplitude", Q_ampl)

    #########Displacement

    d_Theta = []
    d_Phi = []
    d_Q = []
    d_Phi = []
    cos_dPhi = []
    sin_dPhi = []

    for t in range(0, 190, 10):
        d_Theta.append(math.radians(t))

    for phi in range(0, 100, 10):
        d_Phi.append(math.radians(phi))

    for q in np.linspace(float(Qmin), float(Qmax), int(points)):
        d_Q.append(q)

    for ph in d_Phi:
        cos_dPhi.append(math.cos(ph))
        sin_dPhi.append(math.sin(ph))

    print(d_Q, d_Phi, d_Theta)

    z_scan = []

    with open(template, 'r') as temp:
        contents = temp.read()

    for q in d_Q:

        for p in d_Phi:
            for th in d_Theta:
                scan_point = []
                for i in range(0, N):
                    q2 = q * math.sin(p)
                    q3 = q * math.cos(p)
                    print(q3, q2, p, q, th, i)
                    z = (2/N)**(1/2)*q2*math.cos(th+a_m1[i]) + N**(-1/2)*q3*(-1)**i
                    scan_point.append(z)
                # print("writing....")
                print(scan_point, len(scan_point))
                print(scan_point[0], scan_point[1], scan_point[2])
                f = open(scan % (q, th, p), 'w')
                smth = contents.replace('%level%', level).replace('%Q%', str(q)).replace('%TH%', str(th)).replace('%PH%', str(p)).replace('%z1%', str(scan_point[0]))\
                .replace('%z2%', str(scan_point[1])).replace('%z3%', str(scan_point[2])).replace('%z4%', str(scan_point[3])) \
                .replace('%z5%', str(scan_point[4])).replace('%z6%', str(scan_point[5]))
                f.write(smth)
                f.close()
#files
base = "C6_base.xyz"
template='NNC4_g16temp.gjf'
scan='C6_%f_%f_%f.gjf'
dft = 'b3lyp/6-31+G(d,p)'

scanCP(conf_file=base, template=template, scan=scan, level=dft, Qmax='0.8', Qmin='0.0', points='16')