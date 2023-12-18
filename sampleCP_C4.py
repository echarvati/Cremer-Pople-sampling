import numpy as np
import pandas as pd
from decimal import *


def sampleCP(baseconf="coords.xyz", template="g16temp.gjf", scan="outfile.gjf", level='level', q2='q2', points='points'):

    getcontext().prec = 6
    ###########files########

    N = 4
    m = 2

    baseconf = pd.read_csv(conf_file, skiprows=0, engine='python', header=None, delim_whitespace=True)
    df2 = baseconf.rename({0:'atom',1:'X', 2:'Y', 3:'Z'}, axis=1)

    x_conv = df2['X'].values
    y_conv = df2['Y'].values
    z_conv = df2['Z'].values
    z_C = []
    J=[]
    J_k=[]
    for c in range(0, N):
        # print(z_conv[c])
        z_C.append(z_conv[c])
        J.append(c)

    ####Calculate CP coordinates
    for j in J:
        k=(-1)**j
        J_k.append(int(k))
    q_m = (1/N)**(0.5)*(sum((z_C * np.asarray(J_k))))
    # print(q_m)

    print("Cremer-Pople coordinates at m=", m, "with q_m=", q_m)

    ###scan range
    dq_m=[]
    for q in np.linspace(-float(q2),float(q2), int(points)):
        dq_m.append(q)

    dz_scan = []

    for q in dq_m:
        for i in range(0, N):
                z_scan_p = N**(-1/2)*q*(-1)**i
                dz_scan.append(z_scan_p)
        for zsc in range(0, len(dz_scan), N):
            scan_point = dz_scan[zsc:zsc + N]
            print("writing....")
            print(scan_point, len(scan_point))
            with open(template, 'r') as temp:
                contents = temp.read()
                f = open(scan % (q), 'w')
                contents = contents.replace('%level%', level).replace('%q%', str(q)).replace('%q2%', str(q)).replace('%z1%', str(scan_point[0])) \
                .replace('%z2%', str(scan_point[1])).replace('%z3%', str(scan_point[2])).replace('%z4%',str(scan_point[3]))
                f.write(contents)
                f.close()

#files
conf_file = "C4_base.xyz"
template='C4_g16temp.gjf'
scan= 'C4_%f.gjf'
dft = 'b3lyp/cc-pVTZ'

sampleCP(baseconf=conf_file, template=template, scan=scan, level=dft, q2='0.4', points='500')

