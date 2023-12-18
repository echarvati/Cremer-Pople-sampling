import pandas as pd
import numpy as np
import math

def trans_rot(conf_file='conf', N_ring='N_ring', N_all='N_sys', out_conf='outfile' ):

    df = pd.read_csv(conf_file, skiprows=7, engine='python', header=None, delim_whitespace=True)
    system= df[:int(N_all)]
    atoms = system[0].to_list()
    ring= df[:int(N_ring)]
    ring_atoms = ring[0].to_list()


    X_init_all = system[1].to_numpy()
    X_init_ring = system[1][:int(N_ring)].to_numpy()
    Y_init_all = system[2].to_numpy()
    Y_init_ring = system[2][:int(N_ring)].to_numpy()
    Z_init_all = system[3].to_numpy()
    Z_init_ring = system[3][:int(N_ring)].to_numpy()

    Sys_Centroid = [X_init_all.mean(axis=0), Y_init_all.mean(axis=0), Z_init_all.mean(axis=0)]
    Ring_Centroid = [X_init_ring.mean(), Y_init_ring.mean(), Z_init_ring.mean()]
    OC_sys = np.linalg.norm(Sys_Centroid)
    OC_ring = np.linalg.norm(Ring_Centroid)

    print('System center:', Sys_Centroid)
    print('Distance system centroid from origin', OC_sys)
    print('Ring center:', Ring_Centroid)
    print('Distance ring centroid from origin', OC_ring)


    #####Translate the origin
    Trans_X_sys = X_init_all - Sys_Centroid[0]
    Trans_Y_sys = Y_init_all - Sys_Centroid[1]
    Trans_Z_sys = Z_init_all - Sys_Centroid[2]

    print('System after translation:', sum(Trans_X_sys), sum(Trans_Y_sys), sum(Trans_Z_sys))
    print('System centroid after translation:', Trans_X_sys.mean(), Trans_Y_sys.mean(), Trans_Z_sys.mean())

    Trans_Centroid_sys=[Trans_X_sys.mean(), Trans_Y_sys.mean(), Trans_Z_sys.mean()]
    Trans_OC_sys = np.linalg.norm(Trans_Centroid_sys)
    print('Distance system centroid from origin after translation:', Trans_OC_sys)

    ###Check

    J = np.arange(int(N_all), dtype=int)
    check1_z = sum(Trans_Z_sys)*np.cos(2*math.pi*J/int(N_all))
    check2_z = sum(Trans_Z_sys) * np.sin(2 * math.pi * J / int(N_all))

    check1 = np.isclose(check1_z, np.zeros(len(check1_z)), atol=1e-8)
    check2 = np.isclose(check2_z, np.zeros(len(check2_z)), atol=1e-8)

    if  np.all(check1)==True and np.all(check2)==True:
        print('Check passed')
    else:
        print('Check failed ')


    ######Rotate the z axis

    X_prime_sys = np.dot(Trans_X_sys,np.sin((2*math.pi*J)/int(N_all)))
    Y_prime_sys = np.dot(Trans_Y_sys,np.sin((2*math.pi*J)/int(N_all)))
    Z_prime_sys = np.dot(Trans_Z_sys,np.sin((2*math.pi*J)/int(N_all)))

    X_2prime_sys = np.dot(Trans_X_sys,np.cos((2*math.pi*J)/int(N_all)))
    Y_2prime_sys = np.dot(Trans_Y_sys,np.cos((2*math.pi*J)/int(N_all)))
    Z_2prime_sys = np.dot(Trans_Z_sys,np.cos((2*math.pi*J)/int(N_all)))

    R_prime_sys = np.array([X_prime_sys, Y_prime_sys,Z_prime_sys])
    R_2prime_sys = np.array([X_2prime_sys, Y_2prime_sys,Z_2prime_sys])
    unit_sys = (np.cross(R_prime_sys, R_2prime_sys)/np.linalg.norm(np.cross(R_prime_sys, R_2prime_sys)))

    Rot_Z = Trans_Z_sys*unit_sys[2]

    ####write coords

    template_trans=out_conf
    with open(template_trans, 'w') as temp_tr:
        for i in range(len(atoms)):
            r_i = (atoms[i], Trans_X_sys[i], Trans_Y_sys[i],Rot_Z[i])
            temp_tr.write(str(r_i).strip("(,)").replace("'", "").replace(",", ""))
            temp_tr.write("\n")

trans_rot(conf_file="OC6_chair_coords.gjf", N_ring="6", N_all="16", out_conf="OC6_base.xyz")
