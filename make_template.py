import pandas as pd


def write_template(base_coords='base.xyz', base_int='opt_int.xyz', temp='temp.gjf', Nring='N', Nall='sys'):

    Z4 = ['%z1%', '%z2%', '%z3%', '%z4%']
    Z5 = ['%z1%', '%z2%', '%z3%', '%z4%', '%z5%']
    Z6 = ['%z1%', '%z2%', '%z3%', '%z4%', '%z5%', '%z6%']

    ring = []
    with open(base_coords, 'r'):
        df_coords = pd.read_csv(base_coords, header=None, delim_whitespace=True)
        # print(df_coords)
        atoms = df_coords[0].to_list()
        # print(len(atoms))
        for a in range(int(Nring)):
            ring.append(atoms[a])
            X = df_coords[1].to_list()
            Y = df_coords[2].to_list()
            Z = df_coords[3]
    print(ring)
    extra =[]
    with open(base_int, 'r'):
        #TODO fix skiprows thing!
        df_int = pd.read_csv(base_int,skiprows=7+int(Nring), header=None, delim_whitespace=True)
        extra_df = df_int.iloc[0:int(int(Nall)-int(Nring)),0:8].astype(object)
        extra_df[1] = extra_df[1].astype(int)
        extra_df[3] = extra_df[3].astype(int)
        extra_df[5] = extra_df[5].astype(int)
        extra_df[7] = extra_df[7].astype(int)
        print(extra_df)
        extra_string = extra_df.to_string(header=False, index=False)
        # print(extra_df[1])
        bond_vars = extra_df[2].to_list()
        angle_vars = extra_df[4].to_list()
        dh_vars = extra_df[6].to_list()
        vars_df = df_int.iloc[len(extra):len(df_int),0:2]
        variable = vars_df[0].to_list()
        value = vars_df[1].to_list()

    bond_variables=[]
    angle_variables=[]
    dh_variables=[]

    for v in range(len(variable)):
        for b in bond_vars:
            if b == variable[v]:
                bond = (b, value[v])
                bond_variables.append(bond)
        for a in angle_vars:
            if a == variable[v]:
                angle = (a, value[v])
                angle_variables.append(angle)
        for d in dh_vars:
            if d == variable[v]:
                dh = (d, value[v])
                dh_variables.append(dh)

    if int(Nring)==4:
        checkpoint = '%chk=scan_%q%.chk'
        title = 'Cremer-Pople coordinates at { q2=%qm% }'
    elif int(Nring)==5:
       checkpoint= '%chk=scan_%qm%-%thetam%.chk'
       title = 'Cremer-Pople coordinates at { q2=%qm% ,theta2=%thetam% }'
    elif int(Nring)==6:
       checkpoint= '%chk=scan_%Q%-%TH%-%PH%.chk'
       title = 'Cremer-Pople coordinates at { Q=%Q% ,Theta=%TH% ,Phi=%PH% }'

    with open(temp, 'w') as t:
        t.write('%NProcShared=4')
        t.write("\n")
        t.write('%mem=800MB')
        t.write("\n")
        t.write(checkpoint)
        t.write("\n")
        t.write('# opt=readallgic %level% pop=(regular,chelpg) EmpiricalDispersion=GD3')
        t.write("\n")
        t.write("\n")
        t.write(title)
        t.write("\n")
        t.write("\n")
        t.write("0 1")
        t.write("\n")
        for i in range(int(Nring)):
            if int(Nring)==4:
                r_i = (ring[i], X[i], Y[i], Z4[i])
                t.write(str(r_i).strip("(,)").replace("'", "").replace(",", ""))
                t.write("\n")
            if int(Nring)==5:
                r_i = ("", ring[i], X[i], Y[i], Z5[i])
                t.write(str(r_i).strip("(,)").replace("'", "").replace(",", ""))
                t.write("\n")
            if int(Nring)==6:
                r_i = (ring[i], X[i], Y[i], Z6[i])
                t.write(str(r_i).strip("(,)").replace("'", "").replace(",", ""))
                t.write("\n")
        t.write(extra_string)
        t.write("\n")
        t.write("\n")
        for bnd in bond_variables:
            t.write(str(bnd).strip("(,)").replace("'", "").replace(",", ""))
            t.write("\n")
        for ang in angle_variables:
            t.write(str(ang).strip("(,)").replace("'", "").replace(",", ""))
            t.write("\n")
        for h in dh_variables:
            t.write(str(h).strip("(,)").replace("'", "").replace(",", ""))
            t.write("\n")
        t.write("\n")
        for i in range(len(atoms)):
            j=i+1
            print(f"X{j}=X({j})", file=t)
        for m in range(len(atoms)):
            n=m+1
            print(f"Y{n}=Y({n})", file=t)
        for k in range(len(ring)):
            l=k+1
            print(f"Z{l}(freeze)=Z({l})", file=t)
        for q in range(len(ring), len(atoms)):
            o=q+1
            print(f"Z{o}=Z({o})", file=t)

write_template(base_coords="C4_base.xyz", base_int='C4_int.gjf', temp='C4_g16temp.gjf',Nring="4", Nall="12")
