More details about the Potential Energy Surfaces and Cremer-Pople coordinate sampling can be found here: https://pubs.acs.org/doi/10.1021/acs.jpca.3c00095

# Cremer-Pople-sampling
    1. Pre-Optimization
The starting point of CP sampling is an optimized ring structure (mol.log). This step is to provide more reasonable internal coordinates for the extra-ring substituents and the optimization can be done at a low level of theory. From there we extract the atom Cartesian (mol_Cart.gjf) and internal (mol_int.gjf) coordinates. This is done manually from Gausview.

    2. Translation/Rotation: convert.py
The translation and rotation operations are done using only the Cartesian coordinate file (mol_Cart.gjf) as input. At this step it is crucial that the input (mol_Cart.gjf) structure matches the output (mol_base.xyz) structure so translation and rotation are applied to the whole molecule.  Required parameters:
    • N: number of ring atoms
    • N_all: total number of atoms in the molecule

    3. G16 template: make_template.py
This script makes a calculation template for each molecule. Cartesian coordinates for the ring atoms are  taken from the mol_base.xyz file and internal coordinates of the extra-ring substituents from the mol_int.gjf file. Inputs:
    • Cartesian coordinates (mol_coords.gjf)
    • Internal coordinates (mol_int.gjf)
    • N: number of ring atoms
    • N_all: total number of atoms in the molecule
    • template name

    4. CP sampling
At this step, we select the CP sampling parameters and desired level of theory for the sampling. The name of each sampling file is mol_[CPparas].gjf. Operations at this part are perfomed only on the ring atoms. Input files:
    • base coordinates: mol_base.xyz
    • template file
    • dft level
    • sampling range
Scan files naming format:
        ◦ 4-MR: mol_[q2 value].gjf
        ◦ 5-MR: mol_[q2 value]_[theta2 value].gjf
        ◦ 6-MR: mol_[Q value]_[Theta value]_[Phi value].gjf

Auxiliary scripts:
    • run.sh: submit a batch of G16 calculations inside the same directory 
    • control.sh: monitor calculation status finished, unfinished, failed
