# Loading AmberTools 24 environment
module load gcc/9.4.0  openmpi/4.1.2 amber/24
# Loading gaussian software
module load gaussian/G16RevC.02

InitialStructure=/path/to/drug/molecule/structure/pdb/or/mol2/format
WorkDir=/current/working/directory

# Drug molecule parameterization involves 
#(1) Optimizing drug molecule geometry using quantum computations
#(2) Generating electrostatic potential surrounding the drug's molecular surface using quantum computations
#(3) Running Restrained Electrostatic Potential (REsP) fitting approach to place partial charge on each atom of the drug molecule at the right intensity to mimic the electrostatic field outside of the drug's molecular surface that was computed in step 2
#(4) Generate drug molecule simulation parameters in format that are compatible with Gromacs simulation software


# 1. Geometry optimization performed using MP2 theory level with 6-31G* basis set
mkdir $WorkDir/1.Geometry_optimization
cd $WorkDir/1.Geometry_optimization

# Convert .mol2 file into .com file, which is input file for quantum computation using Gaussian software. Output file is named "YT0_C.com"
antechamber -i $InitialStructure -fi mol2 -o $WorkDir/1.Geometry_optimization/YT0_C.com -fo gcrt

# Modify the output file, which should looks like as follows:

	#################################################
--Link1--
%nproc=32
%mem=64GB
%chk=YT0_opt.chk
#P MP2/6-31G* opt

YT0 2nd geometry optimization

1   1
    C   18.6920000000       16.9560000000       22.0250000000     
    C   18.1540000000       18.0720000000       21.1860000000     
    C   17.0000000000       18.8530000000       21.4330000000     
    C   16.8920000000       19.9170000000       20.6080000000     
    C   15.7860000000       22.0400000000       20.0280000000     
    C   14.6770000000       24.0280000000       19.7950000000     
    C   15.4900000000       24.4060000000       18.7520000000     
    N   15.8900000000       20.8370000000       20.6920000000     
    C   15.2750000000       25.7550000000       18.1120000000     
    C   16.5020000000       23.4850000000       18.3780000000     
    C   18.4640000000       22.8040000000       16.8230000000     
    C   18.3090000000       23.5220000000       15.4470000000     
    C   17.7670000000       22.5680000000       14.3980000000     
    C   16.3210000000       22.5880000000       14.1950000000     
    C   19.6310000000       24.9960000000       13.9940000000     
    C   20.9910000000       25.6450000000       13.8430000000     
    N   14.8070000000       22.8650000000       20.4390000000     
    C   23.3270000000       25.2700000000       13.4390000000     
    C   22.0260000000       23.6070000000       14.5610000000     
    C   20.6670000000       22.9810000000       14.7780000000     
    C   17.2670000000       24.5610000000       16.1000000000     
    N   17.3740000000       23.6820000000       17.2730000000     
    N   15.1940000000       22.5790000000       13.9630000000     
    N   19.6680000000       24.0190000000       15.0910000000     
    N   22.0060000000       24.6410000000       13.5290000000     
    N   16.6420000000       22.3160000000       19.0230000000     
    S   18.1150000000       19.8280000000       19.3980000000     
    N   18.7490000000       18.4060000000       20.0360000000     
    H   16.2790000000       18.6220000000       22.2030000000     
    H   13.8960000000       24.7040000000       20.1080000000     
    H   15.1310000000       20.6090000000       21.3180000000     
    H   18.2280000000       21.7400000000       16.8180000000     
    H   19.4140000000       23.0240000000       17.3110000000     
    H   18.2400000000       22.8150000000       13.4480000000     
    H   18.0600000000       21.5550000000       14.6730000000     
    H   18.8870000000       25.7610000000       14.2150000000     
    H   19.3660000000       24.4900000000       13.0660000000     
    H   20.9510000000       26.3780000000       13.0380000000     
    H   21.2560000000       26.1460000000       14.7740000000     
    H   22.7260000000       22.8270000000       14.2610000000     
    H   22.3670000000       24.0480000000       15.4980000000     
    H   20.7250000000       22.2760000000       15.6070000000     
    H   20.3660000000       22.4510000000       13.8750000000     
    H   17.6410000000       25.5750000000       16.2390000000     
    H   16.2830000000       24.5300000000       15.6320000000     
    H   21.7830000000       24.2150000000       12.6410000000     
    H   23.3110000000       26.0410000000       12.6690000000     
    H   23.5810000000       25.7200000000       14.3990000000     
    H   24.0720000000       24.5160000000       13.1830000000     
    H   14.5720000000       26.3340000000       18.7110000000     
    H   16.2260000000       26.2850000000       18.0540000000     
    H   14.8720000000       25.6210000000       17.1080000000     
    H   18.0740000000       16.8390000000       22.9150000000     
    H   18.6790000000       16.0300000000       21.4500000000     
    H   19.7160000000       17.1850000000       22.3210000000     

	#################################################

# Running quantum computations for geometry optimization
g16 < YT0_C.com > YT0_C.log

# check "YT0_C.log" for detailed output from quantum computation executed above
# The output checkpoint file "YT0_opt.chk" can be converted to .pdb file for viewing the geometry of optimized drug molecule
newzmat -ichk -opdb $WorkDir/1.Geometry_optimization/YT0_opt.chk $WorkDir/1.Geometry_optimization/YT0_opt.pdb

# The log output "YT0_C.log" can also be converted to .mol2 file for viewing the geometry:
antechamber -i $WorkDir/1.Geometry_optimization/YT0_C.log -fi gout -o $WorkDir/1.Geometry_optimization/YT0_opt.mol2 -fo mol2


# 2. Electrostatic potential Charges calculation (ESP-charges)
#Using the geometrical optimized structure of drug from step 1, a Hartree-Fock (HF) calculation was performed to obtain its Electrostatic Potential.


mkdir $WorkDir/2.ESP-charges-calculation
cd $WorkDir/2.ESP-charges-calculation

# Convert .mol2 file into .com file, which is input file for quantum computation using Gaussian software. Output file is named "YT0_opt_charges.com"
antechamber -i $WorkDir/1.Geometry_optimization/YT0_opt.mol2 -fi mol2 -o $WorkDir/2.ESP-charges-calculation/YT0_opt_charges.com -fo gcrt

# Modify the output file "YT0_opt_charges.com", which should looks like as follows:

	#################################################
--Link1--
%nproc=32
%mem=64GB
%chk=YT0_opt_charges.chk
#P HF/6-31G*  SCF=Tight Pop=MK IOp(6/33=2)

YT0 ESP charges calculation

1   1
    C    5.8930000000       -3.6820000000        0.0330000000     
    C    4.8770000000       -2.5860000000       -0.0750000000     
    C    5.1810000000       -1.2240000000       -0.3260000000     
    C    4.0340000000       -0.4560000000       -0.3660000000     
    C    2.8590000000        1.6900000000       -0.6560000000     
    C    1.9390000000        3.7180000000       -1.0340000000     
    C    0.6490000000        3.2260000000       -0.8260000000     
    N    3.9830000000        0.9060000000       -0.6000000000     
    C   -0.5680000000        4.0830000000       -1.0210000000     
    C    0.5940000000        1.8660000000       -0.4770000000     
    C   -0.5640000000       -0.1000000000        0.4510000000     
    C   -1.9550000000        0.2190000000        1.0460000000     
    C   -2.1430000000       -0.0510000000        2.5500000000     
    C   -1.1050000000        0.5950000000        3.3530000000     
    C   -4.3550000000        0.1460000000        0.4310000000     
    C   -5.2060000000       -0.2220000000       -0.7720000000     
    N    3.0580000000        2.9850000000       -0.9530000000     
    C   -6.0610000000       -2.1400000000       -2.1150000000     
    C   -3.8470000000       -2.2980000000       -0.9800000000     
    C   -3.0630000000       -1.8340000000        0.2340000000     
    C   -1.6860000000        1.6850000000        0.6190000000     
    N   -0.6050000000        1.1900000000       -0.2580000000     
    N   -0.2580000000        1.1470000000        3.9640000000     
    N   -3.0150000000       -0.3790000000        0.2210000000     
    N   -5.2420000000       -1.7220000000       -0.9330000000     
    N    1.6860000000        1.0910000000       -0.4030000000     
    S    2.6690000000       -1.4890000000       -0.0860000000     
    N    3.5840000000       -2.8810000000        0.0730000000     
    H    6.1820000000       -0.8290000000       -0.4680000000     
    H    2.0810000000        4.7660000000       -1.2960000000     
    H    4.8550000000        1.3940000000       -0.7850000000     
    H    0.2250000000       -0.1410000000        1.2120000000     
    H   -0.4630000000       -0.9700000000       -0.2050000000     
    H   -3.1200000000        0.3150000000        2.8890000000     
    H   -2.1130000000       -1.1300000000        2.7470000000     
    H   -4.3280000000        1.2370000000        0.4990000000     
    H   -4.8380000000       -0.2270000000        1.3550000000     
    H   -6.2410000000        0.1170000000       -0.6720000000     
    H   -4.7690000000        0.1800000000       -1.6890000000     
    H   -3.9470000000       -3.3860000000       -1.0220000000     
    H   -3.3860000000       -1.9310000000       -1.8990000000     
    H   -2.0540000000       -2.2460000000        0.1560000000     
    H   -3.5180000000       -2.2510000000        1.1530000000     
    H   -2.4910000000        2.2130000000        0.0990000000     
    H   -1.3250000000        2.2980000000        1.4550000000     
    H   -5.6940000000       -2.1030000000       -0.0890000000     
    H   -7.0680000000       -1.7370000000       -2.0060000000     
    H   -5.5890000000       -1.7420000000       -3.0130000000     
    H   -6.0910000000       -3.2290000000       -2.1550000000     
    H   -0.2860000000        5.0290000000       -1.4900000000     
    H   -1.2910000000        3.5850000000       -1.6770000000     
    H   -1.0710000000        4.3220000000       -0.0780000000     
    H    6.4740000000       -3.7720000000       -0.8910000000     
    H    6.5940000000       -3.4890000000        0.8510000000     
    H    5.3860000000       -4.6300000000        0.2220000000     

	#################################################
	
# Running quantum computations for electrostatics potential surrounding the drug's molecular surface 
g16 < YT0_opt_charges.com > YT0_opt_charges.log


# The log output "YT0_C.log" can also be converted to mol2 for viewing the geometry:
antechamber -i $WorkDir/2.ESP-charges-calculation/YT0_opt_charges.log -fi gout -o $WorkDir/2.ESP-charges-calculation/YT0_opt_charges.mol2 -fo mol2


# 3.Restrained Fitting of the Partial Charges to the Electrostatic Potentials computed in step 2(RESP Calculation)
mkdir $WorkDir/3.RESP-calculation
cd $WorkDir/3.RESP-calculation

# The electrostatic potential must be in a format that the resp program can understand
#The following command (espgen) aims to extract the electrostatic potential from the Gaussian output
#The output file is named "esp.dat"
espgen -i $WorkDir/2.ESP-charges-calculation/YT0_opt_charges.log -o $WorkDir/3.RESP-calculation/esp.dat

# Then the following command generate estimated partial charge for each atom on drug molecule so that it mimics
#the electrostatics potential on drug molecular surface. The list of partial charges are stored on file "resp.chg"
resp -O -i $WorkDir/3.RESP-calculation/resp.in -o $WorkDir/3.RESP-calculation/resp.out -p $WorkDir/3.RESP-calculation/resp.pch -t $WorkDir/3.RESP-calculation/resp.chg -e $WorkDir/3.RESP-calculation/esp.dat

# Note that: file named "resp.in" needs to be manually generated
#Since atom with the same chemical properties should share similar partial charge.
#We can give resp fitting approach that information through resp.in file
#The description of this file can be further clarified here https://ambermd.org/tutorials/advanced/tutorial1/section1.php



# Now we have partial charges, we only needs to take care of bonded parameters

# The section below will detail how the bonded parameters were derived for drug molecule 
#MD simulation using Gromacs software
mkdir $WorkDir/4.Drug_Gromacs_parameters
cd $WorkDir/4.Drug_Gromacs_parameters


# List of partial charges for each atom is stored in file "resp.chg", generated in step 3
#and  bond parameters of drug molecule is defined with general amber force field 2 (gaff2)
# The gaffs has parameters for almost all organic molecules that are made up of C, N, O, H, S, P, F, Cl, Br, and I.
#Thus it can be used to generate parameters for the simulation of drug or small ligand.

antechamber -fi gout -i $WorkDir/2.ESP-charges-calculation/YT0_opt_charges.log -bk YT0 -fo mol2 -o $WorkDir/4.Drug_Gromacs_parameters/YT0_opt_charges-assigned_gaff2.mol2 -c rc -cf $WorkDir/3.RESP-calculation/resp.chg -nc 1 -at gaff2

# In case any parameters of drug molecule that are not already available on gaff2
#gaff2 will try to fill in these missing parameters by analogy to a similar parameter.

#We can check if all the parameters we need are available + how many parameters have been filled in through analogy 
#by executing the following command:
parmchk2 -i $WorkDir/4.Drug_Gromacs_parameters/YT0_opt_charges-assigned_gaff2.mol2 -f mol2 -o $WorkDir/4.Drug_Gromacs_parameters/YT0_missing_gaff2.frcmod -s gaff2

# If no analogy is found, either a default value that the antechamber tool think it is reasonable will be assigned
#or a place holder (with zeros everywhere) and the comment "ATTN: needs revision" will be inserted for the missing parameters
# In this case you will have to manually parameterise this yourself. 
# In the attached slide, we briefly demonstrated how the manual parameterization should be conducted

# Let's check the content of output file "YT0_missing_gaff2.frcmod"

	#################################################
Remark line goes here
MASS

BOND

ANGLE

DIHE
cc-cd-nu-ca   4    4.200       180.000           2.000      same as X -cd-nh-X , penalty score=  0.0
cc-cd-nu-hn   4    4.200       180.000           2.000      same as X -cd-nh-X , penalty score=  0.0
cc-cd-ss-nd   2    2.200       180.000           2.000      same as X -c2-ss-X , penalty score=232.0
ca-ca-nn-cy   4    4.200       180.000           2.000      same as X -ca-nh-X , penalty score=  0.0
nu-cd-ss-nd   2    2.200       180.000           2.000      same as X -c2-ss-X , penalty score=232.0
c1-c3-cy-cy   9    1.400         0.000           3.000      same as X -c3-c3-X , penalty score= 86.5
hc-c3-cy-cy   1    0.130         0.000           3.000      same as c3-c3-c3-hc, penalty score=173.5
cy-cy-n3-c6   1    0.730         0.000           3.000      same as c3-c3-n3-c3, penalty score=173.5
cy-cy-nn-ca   6    0.000         0.000           2.000      same as X -c3-nh-X , penalty score= 86.5
cy-cy-nn-cy   6    0.000         0.000           2.000      same as X -c3-nh-X , penalty score= 86.5
c3-cy-n3-c6   1    0.730         0.000           3.000      same as c3-c3-n3-c3, penalty score= 86.5
c6-c6-nx-c3   1    0.156         0.000           3.000      same as c3-c3-n4-c3, penalty score=  0.0
c6-c6-nx-c6   1    0.156         0.000           3.000      same as c3-c3-n4-c3, penalty score=  0.0
c6-c6-nx-hn   9    1.400         0.000           3.000      same as X -c3-n4-X , penalty score=  0.0
c6-c6-n3-cy   6    1.800         0.000           3.000      same as X -c3-n3-X , penalty score=  0.0
c6-c6-n3-c6   1    0.730         0.000           3.000      same as c3-c3-n3-c3, penalty score=  0.0
nb-ca-nu-cd   4    4.200       180.000           2.000      same as X -ca-nh-X , penalty score=  0.0
nb-ca-nu-hn   4    4.200       180.000           2.000      same as X -ca-nh-X , penalty score=  0.0
c1-c3-cy-n3   9    1.400         0.000           3.000      same as X -c3-c3-X , penalty score= 86.5
hc-c3-cy-n3   1    0.100         0.000           3.000      same as hc-c3-c3-n3, penalty score= 86.5
n3-c6-c6-nx   9    1.400         0.000           3.000      same as X -c3-c3-X , penalty score=  0.0
hx-c6-c6-n3   9    1.400         0.000           3.000      same as X -c3-c3-X , penalty score=  0.0
h1-c6-c6-nx   9    1.400         0.000           3.000      same as X -c3-c3-X , penalty score=  0.0
nb-ca-nn-cy   4    4.200       180.000           2.000      same as X -ca-nh-X , penalty score=  0.0
ss-cd-nu-ca   4    4.200       180.000           2.000      same as X -cd-nh-X , penalty score=  0.0
ss-cd-nu-hn   4    4.200       180.000           2.000      same as X -cd-nh-X , penalty score=  0.0
h1-cy-nn-ca   1    0.332         0.000           2.000      same as h1-c3-nh-ca, penalty score= 86.5
h1-cy-nn-cy   6    0.000         0.000           2.000      same as X -c3-nh-X , penalty score= 86.5
h1-c6-c6-hx   9    1.400         0.000           3.000      same as X -c3-c3-X , penalty score=  0.0
h1-c6-n3-cy   6    1.800         0.000           3.000      same as X -c3-n3-X , penalty score=  0.0
h1-c6-n3-c6   1    0.225         0.000           3.000      same as h1-c3-n3-c3, penalty score=  0.0
hx-c6-nx-c3   9    1.400         0.000           3.000      same as X -c3-n4-X , penalty score=  0.0
hx-c6-nx-c6   9    1.400         0.000           3.000      same as X -c3-n4-X , penalty score=  0.0
hx-c6-nx-hn   1    0.109         0.000           3.000      same as hx-c3-n4-hn, penalty score=  0.0
hx-c3-nx-c6   9    1.400         0.000           3.000      same as X -c3-n4-X , penalty score=  0.0
hx-c3-nx-hn   1    0.109         0.000           3.000      same as hx-c3-n4-hn, penalty score=  0.0

IMPROPER
c3-cc-cc-nd         1.1          180.0         2.0          Using the default value
cc-cd-cc-ha         1.1          180.0         2.0          Same as X -X -ca-ha, penalty score= 38.9 (use general term))
cc-nu-cd-ss         1.1          180.0         2.0          Using the default value
nb-nb-ca-nu        10.5          180.0         2.0          Same as X -n2-ca-n2, penalty score= 48.6 (use general term))
ca-h4-ca-nb         1.1          180.0         2.0          Same as X -X -ca-ha, penalty score= 44.3 (use general term))
ca-cd-nu-hn         1.1          180.0         2.0          Using the default value
ca-nb-ca-nn         1.1          180.0         2.0          Using the default value
ca-cy-nn-cy         1.1          180.0         2.0          Using the default value

NONBON

	#################################################

# In general, there are 36 missing dihedral terms and 8 missing improper term
# While the analogy can be found in all dihedral term
#5 improper terms are defined using default value 


# => If manual parameterization is too much to achieve in short time
#at least proper validation needs to be done.

# We will convert all the parameters that we have at this point into the format that Gromacs can understand
# Afterwards, the simulation of drug molecule in water box will be conducted to validate the paramters
mkdir $WorkDir/4.Drug_Gromacs_parameters
cd $WorkDir/4.Drug_Gromacs_parameters

tleapfile="source leaprc.gaff2\n
YT0 = loadmol2 $WorkDir/4.Drug_Gromacs_parameters/YT0_opt_charges-assigned_gaff2.mol2\n
loadamberparams $WorkDir/4.Drug_Gromacs_parameters/YT0_missing_gaff2.frcmod\n
saveoff YT0 $WorkDir/4.Drug_Gromacs_parameters/YT0.lib\n
saveamberparm YT0 $WorkDir/4.Drug_Gromacs_parameters/YT0.prmtop $WorkDir/4.Drug_Gromacs_parameters/YT0_opt_charges-assigned_gaff2.inpcrd\n
quit\n"

echo -e $tleapfile > YT0.tleap

# After executing the following command, the output files include "YT0.lib", "YT0_opt_charges-assigned_gaff2.inpcrd", and "YT0.prmtop"
tleap -f YT0.tleap


# Convert amber format to gromacs format using acpype tool, which is free and available here: https://github.com/alanwilter/acpype
acpype -p $WorkDir/4.Drug_Gromacs_parameters/YT0.prmtop -x $WorkDir/4.Drug_Gromacs_parameters/YT0_opt_charges-assigned_gaff2.inpcrd -b YT0

# The output directory named "YT0.amb2gmx" contains all the information we need to conduct an MD simulation of drug molecule
