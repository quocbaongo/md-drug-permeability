#!/bin/bash

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

set -e

WorkDir=/path/to/working/directory
# Solute coordinate and topology
Input_strct_coord=/path/to/drug/molecule/coordinates/YT0_GMX.gro
Input_strct_top=/path/to/drug/molecule/topology/in/Gromacs/format/YT0_GMX.itp

Lig_top_name=YT0

Solvent_coord=/path/to/Octanol/solvent/coordinate/Octanol_box_modified.pdb
Solvent_top=/path/to/octanol/topology/OCT_GMX.itp
Solvent_name=OCT

GMX=gmx_mpi
NCONFS=31


				################################################### Step 1 ###########################################################
				######################################## System establishment and equilibration ######################################

# System topology establishment
mkdir $WorkDir/System_topology
cd $WorkDir/System_topology

# Create system topology file
echo "[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               yes             0.5     0.8333333333

[ atomtypes ]
;name   bond_type     mass     charge   ptype   sigma         epsilon       Amb
 ho       ho          0.00000  0.00000   A     5.37925e-02   1.96648e-02 ; 0.30  0.0047
 hc       hc          0.00000  0.00000   A     2.60018e-01   8.70272e-02 ; 1.46  0.0208
 ha       ha          0.00000  0.00000   A     2.62548e-01   6.73624e-02 ; 1.47  0.0161
 h4       h4          0.00000  0.00000   A     2.53639e-01   6.73624e-02 ; 1.42  0.0161
 hn       hn          0.00000  0.00000   A     1.10650e-01   4.18400e-02 ; 0.62  0.0100
 h1       h1          0.00000  0.00000   A     2.42200e-01   8.70272e-02 ; 1.36  0.0208
 c3       c3          0.00000  0.00000   A     3.39771e-01   4.51035e-01 ; 1.91  0.1078
 cc       cc          0.00000  0.00000   A     3.31521e-01   4.13379e-01 ; 1.86  0.0988
 cd       cd          0.00000  0.00000   A     3.31521e-01   4.13379e-01 ; 1.86  0.0988
 ca       ca          0.00000  0.00000   A     3.31521e-01   4.13379e-01 ; 1.86  0.0988
 nu       nu          0.00000  0.00000   A     3.27904e-01   6.46428e-01 ; 1.84  0.1545
 cy       cy          0.00000  0.00000   A     3.39771e-01   4.51035e-01 ; 1.91  0.1078
 c1       c1          0.00000  0.00000   A     3.47896e-01   6.67766e-01 ; 1.95  0.1596
 c6       c6          0.00000  0.00000   A     3.39771e-01   4.51035e-01 ; 1.91  0.1078
 nb       nb          0.00000  0.00000   A     3.38417e-01   3.93714e-01 ; 1.90  0.0941
 nn       nn          0.00000  0.00000   A     3.18995e-01   8.99560e-01 ; 1.79  0.2150
 n1       n1          0.00000  0.00000   A     3.27352e-01   4.59403e-01 ; 1.84  0.1098
 n3       n3          0.00000  0.00000   A     3.36510e-01   3.58987e-01 ; 1.89  0.0858
 ss       ss          0.00000  0.00000   A     3.53241e-01   1.18156e+00 ; 1.98  0.2824
 nd       nd          0.00000  0.00000   A     3.38417e-01   3.93714e-01 ; 1.90  0.0941
 oh       oh          0.00000  0.00000   A     3.24287e-01   3.89112e-01 ; 1.82  0.0930

; Include solute parameters
#include \"$Input_strct_top\"

; Include solvent topology
#include \"$Solvent_top\"

[ system ]
; Name
solute in octanol solvent box

[ molecules ]
; Compound        #mols" >> $WorkDir/System_topology/system.top


# Establish simulation box
# Randonly placing drug molecule to pre-equilibration octanol solvent box
# This command should create a new file named: conf.gro in $WorkDir/System_topology
echo "System" | $GMX insert-molecules -f $Solvent_coord -ci $Input_strct_coord \
                                      -o $WorkDir/System_topology/conf.gro \
                                      -nmol 1 -try 10000 -scale 0.1

# Update system topology file 
echo 'import MDAnalysis as mda
import sys
from collections import Counter

def list_molecules(input_file):
    # Load the universe from .pdb or .gro
    universe = mda.Universe(input_file)
    
    # Get all unique resnames and count their occurrences (number of residues/molecules)
    resnames = [res.resname for res in universe.residues]
    molecule_counts = Counter(resnames)
    
    # Print each molecule name and its count
    for mol_name, count in molecule_counts.items():
        print(f"{mol_name}		  {count}")

# Usage: python script.py input.pdb (or .gro)
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_file.pdb or .gro>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    list_molecules(input_file)' >> $WorkDir/System_topology/extract_molecule.py

python3 $WorkDir/System_topology/extract_molecule.py $WorkDir/System_topology/conf.gro >> $WorkDir/System_topology/system.top


mkdir $WorkDir/Pre-equilibration
cd $WorkDir/Pre-equilibration

mkdir $WorkDir/Pre-equilibration/EM
cd $WorkDir/Pre-equilibration/EM

echo "; minim.mdp - used as input into grompp to generate em.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 100.0         ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet    ; Buffered neighbor searching
vdw-modifier    = Potential-shift-Verlet
coulombtype     = PME       ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
" > $WorkDir/Pre-equilibration/EM/em.mdp

# This command should create new files named: em.tpr and mdout.mdp in $WorkDir/Pre-equilibration/EM directory
$GMX grompp -f $WorkDir/Pre-equilibration/EM/em.mdp -p $WorkDir/System_topology/system.top \
	    -c $WorkDir/System_topology/conf.gro -o $WorkDir/Pre-equilibration/EM/em.tpr

# Energy minimization execution
cd $WorkDir/Pre-equilibration/EM

# This command should generate multiple new files including: confout.gro, ener.edr, md.log, and traj.trr in $WorkDir/Pre-equilibration/EM directory
srun $GMX mdrun -s em.tpr

# System equilibration with position restraints imposed on drug molecule
mkdir $WorkDir/Pre-equilibration/Equilibration
cd $WorkDir/Pre-equilibration/Equilibration

# Impose position restraint on drug molecule during equilibration
# LIG restraint
# This command should generate a new file named index.ndx in $WorkDir/Pre-equilibration/Equilibration directory
$GMX make_ndx -f $Input_strct_coord -o $WorkDir/Pre-equilibration/Equilibration/index.ndx <<EOF
0 & ! a H*
q
EOF

# This command should generate a new file named posre_LIG.itp in $WorkDir/System_topology directory
echo 'System_&_!H*' | $GMX genrestr -f $Input_strct_coord -n $WorkDir/Pre-equilibration/Equilibration/index.ndx \
				    -o $WorkDir/System_topology/posre_LIG.itp \
				    -fc 1000 1000 1000

# Edit existing system.top file in $WorkDir/System_topology
sed -i "\|#include \"$Input_strct_top\"|a\ \n; Ligand position restraints\n#ifdef POSRES\n#include \"$WorkDir/System_topology/posre_LIG.itp\"\n#endif" $WorkDir/System_topology/system.top


echo ";====================================================
; Equilibrium Simulations
;====================================================

; RUN CONTROL
;----------------------------------------------------
define       = -DPOSRES      ; use position restraints
integrator   = sd            ; stochastic leap-frog integrator
nsteps       = 500000        ; 2 * 500,000 fs = 1000 ps
dt           = 0.002         ; 2 fs
comm-mode    = Linear        ; remove center of mass translation
nstcomm      = 100           ; frequency for center of mass motion removal

; OUTPUT CONTROL
;----------------------------------------------------
nstxout                = 0         ; save coordinates to .trr every 100 ps
nstvout                = 0         ; save velocities to .trr every 100 ps
nstfout                = 0         ; save forces to .trr every 100 ps
nstxout-compressed     = 50000     ; xtc compressed trajectory output every 10 ps
compressed-x-precision = 1000      ; precision with which to write to the compressed trajectory file
nstlog                 = 10000     ; update log file every 20 ps
nstenergy              = 10000     ; save energies every 20 ps
nstcalcenergy          = 100       ; calculate energies every 100 steps (default=100)

; BONDS
;----------------------------------------------------
constraint_algorithm   = lincs      ; holonomic constraints
constraints            = h-bonds    ; constrain H-bonds
lincs-order            = 4
lincs-iter             = 1
lincs-warnangle        = 30         ; maximum angle that a bond can rotate before LINCS will complain (30 is default)
continuation           = no         ; formerly known as 'unconstrained-start' - useful for exact continuations and reruns

; NEIGHBOR SEARCHING
;----------------------------------------------------
cutoff-scheme         = Verlet ; group or Verlet
verlet-buffer-tolerance = 0.005
ns-type               = grid   ; search neighboring grid cells
nstlist               = 10     ; 20 fs (default is 10)
rlist                 = 1.0    ; short-range neighborlist cutoff (in nm)
pbc                   = xyz    ; 3D PBC
periodic-molecules    = no

; ELECTROSTATICS & EWALD
;----------------------------------------------------
coulombtype      = PME                       ; Particle Mesh Ewald for long-range electrostatics
coulomb-modifier = Potential-shift-Verlet
rcoulomb         = 1.0                       ; short-range electrostatic cutoff (in nm)
rcoulomb-switch  = 0
ewald-geometry   = 3d                        ; Ewald sum is performed in all three dimensions
pme-order        = 4                         ; interpolation order for PME (default is 4)
fourierspacing   = 0.12                      ; grid spacing for FFT
ewald-rtol       = 1e-5                      ; relative strength of the Ewald-shifted direct potential at rcoulomb

; VAN DER WAALS
;----------------------------------------------------
vdw-type          = Cut-off      ; potential switched off at rvdw-switch to reach zero at rvdw
vdw-modifier      = Potential-shift-Verlet
rvdw              = 1.0          ; van der Waals cutoff (in nm)
DispCorr          = EnerPres     ; apply analytical long range dispersion corrections for Energy and Pressure

; TEMPERATURE COUPLING (Langevin)
;----------------------------------------------------
tcoupl            = no
tc-grps           = System
tau-t             = 1.0
ref-t             = 300
gen-vel           = yes           ; Velocity generation (if gen-vel is 'yes', continuation should be 'no')
gen-temp          = 300
gen-seed          = -1

; PRESSURE COUPLING
;----------------------------------------------------
pcoupl           = C-rescale
pcoupltype       = isotropic
nstpcouple       = 10
tau_p            = 1.0                  ; time constant (ps)
ref_p            = 1.0
compressibility  = 4.5e-5
refcoord-scaling = com	
" >> $WorkDir/Pre-equilibration/Equilibration/equil.mdp

# This command should create new files named: npt.tpr and mdout.mdp in $WorkDir/Pre-equilibration/Equilibration directory
$GMX grompp -f $WorkDir/Pre-equilibration/Equilibration/equil.mdp -p $WorkDir/System_topology/system.top \
	    -c $WorkDir/Pre-equilibration/EM/confout.gro -r $WorkDir/Pre-equilibration/EM/confout.gro \
	    -o $WorkDir/Pre-equilibration/Equilibration/npt.tpr


cd $WorkDir/Pre-equilibration/Equilibration

# This command should generate multiple new files including: confout.gro, md.log, traj_comp.xtc, ener.edr, and state.cpt in $WorkDir/Pre-equilibration/Equilibration
srun $GMX mdrun -s npt.tpr 


# System equilibration without position restraints
mkdir $WorkDir/Pre-equilibration/Equilibration_no_restraints
cd $WorkDir/Pre-equilibration/Equilibration_no_restraints


echo ";====================================================
; Equilibrium Simulations
;====================================================

; RUN CONTROL
;----------------------------------------------------
integrator   = sd            ; stochastic leap-frog integrator
nsteps       = 50000000      ; 1 * 50,000,000 fs = 100000 ps
dt           = 0.002         ; 2 fs
comm-mode    = Linear        ; remove center of mass translation
nstcomm      = 100           ; frequency for center of mass motion removal

; OUTPUT CONTROL
;----------------------------------------------------
nstxout                = 0         ; save coordinates to .trr every 100 ps
nstvout                = 0         ; save velocities to .trr every 100 ps
nstfout                = 0         ; save forces to .trr every 100 ps
nstxout-compressed     = 50000     ; xtc compressed trajectory output every 10 ps
compressed-x-precision = 1000      ; precision with which to write to the compressed trajectory file
nstlog                 = 10000     ; update log file every 20 ps
nstenergy              = 10000     ; save energies every 20 ps
nstcalcenergy          = 100       ; calculate energies every 100 steps (default=100)

; BONDS
;----------------------------------------------------
constraint_algorithm   = lincs      ; holonomic constraints
constraints            = h-bonds    ; constrain H-bonds
lincs-order            = 4
lincs-iter             = 1
lincs-warnangle        = 30         ; maximum angle that a bond can rotate before LINCS will complain (30 is default)
continuation           = yes        ; formerly known as 'unconstrained-start' - useful for exact continuations and reruns

; NEIGHBOR SEARCHING
;----------------------------------------------------
cutoff-scheme         = Verlet ; group or Verlet
verlet-buffer-tolerance = 0.005
ns-type               = grid   ; search neighboring grid cells
nstlist               = 10     ; 20 fs (default is 10)
rlist                 = 1.0    ; short-range neighborlist cutoff (in nm)
pbc                   = xyz    ; 3D PBC
periodic-molecules    = no

; ELECTROSTATICS & EWALD
;----------------------------------------------------
coulombtype      = PME                       ; Particle Mesh Ewald for long-range electrostatics
coulomb-modifier = Potential-shift-Verlet
rcoulomb         = 1.0                       ; short-range electrostatic cutoff (in nm)
rcoulomb-switch  = 0
ewald-geometry   = 3d                        ; Ewald sum is performed in all three dimensions
pme-order        = 4                         ; interpolation order for PME (default is 4)
fourierspacing   = 0.12                      ; grid spacing for FFT
ewald-rtol       = 1e-5                      ; relative strength of the Ewald-shifted direct potential at rcoulomb

; VAN DER WAALS
;----------------------------------------------------
vdw-type          = Cut-off      ; potential switched off at rvdw-switch to reach zero at rvdw
vdw-modifier      = Potential-shift-Verlet
rvdw              = 1.0          ; van der Waals cutoff (in nm)
DispCorr          = EnerPres     ; apply analytical long range dispersion corrections for Energy and Pressure

; TEMPERATURE COUPLING (Langevin)
;----------------------------------------------------
tcoupl            = no
tc-grps           = System
tau-t             = 1.0
ref-t             = 300
gen-vel           = no           ; Velocity generation (if gen-vel is 'yes', continuation should be 'no')
gen-temp          = 300
gen-seed          = -1

; PRESSURE COUPLING
;----------------------------------------------------
pcoupl           = C-rescale
pcoupltype       = isotropic
nstpcouple       = 10
tau_p            = 1.0                  ; time constant (ps)
ref_p            = 1.0
compressibility  = 4.5e-5
refcoord-scaling = no	
" >> $WorkDir/Pre-equilibration/Equilibration_no_restraints/equil.mdp

# This command should create new files named: md.tpr and mdout.mdp in $WorkDir/Pre-equilibration/Equilibration_no_restraints
$GMX grompp -f $WorkDir/Pre-equilibration/Equilibration_no_restraints/equil.mdp \
	    -p $WorkDir/System_topology/system.top \
	    -c $WorkDir/Pre-equilibration/Equilibration/confout.gro \
	    -t $WorkDir/Pre-equilibration/Equilibration/state.cpt \
	    -o $WorkDir/Pre-equilibration/Equilibration_no_restraints/md.tpr

cd $WorkDir/Pre-equilibration/Equilibration_no_restraints

# This command should generate multiple new files including: confout.gro, md.log, traj_comp.xtc, ener.edr, state.cpt, and state_prev.cpt in $WorkDir/Pre-equilibration/Equilibration_no_restraints
srun $GMX mdrun -s md.tpr -e ener.edr -g md.log -cpi state.cpt -x traj_comp.xtc

# TThis command should create new file named: confout100ns_whole_molecule.gro in $WorkDir/Pre-equilibration/Equilibration_no_restraints
echo "System" | $GMX trjconv -f $WorkDir/Pre-equilibration/Equilibration_no_restraints/traj_comp.xtc \
			     -s $WorkDir/Pre-equilibration/EM/em.tpr \
			     -o $WorkDir/Pre-equilibration/Equilibration_no_restraints/confout100ns_whole_molecule.gro \
			     -ur compact -pbc mol -dump 100000


# Fix the PBC
# The following command aims to bringing all the molecules in the simulated system to the center of the simulation box only for visualization purpose

mkdir $WorkDir/Pre-equilibration/Equilibration_no_restraints/FixPBC
cd $WorkDir/Pre-equilibration/Equilibration_no_restraints/FixPBC

$GMX make_ndx -f $Input_strct_coord -o $WorkDir/Pre-equilibration/Equilibration_no_restraints/FixPBC/index.ndx <<EOF
q
EOF

echo "$Lig_top_name" "$Lig_top_name" | $GMX trjconv -s $WorkDir/Pre-equilibration/EM/em.tpr \
				-f $WorkDir/Pre-equilibration/Equilibration_no_restraints/traj_comp.xtc \
				-o $WorkDir/Pre-equilibration/Equilibration_no_restraints/FixPBC/md.xtc \
				-center -pbc mol -ur compact

echo "$Lig_top_name" "$Lig_top_name" | $GMX trjconv -f $WorkDir/Pre-equilibration/Equilibration_no_restraints/FixPBC/md.xtc \
				-s $WorkDir/Pre-equilibration/EM/em.tpr \
				-o $WorkDir/Pre-equilibration/Equilibration_no_restraints/FixPBC/md.fit.xtc \
				-n $WorkDir/Pre-equilibration/Equilibration_no_restraints/FixPBC/index.ndx \
				-fit rot+trans



cd $WorkDir/Pre-equilibration/Equilibration_no_restraints/FixPBC
echo "$Lig_top_name" | $GMX mindist -s $WorkDir/Pre-equilibration/EM/em.tpr -f $WorkDir/Pre-equilibration/Equilibration_no_restraints/FixPBC/md.fit.xtc \
			  -pi -od $WorkDir/Pre-equilibration/Equilibration_no_restraints/FixPBC/mindist.xvg \
			  -n $WorkDir/Pre-equilibration/Equilibration_no_restraints/FixPBC/index.ndx 



				#################################################### Step 2 ##############################################################
				######################################## Gradually turning off electrostatics + vdw ######################################


InputStructure=$WorkDir/Pre-equilibration/Equilibration_no_restraints/confout100ns_whole_molecule.gro
InputState=$WorkDir/Pre-equilibration/Equilibration_no_restraints/state.cpt

# Gradually turning off electrostatics
mkdir $WorkDir/Starting_conformations_generation
cd $WorkDir/Starting_conformations_generation

mkdir $WorkDir/Starting_conformations_generation/Electrostatics_decoupling
cd $WorkDir/Starting_conformations_generation/Electrostatics_decoupling

echo ";====================================================
; Equilibrium Simulations
;====================================================

; RUN CONTROL
;----------------------------------------------------
integrator   = sd            ; stochastic leap-frog integrator
nsteps       = 50000000      ; 2 * 50,000,000 fs = 100000 ps
dt           = 0.002         ; 2 fs
comm-mode    = Linear        ; remove center of mass translation
nstcomm      = 100           ; frequency for center of mass motion removal

; OUTPUT CONTROL
;----------------------------------------------------
nstxout                = 10000     ; save coordinates to .trr every 100 ps
nstvout                = 10000     ; save velocities to .trr every 100 ps
nstfout                = 10000     ; save forces to .trr every 100 ps
nstxout-compressed     = 50000     ; xtc compressed trajectory output every 10 ps
compressed-x-precision = 1000      ; precision with which to write to the compressed trajectory file
nstlog                 = 10000     ; update log file every 20 ps
nstenergy              = 10000     ; save energies every 20 ps
nstcalcenergy          = 100       ; calculate energies every 100 steps (default=100)

; BONDS
;----------------------------------------------------
constraint_algorithm   = lincs      ; holonomic constraints
constraints            = h-bonds    ; constrain H-bonds
lincs-order            = 4
lincs-iter             = 1
lincs-warnangle        = 30         ; maximum angle that a bond can rotate before LINCS will complain (30 is default)
continuation           = yes        ; formerly known as 'unconstrained-start' - useful for exact continuations and reruns

; NEIGHBOR SEARCHING
;----------------------------------------------------
cutoff-scheme         = Verlet ; group or Verlet
verlet-buffer-tolerance = 0.005
ns-type               = grid   ; search neighboring grid cells
nstlist               = 10     ; 20 fs (default is 10)
rlist                 = 1.0    ; short-range neighborlist cutoff (in nm)
pbc                   = xyz    ; 3D PBC
periodic-molecules    = no

; ELECTROSTATICS & EWALD
;----------------------------------------------------
coulombtype      = PME                       ; Particle Mesh Ewald for long-range electrostatics
coulomb-modifier = Potential-shift-Verlet
rcoulomb         = 1.0                       ; short-range electrostatic cutoff (in nm)
rcoulomb-switch  = 0
ewald-geometry   = 3d                        ; Ewald sum is performed in all three dimensions
pme-order        = 4                         ; interpolation order for PME (default is 4)
fourierspacing   = 0.12                      ; grid spacing for FFT
ewald-rtol       = 1e-5                      ; relative strength of the Ewald-shifted direct potential at rcoulomb

; VAN DER WAALS
;----------------------------------------------------
vdw-type          = Cut-off      ; potential switched off at rvdw-switch to reach zero at rvdw
vdw-modifier      = Potential-shift-Verlet
rvdw              = 1.0          ; van der Waals cutoff (in nm)
DispCorr          = EnerPres     ; apply analytical long range dispersion corrections for Energy and Pressure

; TEMPERATURE COUPLING (Langevin)
;----------------------------------------------------
tcoupl            = no
tc-grps           = System
tau-t             = 1.0
ref-t             = 300
gen-vel           = no           ; Velocity generation (if gen-vel is 'yes', continuation should be 'no')
gen-temp          = 300
gen-seed          = -1

; PRESSURE COUPLING
;----------------------------------------------------
pcoupl           = C-rescale
pcoupltype       = isotropic
nstpcouple       = 10
tau_p            = 1.0                  ; time constant (ps)
ref_p            = 1.0
compressibility  = 4.5e-5
refcoord-scaling = no

; FREE ENERGY
;----------------------------------------------------
free-energy	= yes
init-lambda	= 0       
delta-lambda    = 2e-08
couple-moltype  = $Lig_top_name
couple-lambda0  = vdw-q
couple-lambda1  = vdw
couple-intramol = yes
nstdhdl         = 100
sc-alpha        = 0.0" >> $WorkDir/Starting_conformations_generation/Electrostatics_decoupling/md.mdp

# This command should create new files named: tpr.tpr and mdout.mdp in $WorkDir/Starting_conformations_generation/Electrostatics_decoupling directory
$GMX grompp -f $WorkDir/Starting_conformations_generation/Electrostatics_decoupling/md.mdp \
	    -c $InputStructure \
	    -t $InputState \
	    -p $WorkDir/System_topology/system.top \
	    -o $WorkDir/Starting_conformations_generation/Electrostatics_decoupling/tpr.tpr

# This command should generate multiple new files including: confout.gro, md.log, traj_comp.xtc, traj.trr, dhdl.xvg ener.edr, state_prev.cpt, and state.cpt in $WorkDir/Starting_conformations_generation/Electrostatics_decoupling directory
srun $GMX mdrun -s tpr.tpr

# Gradually turning off vdw
mkdir $WorkDir/Starting_conformations_generation/vdw_decoupling
cd $WorkDir/Starting_conformations_generation/vdw_decoupling


echo ";====================================================
; Equilibrium Simulations
;====================================================

; RUN CONTROL
;----------------------------------------------------
integrator   = sd            ; stochastic leap-frog integrator
nsteps       = 50000000      ; 2 * 50,000,000 fs = 100000 ps
dt           = 0.002         ; 2 fs
comm-mode    = Linear        ; remove center of mass translation
nstcomm      = 100           ; frequency for center of mass motion removal

; OUTPUT CONTROL
;----------------------------------------------------
nstxout                = 10000     ; save coordinates to .trr every 100 ps
nstvout                = 10000     ; save velocities to .trr every 100 ps
nstfout                = 10000     ; save forces to .trr every 100 ps
nstxout-compressed     = 50000     ; xtc compressed trajectory output every 10 ps
compressed-x-precision = 1000      ; precision with which to write to the compressed trajectory file
nstlog                 = 10000     ; update log file every 20 ps
nstenergy              = 10000     ; save energies every 20 ps
nstcalcenergy          = 100       ; calculate energies every 100 steps (default=100)

; BONDS
;----------------------------------------------------
constraint_algorithm   = lincs      ; holonomic constraints
constraints            = h-bonds    ; constrain H-bonds
lincs-order            = 4
lincs-iter             = 1
lincs-warnangle        = 30         ; maximum angle that a bond can rotate before LINCS will complain (30 is default)
continuation           = yes        ; formerly known as 'unconstrained-start' - useful for exact continuations and reruns

; NEIGHBOR SEARCHING
;----------------------------------------------------
cutoff-scheme         = Verlet ; group or Verlet
verlet-buffer-tolerance = 0.005
ns-type               = grid   ; search neighboring grid cells
nstlist               = 10     ; 20 fs (default is 10)
rlist                 = 1.0    ; short-range neighborlist cutoff (in nm)
pbc                   = xyz    ; 3D PBC
periodic-molecules    = no

; ELECTROSTATICS & EWALD
;----------------------------------------------------
coulombtype      = PME                       ; Particle Mesh Ewald for long-range electrostatics
coulomb-modifier = Potential-shift-Verlet
rcoulomb         = 1.0                       ; short-range electrostatic cutoff (in nm)
rcoulomb-switch  = 0
ewald-geometry   = 3d                        ; Ewald sum is performed in all three dimensions
pme-order        = 4                         ; interpolation order for PME (default is 4)
fourierspacing   = 0.12                      ; grid spacing for FFT
ewald-rtol       = 1e-5                      ; relative strength of the Ewald-shifted direct potential at rcoulomb

; VAN DER WAALS
;----------------------------------------------------
vdw-type          = Cut-off      ; potential switched off at rvdw-switch to reach zero at rvdw
vdw-modifier      = Potential-shift-Verlet
rvdw              = 1.0          ; van der Waals cutoff (in nm)
DispCorr          = EnerPres     ; apply analytical long range dispersion corrections for Energy and Pressure

; TEMPERATURE COUPLING (Langevin)
;----------------------------------------------------
tcoupl            = no
tc-grps           = System
tau-t             = 1.0
ref-t             = 300
gen-vel           = no           ; Velocity generation (if gen-vel is 'yes', continuation should be 'no')
gen-temp          = 300
gen-seed          = -1

; PRESSURE COUPLING
;----------------------------------------------------
pcoupl           = C-rescale
pcoupltype       = isotropic
nstpcouple       = 10
tau_p            = 1.0                  ; time constant (ps)
ref_p            = 1.0
compressibility  = 4.5e-5
refcoord-scaling = no

; FREE ENERGY
;----------------------------------------------------
free-energy	= yes
init-lambda	= 0       
delta-lambda    = 2e-08
couple-moltype  = $Lig_top_name
couple-lambda0  = vdw
couple-lambda1  = none
couple-intramol = yes
nstdhdl         = 100
sc-alpha        = 0.5
sc-power        = 1
sc-sigma        = 0.3" >> $WorkDir/Starting_conformations_generation/vdw_decoupling/md.mdp

# This command should create new files named: tpr.tpr and mdout.mdp in $WorkDir/Starting_conformations_generation/vdw_decoupling directory
$GMX grompp -f $WorkDir/Starting_conformations_generation/vdw_decoupling/md.mdp \
	    -c $WorkDir/Starting_conformations_generation/Electrostatics_decoupling/confout.gro \
	    -t $WorkDir/Starting_conformations_generation/Electrostatics_decoupling/state.cpt \
	    -p $WorkDir/System_topology/system.top \
	    -o $WorkDir/Starting_conformations_generation/vdw_decoupling/tpr.tpr


# This command should generate multiple new files including: confout.gro, md.log, traj_comp.xtc, traj.trr, dhdl.xvg ener.edr, state_prev.cpt, and state.cpt in $WorkDir/Starting_conformations_generation/vdw_decoupling directory
srun $GMX mdrun -s tpr.tpr

# Extracting corresponding conformations
mkdir $WorkDir/Starting_conformations_generation/Extracted_conformations
cd $WorkDir/Starting_conformations_generation/Extracted_conformations

delta_lambda='2*10^-8'

# Decoupling electrostatics conformations
vdw_lambdas=(0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00)
coul_lambdas=(0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00)

for i in "${!coul_lambdas[@]}"
do
	coul_lambda=${coul_lambdas[$i]}
	vdw_lambda=${vdw_lambdas[$i]}
	time_step=$(echo "scale=10; $coul_lambda / ($delta_lambda) * 0.002" | bc)
	
	
	# This command will generate multiple new files with name starting with conf_ and ending with .gro in $WorkDir/Starting_conformations_generation/Extracted_conformations directory
	echo "System" | $GMX trjconv -f $WorkDir/Starting_conformations_generation/Electrostatics_decoupling/traj_comp.xtc \
				     -s $WorkDir/Starting_conformations_generation/Electrostatics_decoupling/tpr.tpr \
				     -dump $time_step \
				     -o $WorkDir/Starting_conformations_generation/Extracted_conformations/conf_${coul_lambda}_${vdw_lambda}.gro

done

# Decoupling vdw conformations
vdw_lambdas=(0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00)
coul_lambdas=(1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00)


for i in "${!vdw_lambdas[@]}"
do
	coul_lambda=${coul_lambdas[$i]}
	vdw_lambda=${vdw_lambdas[$i]}
	
	time_step=$(echo "scale=10; $vdw_lambda / ($delta_lambda) * 0.002" | bc)
	
	# This command will generate multiple new files with name starting with conf_ and ending with .gro in $WorkDir/Starting_conformations_generation/Extracted_conformations directory
	echo "System" | $GMX trjconv -f $WorkDir/Starting_conformations_generation/vdw_decoupling/traj_comp.xtc \
				     -s $WorkDir/Starting_conformations_generation/vdw_decoupling/tpr.tpr \
				     -dump $time_step \
				     -o $WorkDir/Starting_conformations_generation/Extracted_conformations/conf_${coul_lambda}_${vdw_lambda}.gro

done



				################################################### Step 3 ###########################################################
				########################## Thorough sampling for solvation free energy compuations ###################################


# Prepare simulation parameters fles (.mdp)
mkdir $WorkDir/Solvation_free_energy
mkdir $WorkDir/Solvation_free_energy/EM
cd $WorkDir/Solvation_free_energy/EM

# EM
echo "integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 100.0         ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet    ; Buffered neighbor searching
vdw-modifier    = Potential-shift-Verlet
coulombtype     = PME       ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions

free-energy              = yes
init-lambda-state        = MYLAMBDA
calc-lambda-neighbors    = -1
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30
vdw_lambdas              = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00
coul_lambdas             = 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00
couple-moltype           = $Lig_top_name
couple-lambda0           = vdw-q    ; only van der Waals interactions
couple-lambda1           = none     ; turn off everything, in this case only vdW
couple-intramol          = yes
separate-dhdl-file       = yes
nstdhdl                  = 100
sc-alpha                 = 0.5
sc-coul                  = no
sc-power                 = 1
sc-sigma                 = 0.3
" >> $WorkDir/Solvation_free_energy/EM/em.mdp


vdw_lambdas=(0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00)
coul_lambdas=(0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00)


# The following commands will create multiple sub-directories with name starting with "Lambda" in $WorkDir/Solvation_free_energy/EM directory
# Each of those sub-directories will contain two new files named em.tpr and mdout.mdp
for i in "${!coul_lambdas[@]}"
do
	coul_lambda=${coul_lambdas[$i]}
	vdw_lambda=${vdw_lambdas[$i]}

	mkdir $WorkDir/Solvation_free_energy/EM/Lambda$i
	cd $WorkDir/Solvation_free_energy/EM/Lambda$i
	
	sed 's/MYLAMBDA/'$i'/g' $WorkDir/Solvation_free_energy/EM/em.mdp > $WorkDir/Solvation_free_energy/EM/Lambda$i/em.mdp
	
	$GMX grompp -f $WorkDir/Solvation_free_energy/EM/Lambda$i/em.mdp \
		    -c $WorkDir/Starting_conformations_generation/Extracted_conformations/conf_${coul_lambda}_${vdw_lambda}.gro \
		    -p $WorkDir/System_topology/system.top \
		    -o $WorkDir/Solvation_free_energy/EM/Lambda$i/em.tpr
	
done

# Energy minimization on parallel
cd $WorkDir/Solvation_free_energy/EM

# This command will generate multiple new files including confout.gro, ener.edr, md.log and traj.trr in each of the "Lambda" sub-directories in $WorkDir/Solvation_free_energy/EM directory
srun $GMX mdrun -s em.tpr -multidir $WorkDir/Solvation_free_energy/EM/Lambda*


# Equilibration
mkdir $WorkDir/Solvation_free_energy/Equilibration
cd $WorkDir/Solvation_free_energy/Equilibration

echo ";====================================================
; Equilibrium Simulations
;====================================================

; RUN CONTROL
;----------------------------------------------------
define       = -DPOSRES      ; use position restraints
integrator   = sd            ; stochastic leap-frog integrator
nsteps       = 500000        ; 2 * 500,000 fs = 1000 ps
dt           = 0.002         ; 2 fs
comm-mode    = Linear        ; remove center of mass translation
nstcomm      = 100           ; frequency for center of mass motion removal

; OUTPUT CONTROL
;----------------------------------------------------
nstxout                = 0         ; save coordinates to .trr every 100 ps
nstvout                = 0         ; save velocities to .trr every 100 ps
nstfout                = 0         ; save forces to .trr every 100 ps
nstxout-compressed     = 50000     ; xtc compressed trajectory output every 10 ps
compressed-x-precision = 1000      ; precision with which to write to the compressed trajectory file
nstlog                 = 10000     ; update log file every 20 ps
nstenergy              = 10000     ; save energies every 20 ps
nstcalcenergy          = 100       ; calculate energies every 100 steps (default=100)

; BONDS
;----------------------------------------------------
constraint_algorithm   = lincs      ; holonomic constraints
constraints            = h-bonds    ; constrain H-bonds
lincs-order            = 4
lincs-iter             = 1
lincs-warnangle        = 30         ; maximum angle that a bond can rotate before LINCS will complain (30 is default)
continuation           = no         ; formerly known as 'unconstrained-start' - useful for exact continuations and reruns

; NEIGHBOR SEARCHING
;----------------------------------------------------
cutoff-scheme         = Verlet ; group or Verlet
verlet-buffer-tolerance = 0.005
ns-type               = grid   ; search neighboring grid cells
nstlist               = 10     ; 20 fs (default is 10)
rlist                 = 1.0    ; short-range neighborlist cutoff (in nm)
pbc                   = xyz    ; 3D PBC
periodic-molecules    = no

; ELECTROSTATICS & EWALD
;----------------------------------------------------
coulombtype      = PME                       ; Particle Mesh Ewald for long-range electrostatics
coulomb-modifier = Potential-shift-Verlet
rcoulomb         = 1.0                       ; short-range electrostatic cutoff (in nm)
rcoulomb-switch  = 0
ewald-geometry   = 3d                        ; Ewald sum is performed in all three dimensions
pme-order        = 4                         ; interpolation order for PME (default is 4)
fourierspacing   = 0.12                      ; grid spacing for FFT
ewald-rtol       = 1e-5                      ; relative strength of the Ewald-shifted direct potential at rcoulomb

; VAN DER WAALS
;----------------------------------------------------
vdw-type          = Cut-off      ; potential switched off at rvdw-switch to reach zero at rvdw
vdw-modifier      = Potential-shift-Verlet
rvdw              = 1.0          ; van der Waals cutoff (in nm)
DispCorr          = EnerPres     ; apply analytical long range dispersion corrections for Energy and Pressure

; TEMPERATURE COUPLING (Langevin)
;----------------------------------------------------
tcoupl            = no
tc-grps           = System
tau-t             = 1.0
ref-t             = 300
gen-vel           = yes           ; Velocity generation (if gen-vel is 'yes', continuation should be 'no')
gen-temp          = 300
gen-seed          = -1

; PRESSURE COUPLING
;----------------------------------------------------
pcoupl           = C-rescale
pcoupltype       = isotropic
nstpcouple       = 10
tau_p            = 1.0                  ; time constant (ps)
ref_p            = 1.0
compressibility  = 4.5e-5
refcoord-scaling = com

; FREE ENERGY
;----------------------------------------------------
free-energy              = yes
init-lambda-state        = MYLAMBDA
calc-lambda-neighbors    = -1
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30
vdw_lambdas              = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00
coul_lambdas             = 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00
couple-moltype           = $Lig_top_name
couple-lambda0           = vdw-q
couple-lambda1           = none
couple-intramol          = yes
separate-dhdl-file       = yes
nstdhdl                  = 100
sc-alpha                 = 0.5
sc-coul                  = no
sc-power                 = 1
sc-sigma                 = 0.3
" >> $WorkDir/Solvation_free_energy/Equilibration/npt.mdp

for ((i = 0 ; i < $NCONFS ; i++))
do
	mkdir $WorkDir/Solvation_free_energy/Equilibration/Lambda$i
	cd $WorkDir/Solvation_free_energy/Equilibration/Lambda$i
	
	# The following two commands will generate new files including npt.mdp, mdout.mdp and npt.tpr in each of Lambda subdirectories within $WorkDir/Solvation_free_energy/Equilibration directory
	sed 's/MYLAMBDA/'$i'/g' $WorkDir/Solvation_free_energy/Equilibration/npt.mdp > $WorkDir/Solvation_free_energy/Equilibration/Lambda$i/npt.mdp
	
	$GMX grompp -f $WorkDir/Solvation_free_energy/Equilibration/Lambda$i/npt.mdp \
		    -c $WorkDir/Solvation_free_energy/EM/Lambda$i/confout.gro \
		    -r $WorkDir/Solvation_free_energy/EM/Lambda$i/confout.gro \
		    -p $WorkDir/System_topology/system.top \
		    -o $WorkDir/Solvation_free_energy/Equilibration/Lambda$i/npt.tpr
	
done

# Equilibration on parallel
cd $WorkDir/Solvation_free_energy/Equilibration

# This command will generate multiple new files including confout.gro, ener.edr, traj_comp.xtc, dhdl.xvg, md.log, and state.cpt in each of the "Lambda" sub-directories in $WorkDir/Solvation_free_energy/Equilibration directory
srun $GMX mdrun -s npt.tpr -multidir $WorkDir/Solvation_free_energy/Equilibration/Lambda*


# Production
mkdir $WorkDir/Solvation_free_energy/Production
cd $WorkDir/Solvation_free_energy/Production

echo ";====================================================
; Equilibrium Simulations
;====================================================

; RUN CONTROL
;----------------------------------------------------
integrator   = sd            ; stochastic leap-frog integrator
nsteps       = 50000000      ; 2 * 50,000,000 fs = 100000 ps
dt           = 0.002         ; 2 fs
comm-mode    = Linear        ; remove center of mass translation
nstcomm      = 100           ; frequency for center of mass motion removal

; OUTPUT CONTROL
;----------------------------------------------------
nstxout                = 0         ; save coordinates to .trr every 100 ps
nstvout                = 0         ; save velocities to .trr every 100 ps
nstfout                = 0         ; save forces to .trr every 100 ps
nstxout-compressed     = 50000     ; xtc compressed trajectory output every 10 ps
compressed-x-precision = 1000      ; precision with which to write to the compressed trajectory file
nstlog                 = 10000     ; update log file every 20 ps
nstenergy              = 10000     ; save energies every 20 ps
nstcalcenergy          = 100       ; calculate energies every 100 steps (default=100)

; BONDS
;----------------------------------------------------
constraint_algorithm   = lincs      ; holonomic constraints
constraints            = h-bonds    ; constrain H-bonds
lincs-order            = 4
lincs-iter             = 1
lincs-warnangle        = 30         ; maximum angle that a bond can rotate before LINCS will complain (30 is default)
continuation           = yes        ; formerly known as 'unconstrained-start' - useful for exact continuations and reruns

; NEIGHBOR SEARCHING
;----------------------------------------------------
cutoff-scheme         = Verlet ; group or Verlet
verlet-buffer-tolerance = 0.005
ns-type               = grid   ; search neighboring grid cells
nstlist               = 10     ; 20 fs (default is 10)
rlist                 = 1.0    ; short-range neighborlist cutoff (in nm)
pbc                   = xyz    ; 3D PBC
periodic-molecules    = no

; ELECTROSTATICS & EWALD
;----------------------------------------------------
coulombtype      = PME                       ; Particle Mesh Ewald for long-range electrostatics
coulomb-modifier = Potential-shift-Verlet
rcoulomb         = 1.0                       ; short-range electrostatic cutoff (in nm)
rcoulomb-switch  = 0
ewald-geometry   = 3d                        ; Ewald sum is performed in all three dimensions
pme-order        = 4                         ; interpolation order for PME (default is 4)
fourierspacing   = 0.12                      ; grid spacing for FFT
ewald-rtol       = 1e-5                      ; relative strength of the Ewald-shifted direct potential at rcoulomb

; VAN DER WAALS
;----------------------------------------------------
vdw-type          = Cut-off      ; potential switched off at rvdw-switch to reach zero at rvdw
vdw-modifier      = Potential-shift-Verlet
rvdw              = 1.0          ; van der Waals cutoff (in nm)
DispCorr          = EnerPres     ; apply analytical long range dispersion corrections for Energy and Pressure

; TEMPERATURE COUPLING (Langevin)
;----------------------------------------------------
tcoupl            = no
tc-grps           = System
tau-t             = 1.0
ref-t             = 300
gen-vel           = no           ; Velocity generation (if gen-vel is 'yes', continuation should be 'no')
gen-temp          = 300
gen-seed          = -1

; PRESSURE COUPLING
;----------------------------------------------------
pcoupl           = C-rescale
pcoupltype       = isotropic
nstpcouple       = 10
tau_p            = 1.0                  ; time constant (ps)
ref_p            = 1.0
compressibility  = 4.5e-5
refcoord-scaling = no

; FREE ENERGY
;----------------------------------------------------
free-energy              = yes
init-lambda-state        = MYLAMBDA
calc-lambda-neighbors    = -1
; init_lambda_state        0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30
vdw_lambdas              = 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00
coul_lambdas             = 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00
couple-moltype           = $Lig_top_name
couple-lambda0           = vdw-q
couple-lambda1           = none
couple-intramol          = yes
separate-dhdl-file       = yes
nstdhdl                  = 100
sc-alpha                 = 0.5
sc-coul                  = no
sc-power                 = 1
sc-sigma                 = 0.3
" >> $WorkDir/Solvation_free_energy/Production/md.mdp

for ((i = 0 ; i < $NCONFS ; i++))
do
	mkdir $WorkDir/Solvation_free_energy/Production/Lambda$i
	cd $WorkDir/Solvation_free_energy/Production/Lambda$i
	
	# The following two commands will generate new files including md.mdp, mdout.mdp and md.tpr in each of Lambda subdirectories within $WorkDir/Solvation_free_energy/Production directory
	sed 's/MYLAMBDA/'$i'/g' $WorkDir/Solvation_free_energy/Production/md.mdp > $WorkDir/Solvation_free_energy/Production/Lambda$i/md.mdp
	
	$GMX grompp -f $WorkDir/Solvation_free_energy/Production/Lambda$i/md.mdp \
		    -c $WorkDir/Solvation_free_energy/Equilibration/Lambda$i/confout.gro \
		    -t $WorkDir/Solvation_free_energy/Equilibration/Lambda$i/state.cpt \
		    -p $WorkDir/System_topology/system.top \
		    -o $WorkDir/Solvation_free_energy/Production/Lambda$i/md.tpr

done


cd $WorkDir/Solvation_free_energy/Production

# This command will generate multiple new files including confout.gro, ener.edr, traj_comp.xtc, dhdl.xvg, md.log, state_prev.cpt and state.cpt in each of the "Lambda" sub-directories in $WorkDir/Solvation_free_energy/Production directory
srun $GMX mdrun -s md.tpr -multidir $WorkDir/Solvation_free_energy/Production/Lambda*



# Post-simulation analysis
cd $WorkDir/Solvation_free_energy/Production

# The command below will generate two new files named bar.xvg and barint.xvg in $WorkDir/Solvation_free_energy/Production directory
$GMX bar -f Production/Lambda0/dhdl.xvg Production/Lambda1/dhdl.xvg Production/Lambda2/dhdl.xvg Production/Lambda3/dhdl.xvg Production/Lambda4/dhdl.xvg Production/Lambda5/dhdl.xvg Production/Lambda6/dhdl.xvg Production/Lambda7/dhdl.xvg Production/Lambda8/dhdl.xvg Production/Lambda9/dhdl.xvg Production/Lambda10/dhdl.xvg Production/Lambda11/dhdl.xvg Production/Lambda12/dhdl.xvg Production/Lambda13/dhdl.xvg Production/Lambda14/dhdl.xvg Production/Lambda15/dhdl.xvg Production/Lambda16/dhdl.xvg Production/Lambda17/dhdl.xvg Production/Lambda18/dhdl.xvg Production/Lambda19/dhdl.xvg Production/Lambda20/dhdl.xvg Production/Lambda21/dhdl.xvg Production/Lambda22/dhdl.xvg Production/Lambda23/dhdl.xvg Production/Lambda24/dhdl.xvg Production/Lambda25/dhdl.xvg Production/Lambda26/dhdl.xvg Production/Lambda27/dhdl.xvg Production/Lambda28/dhdl.xvg Production/Lambda29/dhdl.xvg Production/Lambda30/dhdl.xvg -o -oi -temp 300











				##########################################################	
				# The output of the command above should look like this: #
				##########################################################

Production/Lambda0/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda0/dhdl.xvg: 0.0 - 100000.0; lambda = (0, 0)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda1/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda1/dhdl.xvg: 0.0 - 100000.0; lambda = (0.1, 0)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda2/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda2/dhdl.xvg: 0.0 - 100000.0; lambda = (0.2, 0)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda3/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda3/dhdl.xvg: 0.0 - 100000.0; lambda = (0.3, 0)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda4/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda4/dhdl.xvg: 0.0 - 100000.0; lambda = (0.4, 0)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda5/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda5/dhdl.xvg: 0.0 - 100000.0; lambda = (0.5, 0)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda6/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda6/dhdl.xvg: 0.0 - 100000.0; lambda = (0.6, 0)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda7/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda7/dhdl.xvg: 0.0 - 100000.0; lambda = (0.7, 0)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda8/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda8/dhdl.xvg: 0.0 - 100000.0; lambda = (0.8, 0)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda9/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda9/dhdl.xvg: 0.0 - 100000.0; lambda = (0.9, 0)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda10/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda10/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 0)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda11/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda11/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 0.05)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda12/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda12/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 0.1)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda13/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda13/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 0.15)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda14/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda14/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 0.2)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda15/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda15/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 0.25)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda16/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda16/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 0.3)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda17/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda17/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 0.35)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda18/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda18/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 0.4)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda19/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda19/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 0.45)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda20/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda20/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 0.5)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda21/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda21/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 0.55)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda22/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda22/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 0.6)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda23/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda23/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 0.65)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda24/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda24/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 0.7)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda25/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda25/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 0.75)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda26/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda26/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 0.8)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda27/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda27/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 0.85)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda28/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda28/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 0.9)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda29/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda29/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 0.95)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)


Production/Lambda30/dhdl.xvg: Ignoring set 'pV (kJ/mol)'.
Production/Lambda30/dhdl.xvg: 0.0 - 100000.0; lambda = (1, 1)
    dH/dl & foreign lambdas:
        dH/dl (coul-lambda) (500001 pts)
        dH/dl (vdw-lambda) (500001 pts)
        delta H to (0, 0) (500001 pts)
        delta H to (0.1, 0) (500001 pts)
        delta H to (0.2, 0) (500001 pts)
        delta H to (0.3, 0) (500001 pts)
        delta H to (0.4, 0) (500001 pts)
        delta H to (0.5, 0) (500001 pts)
        delta H to (0.6, 0) (500001 pts)
        delta H to (0.7, 0) (500001 pts)
        delta H to (0.8, 0) (500001 pts)
        delta H to (0.9, 0) (500001 pts)
        delta H to (1, 0) (500001 pts)
        delta H to (1, 0.05) (500001 pts)
        delta H to (1, 0.1) (500001 pts)
        delta H to (1, 0.15) (500001 pts)
        delta H to (1, 0.2) (500001 pts)
        delta H to (1, 0.25) (500001 pts)
        delta H to (1, 0.3) (500001 pts)
        delta H to (1, 0.35) (500001 pts)
        delta H to (1, 0.4) (500001 pts)
        delta H to (1, 0.45) (500001 pts)
        delta H to (1, 0.5) (500001 pts)
        delta H to (1, 0.55) (500001 pts)
        delta H to (1, 0.6) (500001 pts)
        delta H to (1, 0.65) (500001 pts)
        delta H to (1, 0.7) (500001 pts)
        delta H to (1, 0.75) (500001 pts)
        delta H to (1, 0.8) (500001 pts)
        delta H to (1, 0.85) (500001 pts)
        delta H to (1, 0.9) (500001 pts)
        delta H to (1, 0.95) (500001 pts)
        delta H to (1, 1) (500001 pts)



Temperature: 300 K

Detailed results in kT (see help for explanation):

 lam_A  lam_B      DG   +/-     s_A   +/-     s_B   +/-   stdev   +/- 
     0      1   -3.18  0.19    0.92  0.04    0.91  0.05    1.41  0.04
     1      2   -5.02  0.11    0.94  0.13    0.82  0.11    1.34  0.09
     2      3   -6.31  0.05    0.46  0.11    0.37  0.10    0.93  0.04
     3      4   -6.96  0.04    0.29  0.05    0.24  0.04    0.69  0.03
     4      5   -7.32  0.02    0.12  0.03    0.10  0.02    0.53  0.02
     5      6   -7.53  0.00    0.11  0.02    0.10  0.02    0.46  0.01
     6      7   -7.76  0.01    0.13  0.01    0.13  0.01    0.41  0.01
     7      8   -7.92  0.01    0.03  0.01    0.03  0.01    0.37  0.01
     8      9   -8.03  0.02    0.08  0.00    0.08  0.00    0.36  0.01
     9     10   -8.18  0.01    0.07  0.02    0.08  0.02    0.37  0.01
    10     11    2.66  0.01    0.41  0.01    0.50  0.01    0.95  0.00
    11     12    2.65  0.01    0.43  0.01    0.53  0.01    0.97  0.00
    12     13    2.61  0.02    0.44  0.02    0.56  0.03    0.99  0.01
    13     14    2.60  0.02    0.43  0.02    0.55  0.02    1.00  0.00
    14     15    2.58  0.01    0.47  0.01    0.59  0.01    1.02  0.00
    15     16    2.54  0.01    0.48  0.01    0.59  0.01    1.03  0.00
    16     17    2.52  0.01    0.48  0.01    0.59  0.01    1.04  0.00
    17     18    2.49  0.01    0.51  0.01    0.62  0.01    1.06  0.01
    18     19    2.45  0.01    0.50  0.01    0.61  0.01    1.06  0.00
    19     20    2.45  0.01    0.49  0.00    0.61  0.01    1.06  0.00
    20     21    2.45  0.01    0.52  0.01    0.63  0.01    1.08  0.01
    21     22    2.44  0.00    0.53  0.01    0.66  0.01    1.10  0.00
    22     23    2.37  0.00    0.58  0.01    0.73  0.01    1.15  0.01
    23     24    2.19  0.01    0.69  0.01    0.89  0.01    1.26  0.01
    24     25    1.85  0.02    0.82  0.01    1.09  0.01    1.41  0.01
    25     26    1.25  0.02    1.08  0.02    1.49  0.02    1.67  0.01
    26     27    0.04  0.03    1.59  0.02    2.12  0.02    2.09  0.02
    27     28   -1.38  0.02    1.74  0.02    1.91  0.02    2.12  0.02
    28     29   -1.50  0.02    1.10  0.01    1.06  0.01    1.54  0.01
    29     30   -0.23  0.00    0.53  0.00    0.50  0.00    1.03  0.00


Final results in kJ/mol:

point      0 -      1,   DG -7.93 +/-  0.46
point      1 -      2,   DG -12.53 +/-  0.28
point      2 -      3,   DG -15.73 +/-  0.12
point      3 -      4,   DG -17.37 +/-  0.11
point      4 -      5,   DG -18.27 +/-  0.05
point      5 -      6,   DG -18.78 +/-  0.01
point      6 -      7,   DG -19.35 +/-  0.04
point      7 -      8,   DG -19.75 +/-  0.04
point      8 -      9,   DG -20.02 +/-  0.05
point      9 -     10,   DG -20.39 +/-  0.03
point     10 -     11,   DG  6.63 +/-  0.02
point     11 -     12,   DG  6.61 +/-  0.01
point     12 -     13,   DG  6.52 +/-  0.05
point     13 -     14,   DG  6.48 +/-  0.05
point     14 -     15,   DG  6.43 +/-  0.02
point     15 -     16,   DG  6.34 +/-  0.01
point     16 -     17,   DG  6.28 +/-  0.02
point     17 -     18,   DG  6.20 +/-  0.01
point     18 -     19,   DG  6.11 +/-  0.02
point     19 -     20,   DG  6.11 +/-  0.01
point     20 -     21,   DG  6.12 +/-  0.03
point     21 -     22,   DG  6.08 +/-  0.01
point     22 -     23,   DG  5.91 +/-  0.00
point     23 -     24,   DG  5.46 +/-  0.03
point     24 -     25,   DG  4.63 +/-  0.04
point     25 -     26,   DG  3.11 +/-  0.06
point     26 -     27,   DG  0.10 +/-  0.06
point     27 -     28,   DG -3.44 +/-  0.05
point     28 -     29,   DG -3.75 +/-  0.04
point     29 -     30,   DG -0.58 +/-  0.00

total      0 -     30,   DG -82.77 +/-  0.72


