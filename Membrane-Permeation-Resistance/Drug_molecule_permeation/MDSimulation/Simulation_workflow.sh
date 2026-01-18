#!/bin/bash

GMX=gmx_mpi

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

set -e

WorkDir=/path/to/working/directory

# System coordinate and topology
System_coord_dir=/path/to/System_coordinates/folder
System_topology_dir=/path/to/System_topology/folder

NCONFS=21

cd $WorkDir


# Generate .mdp file
echo ";====================================================
; AWH Simulations
;====================================================

; RUN CONTROL
;----------------------------------------------------
integrator   = sd            ; stochastic leap-frog integrator
nsteps       = 500000000     ; 2 * 500,000,000 fs = 1000000 ps
dt           = 0.002         ; 2 fs
tinit        = 0.0
comm-mode    = Linear        ; remove center of mass translation
nstcomm      = 100           ; frequency for center of mass motion removal
comm_grps    = system

; OUTPUT CONTROL
;----------------------------------------------------
nstxout                = 0         ; save coordinates to .trr every 100 ps
nstvout                = 0         ; save velocities to .trr every 100 ps
nstfout                = 0         ; save forces to .trr every 100 ps
nstxout-compressed     = 5000      ; xtc compressed trajectory output every 10 ps
compressed-x-precision = 1000      ; precision with which to write to the compressed trajectory file
nstlog                 = 50000     ; update log file every 100 ps
nstenergy              = 5000      ; save energies every 10 ps
nstcalcenergy          = 100       ; calculate energies every 100 steps (default=100)
compressed-x-grps      =
energygrps             =

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
ref-t             = 303
gen-vel           = yes           ; Velocity generation (if gen-vel is 'yes', continuation should be 'no')
gen-temp          = 303
gen-seed          = -1

; PRESSURE COUPLING
;----------------------------------------------------
pcoupl           = c-rescale
pcoupltype       = semiisotropic
nstpcouple       = 10
tau_p            = 1.0                  ; time constant (ps)
ref_p            = 1.0 1.0
compressibility  = 4.5e-5 0
refcoord-scaling = all

; FREE ENERGY
;----------------------------------------------------
pull = yes
pull_ngroups = 2
pull_ncoords = 1
pull_nstfout = 100000
pull_nstxout = 50000
pull_fout_average = no
pull-pbc-ref-prev-step-com = yes

pull_group1_name = PullRef
pull_group2_name = YT0

pull_group1_pbcatom = 11451
pull_coord1_groups = 1 2
pull_coord1_type = external-potential
pull_coord1_potential_provider = awh
pull_coord1_dim = N N Y
pull_coord1_geometry = direction
pull_coord1_vec = 0 0 1
pull_coord1_start = yes

free-energy       = yes
couple-lambda0    = none
couple-lambda1    = vdwq
couple-moltype    = YT0
couple-intramol   = yes
init-lambda-state = 20
fep-lambdas       = 1.00 0.95 0.90 0.85 0.80 0.75 0.70 0.65 0.60 0.55 0.50 0.45 0.40 0.35 0.30 0.25 0.20 0.15 0.10 0.05 0.00

calc-lambda-neighbors    = -1
separate-dhdl-file       = no
sc_alpha                 = 0.5
sc_sigma                 = 0.3
sc_power                 = 1
sc_coul                  = yes
nstdhdl                  = 100

awh                      = yes
awh_nstout               = 100000
awh_potential            = umbrella
awh_share_multisim       = yes
awh_nbias                = 1
awh_nstsample            = 100
awh_nsamples-update      = 10
awh1_growth                = exp-linear
awh1_equilibrate_histogram = yes
awh1_target                = constant
awh1_user_data             = no
awh1_error_init            = 20
awh1_share_group           = 1
awh1_ndim                  = 2

awh1_dim1_coord_provider   = pull
awh1_dim1_coord_index      = 1
awh1_dim1_force_constant   = 25000
awh1_dim1_start            = -6.0
awh1_dim1_end              = 6.0
awh1_dim1_diffusion        = 1e-3 ; In nm^2/ps, i.e. multiply cm^2/s by 1e2. 5e-5 nm^2/ps = 5e-7 cm^2/s
awh1_dim1_cover-diameter   = 0.8

awh1-dim2-coord-provider   = fep-lambda
awh1-dim2-start            = 0
awh1-dim2-end              = 20
awh1-dim2-diffusion        = 1e-4
awh1-dim2-cover-diameter   = 21" > $WorkDir/awh.mdp

for ((i=1; i<$NCONFS; i++))
do
	mkdir $WorkDir/AWH_$i
	cd $WorkDir/AWH_$i
	
	# This command will generate file named awh.tpr and mdout.mdp within each subdirectory "AWH" within $WorkDir
	$GMX grompp -f $WorkDir/awh.mdp -p $System_topology_dir/system.top \
		    -c $System_coord_dir/conf_${i}.gro -maxwarn 1 \
		    -o $WorkDir/AWH_$i/awh.tpr -n $System_topology_dir/index.ndx
	
done

cd $WorkDir

# Executing the simulation
srun $GMX mdrun -s awh.tpr -cpi state.cpt -x traj_comp.xtc -c confout.gro -e ener.edr -g md.log -dhdl dhdl.xvg -multidir $WorkDir/AWH_*


# Post-simulation analysis procedure
# The analysis requires python environment
calFile=/path/to/water/solvation/free/energy/file/named/awh_t200000.xvg
time=1000000	# Specify the time where permeation simulation stops in ps
dimToUse=2
temp=303
diffusionCmd=/path/to/Analysis_scripts/directory/awh_diffusion.py
extractionCmd=/path/to/Analysis_scripts/directory/xvg_extract_dimension.py
permeabilityCmd=/path/to/Analysis_scripts/directory/awh_fep_calc_permeability.py

i=1
mkdir $WorkDir/AWH_$i/AWH-related-data_temp
cd $WorkDir/AWH_$i/AWH-related-data_temp 	

gmx_mpi awh -s $WorkDir/AWH_$i/awh.tpr -f $WorkDir/AWH_$i/ener.edr -quiet -more -fric $WorkDir/AWH_$i/AWH-related-data_temp/friction.xvg -o $WorkDir/AWH_$i/AWH-related-data_temp/awh.xvg

$extractionCmd -i $WorkDir/AWH_$i/AWH-related-data_temp/friction_t$time.xvg -c 1 -v 0 -o $WorkDir/AWH_$i/AWH-related-data_temp/extr_friction_t$time.xvg



for ((i=2; i<21; i++))
do
	mkdir $WorkDir/AWH_$i/AWH-related-data_temp
	cd $WorkDir/AWH_$i/AWH-related-data_temp
	gmx_mpi awh -s $WorkDir/AWH_$i/awh.tpr -f $WorkDir/AWH_$i/ener.edr -fric -b $time -e $time
	$extractionCmd -i $WorkDir/AWH_$i/AWH-related-data_temp/friction_t$time.xvg -c 1 -v 0 -o $WorkDir/AWH_$i/AWH-related-data_temp/extr_friction_t$time.xvg
	
done


cd $WorkDir
# calibrated PMF
i=1
$extractionCmd -i $WorkDir/AWH_$i/AWH-related-data_temp/awh_t$time.xvg -c 1 -v 0 -e $dimToUse --calibration_index 20 --calibrate_to_pmf_file $calFile -o $WorkDir/AWH_$i/AWH-related-data_temp/extr_pmf_watercal_t$time.xvg

$diffusionCmd -f $WorkDir/AWH_*/AWH-related-data_temp/extr_friction_t$time.xvg -t $temp -c 1 -s 0.2 -o $WorkDir/AWH_$i/AWH-related-data_temp/diffusion_t$time

$permeabilityCmd -p $WorkDir/AWH_$i/AWH-related-data_temp/extr_pmf_watercal_t$time.xvg -t $temp -l 14 -d $WorkDir/AWH_$i/AWH-related-data_temp/diffusion_t$time.xvg -o $WorkDir/AWH_$i/AWH-related-data_temp/permeability_watercal_t$time



















