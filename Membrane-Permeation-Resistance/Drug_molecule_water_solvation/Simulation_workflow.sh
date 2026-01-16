#!/bin/bash -l

set -e

WorkDir=/path/to/working/directory

# Solute coordinate and topology
Input_strct_coord=/path/to/drug/molecule/structure/YT0_GMX.gro
System_topology=/path/to/system/topology/system.top
# Solvent box and topology
Solvent_box=/path/to/pre-equilibrated/water/box/PreEquilibrated-water-box_100ns.gro

GMX=gmx_mpi
NCONFS=21


cd $WorkDir
# Generate simulation .mdp file
echo ";====================================================
; Equilibrium Simulations
;====================================================

; RUN CONTROL
;----------------------------------------------------
integrator   = sd            ; stochastic leap-frog integrator
nsteps       = 100000000     ; 2 * 100,000,000 fs = 200000 ps
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
free-energy              = yes		
couple-lambda0           = none		
couple-lambda1           = vdwq		
couple-moltype           = YT0
couple-intramol          = yes		
init-lambda-state        = 20		
fep-lambdas              = 1.00 0.95 0.90 0.85 0.80 0.75 0.70 0.65 0.60 0.55 0.50 0.45 0.40 0.35 0.30 0.25 0.20 0.15 0.10 0.05 0.00
calc-lambda-neighbors    = -1		
separate-dhdl-file       = no		
sc_alpha                 = 0.5
sc_sigma                 = 0.3
sc_power                 = 1
sc_coul                  = yes
nstdhdl                  = 2500000

awh                        = yes	
awh_nstout                 = 2500000	
awh_potential              = umbrella	
awh_share_multisim         = yes	
awh_nbias                  = 1		
awh_nstsample              = 100	
awh_nsamples-update        = 10		
awh1_growth                = exp-linear	
awh1_equilibrate_histogram = yes	
awh1_target                = constant	
awh1_user_data             = no	
awh1_error_init            = 10
awh1_share_group           = 1		
awh1_ndim                  = 1		
awh1-dim1-coord-provider   = fep-lambda	
awh1-dim1-start            = 0		
awh1-dim1-end              = 20		
awh1-dim1-diffusion        = 1e-4
awh1-dim1-cover-diameter   = 21		
" > $WorkDir/awh.mdp


for ((i=1; i<$NCONFS; i++))
do
	mkdir $WorkDir/AWH_$i
	cd $WorkDir/AWH_$i
	
	$GMX insert-molecules -f $Solvent_box -ci $Input_strct_coord -o $WorkDir/AWH_$i/conf_${i}.gro -nmol 1 -scale 0.1 -try 10000
	$GMX grompp -f $WorkDir/awh.mdp -p $System_topology -c $WorkDir/AWH_$i/conf_${i}.gro -maxwarn 1 -o $WorkDir/AWH_$i/awh.tpr
	
done

# Run production
cd $WorkDir

srun $GMX mdrun -s awh.tpr -cpi state.cpt -x traj_comp.xtc -c confout.gro -e ener.edr -g md.log -dhdl dhdl.xvg -multidir $WorkDir/AWH_*






















