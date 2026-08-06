# Drug Permeability Estimation via Molecular Dynamics Simulations

This repository provides a reproducible computational workflow for theoretically estimating the **permeability coefficient** of small drug-like molecules across a lipid bilayer membrane structure using **molecular dynamics (MD) simulations**. The workflow follows a pipeline of **Deriving drug molecule simulation parameters → validating simulation parameters → computing drug's resistance permeation**. The Inhomogeneous Solubility-Diffusion (ISD) model was used to compute permeation resistance of the drug molecule across the lipid bilayer.

## Software & Dependencies

**MD simulation & Quantum Mechanics:**
* **GROMACS:** v2022.4
* **AmberTools:** v2024
* **Gaussian:** v16
* * **PyMol:** v3.1.8

**Python v3.9+ Environment:**
* `numpy`: 1.24.3
* `scipy`: 1.5.3
* `sklearn`: 1.3.0 
* `statsmodels`: 0.14.1
* `uncertainties`: 0.14.1

## Repository Structure (suggested / typical layout)

```text
├── Drug-molecule-parameters/                            # Generating simulation parameters for researched drug molecule
│   ├── Drug_molecule_structure/                         # Structure of the drug molecule in .mol2 format
│   ├── Parameterization_workflow/
│   	├── Parameterization-workflow.sh                 # Description of generating simulation parameters for drug molecule using Ambertools and Gaussian softwares
│   	├── 1.Geometry_optimization/
│   	├── 2.ESP-charges-calculation/
│   	├── 3.RESP-calculation/
│   	└── 4.Drug_Gromacs_parameters/
│
├── Water-octanol-partitioning/                           # Validating generated simulation parameters   	
│   ├── logP-thermodynamics-cycle.pptx                    # Theoretical background underlying logP computation
│   ├── Water_solvation/
│   	└── Simulation_workflow.sh			              # Detailed illustration of the solvation free energy for the drug in water 
│   └── Octanol_solvation/
│   	└── Simulation_workflow.sh			              # Detailed illustration of the solvation free energy for the drug in octanol
│
├── Resistance-permeation-computation/			          # Estimating permeability coefficient of the researched drug molecule
│   ├── solvent-to-membrane-thermodynamics-cycle.pptx	  # Theoretical background underlying the resistance permation computation
│   ├── Drug_molecule_water_solvation/
│   	├── Simulation_workflow.sh			              # Detailed illustration of the solvation free energy for the drug in water using AWH method
│	    └── PreEquilibrated-water-box_100ns.gro		      # Pre-equilibrated box of water molecules
│   └── Drug_molecule_permeation/
│   	├── System_topology/				              # Drug molecule simulation's parameters
│   	├── System_coordinates/				              # Drug molecule starting coordinates
│   	└── MDSimulation/
   	        ├── Simulation_workflow.sh			          # Detailed illustration of the drug molecule's permeation through lipid bilayer
	    	└── Analysis_scripts/			              # Python scripts for post-simulation analysis written by Lundborg et al. (2024)
```
## Detailed Workflow

### 1. Drug Molecule Parameterization

See [`Parameterization-workflow.sh`](https://github.com/quocbaongo/md-drug-permeability/blob/main/Drug-molecule-parameters/Parameterization_workflow/Parameterization-workflow.sh)

This script guides through generating simulation parameters for the researched drug molecule. The procedure involves:
- **Charge Derivation**: Assigning partial atomic charges using Quantum Mechanics (QM) methods.

**QM Geometry Optimization** (Gaussian) **->** **Electrostatic Potential (ESP) Calculation** on the optimized structure (Gaussian) **->** **RESP Charge Fitting** (AmberTools)

The drug molecule was first geometry-optimized at the quantum mechanical level,
after which the electrostatic potential surrounding the optimized structure was
computed, both in Gaussian. The resulting ESP was then used in AmberTools to
fit RESP atomic partial charges that reproduce the molecule's electrostatic
potential at its molecular surface.

- **Topology Generation**: Deriving bonded parameters using the General Amber Force Field (GAFF).


<div align="center">
	<video src="https://github.com/user-attachments/assets/184d6711-5ed7-45d3-8e08-ce8fa22dbfb9" controls width="600"></video>
</div>









### 2. Water–Octanol Partitioning

See [`Drug-solvation-free-energy-in-water-workflow.sh`](https://github.com/quocbaongo/md-drug-permeability/blob/main/Water-octanol-partitioning/Water_solvation/Simulation_workflow.sh) and [`Drug-solvation-free-energy-in-octanol-workflow.sh`](https://github.com/quocbaongo/md-drug-permeability/blob/main/Water-octanol-partitioning/Octanol_solvation/Simulation_workflow.sh)

Basic validation of the simulation parameters can be performed by reproducing the experimental water-octanol partition coefficient (logP). The estimation of drug molecule's logP value requires the computation of its solvation free energies in water ($$ΔG_{hydration}$$​) and octanol ($$ΔG_{octanol, solvation}$$​). The partition coefficient is calculated as:

$$
logP = \frac{ΔG_{hydration}​ - ΔG_{octanol, solvation}}{RT \ln(10)}
$$

where  
- $R = 8.31446261815 \times 10^{-3}$ kJ mol⁻¹ K⁻¹  
- $T = 298$ K (standard)  
- $\ln(10) \approx 2.302585$

### 3. Membrane Permeation Resistance

See [`Drug-membrane-permeation-workflow.sh`](https://github.com/quocbaongo/md-drug-permeability/blob/main/Membrane-Permeation-Resistance/Drug_molecule_permeation/MDSimulation/Simulation_workflow.sh)

This section provides the workflow for computing the permeability coefficient of the drug molecule passing through a lipid bilayer membrane. Here, we utilize the Inhomogeneous Solubility-Diffusion (ISD) model. The total resistance to permeation (R) is derived by integrating the local resistance across the membrane, which is placed perpendicular to z-axis:

$$
R = \frac{1}{K_{p}} = N \times \int_{z_1}^{z_2} \frac{e^{\beta \Delta G_\text{rel, water}(z)}}{D(z)}\ \text{d}z
$$

where  
- $\beta = 1/(k_B T)$  
- $\Delta G_\text{rel, water}(z)$ is the Potential of Mean Force (free energy profile) relative to the hydration free energy.
- $D(z)$ is the local diffusion coefficient profile across the membrane. 
- N is the number of the lipid bilayers.
- $$K_{p}$$ is the permeability coefficient.
<br>
<br>
Below present the evolution of free energy profile and local diffusion coefficient profile throughout the simulation time


https://github.com/user-attachments/assets/80f57717-cf14-4ff6-ad3f-0827a14a3e99




https://github.com/user-attachments/assets/c3510ea3-05b1-4ffa-8d89-808f708c592b

<br>
Sampled drug conformations when its non-bonded interactions are fully ON are shown at different positions relative to the membrane: bulk water, membrane–water interface (adsorption), and inside membrane. Note that the simulation box is ~12 nm along z-direction, with the membrane spanning between z = 4–8 nm. Structures shown are the dominant conformations from RMSD clustering pooled across 20 AWH walkers.



https://github.com/user-attachments/assets/6f6e666c-0341-4989-889e-655da50b9ea5


