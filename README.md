# Drug Permeability Estimation via Molecular Dynamics Simulations

This repository provides a reproducible computational workflow for theoretically estimating the **permeability coefficient** of small drug-like molecules across a lipid bilayer membrane structure using **molecular dynamics (MD) simulations**. The workflow follows a pipeline of **Deriving drug molecule simulation parameters → validating simulation parameters → computing drug's resistance permeation**. The Inhomogeneous Solubility-Diffusion (ISD) model was used to compute permeation resistance across the lipid bilayer, using the drug molecule **YT0** (RCSB Ligand ID: [YT0](https://www.rcsb.org/ligand/YT0)) as the primary case study.

## Software & Dependencies

**MD simulation & Quantum Mechanics:**
* **GROMACS:** v2022.4
* **AmberTools:** v2024
* **Gaussian:** v16

**Python v3.9+ Environment:**
* `numpy`: 1.24.3
* `scipy`: 1.5.3
* `sklearn`: 1.3.0 
* `statsmodels`: 0.14.1
* `uncertainties`: 0.14.1

## Repository Structure (suggested / typical layout)

```text
├── Drug-molecule-parameters/                        # Generating simulation parameters for researched drug molecule
│   ├── Drug_molecule_structure/                     # Structure of drug molecule "YT0_H.mol2" in .mol2 format
│   ├── Parameterization_workflow/
│   	├── Parameterization-workflow.sh               # Description of generating simulation parameters for drug molecule using Ambertools and Gaussian softwares
│   	├── 1.Geometry_optimization/
│   	├── 2.ESP-charges-calculation/
│   	├── 3.RESP-calculation/
│   	└── 4.Drug_Gromacs_parameters/
│
├── Water-octanol-partitioning/                       # Validating generated simulation parameters   	
│   ├── logP-thermodynamics-cycle.pptx                # Theoretical background underlying logP computation
│   ├── Water_solvation/
│   	└── Simulation_workflow.sh			               # Detailed illustration of the solvation free energy for the drug in water 
│   └── Octanol_solvation/
│   	└── Simulation_workflow.sh			               # Detailed illustration of the solvation free energy for the drug in octanol
│
├── Resistance-permeation-computation/			         # Estimating permeability coefficient of the researched drug molecule
│   ├── solvent-to-membrane-thermodynamics-cycle.pptx	# Theoretical background underlying the resistance permation computation
│   ├── Drug_molecule_water_solvation/
│   	├── Simulation_workflow.sh			               # Detailed illustration of the solvation free energy for the drug in water using AWH method
│	└── PreEquilibrated-water-box_100ns.gro		      # Pre-equilibrated box of water molecules
│   └── Drug_molecule_permeation/
│   	├── System_topology/				                  # Drug molecule simulation's parameters
│   	├── System_coordinates/				               # Drug molecule starting coordinates
│   	└── MDSimulation/
   	        ├── Simulation_workflow.sh			         # Detailed illustration of the drug molecule's permeation through lipid bilayer
			├── Analysis_procedure.sh			               # Post-simulation analysis procedure
	    	└── Analysis_scripts/			                  # Python scripts for post-simulation analysis written by Lundborg et al. (2024)
```
## Detailed Workflow

### 1. Drug Molecule Parameterization

See [`Parameterization-workflow.sh`](Parameterization-workflow.sh)

This script guides through generating simulation parameters for the researched drug molecule. The procedure involves:
- **Charge Derivation**: Assigning partial atomic charges using Quantum Mechanical (QM) methods. 
- **Topology Generation**: Deriving bonded parameters using the General Amber Force Field (GAFF).

### 2. Water–Octanol Partitioning

Basic validation of the simulation parameters can be performed by reproducing the experimental water-octanol partition coefficient (logP). The estimation of drug molecule's logP value requires the computation of its solvation free energies in water (ΔGw​) and octanol (ΔGo​). The partition coefficient is calculated as:

$$
\log_{10} P_{ow} = \frac{\Delta G^\circ_\text{water} - \Delta G^\circ_\text{octanol}}{RT \ln(10)}
$$

where  
- $R = 8.31446261815 \times 10^{-3}$ kJ mol⁻¹ K⁻¹  
- $T = 298$ K (standard)  
- $\ln(10) \approx 2.302585$

### 3. Membrane Resistance Permeation

This section provides the workflow for computing the permeability coefficient of the drug molecule passing through a lipid bilayer membrane. Here, we utilize the Inhomogeneous Solubility-Diffusion (ISD) model. The total resistance to permeation (R) is derived by integrating the local resistance across the membrane, which is placed perpendicular to z-axis:

$$
R = \frac{1}{P} = N \times \int_{z_1}^{z_2} \frac{e^{\beta \Delta G_\text{rel, water}(z)}}{D(z)}\ \text{d}z
$$

where  
- $\beta = 1/(k_B T)$  
- $\Delta G_\text{rel, water}(z)$ is the Potential of Mean Force (free energy profile) relative to the hydration free energy.
- $D(z)$ is the local diffusion coefficient profile across the membrane. 
- N is the number of the lipid bilayers.
- P is the permeability coefficient.

