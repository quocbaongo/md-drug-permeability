# Drug Permeability Estimation via Molecular Dynamics Simulations

This repository provides a reproducible computational workflow for theoretically estimating the **permeability coefficient** of small drug-like molecules across a lipid bilayer membrane structure using **molecular dynamics (MD) simulations**. The workflow follows a pipeline of **Deriving drug molecule simulation parameters → validating simulation parameters → computing drug's resistance permeation**. The Inhomogeneous Solubility-Diffusion (ISD) model was used to compute permeation resistance of the drug molecule across the lipid bilayer.

## Software & Dependencies

**MD simulation & Quantum Mechanics:**
* **GROMACS:** v2022.4
* **AmberTools:** v2024
* **Gaussian:** v16
* **PyMol:** v3.1.8

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
- **Charge Derivation**: Atomic partial charges were assigned using the **Restrained Electrostatic Potential (RESP)** method:

  **QM Geometry Optimization** (Gaussian) **->** **Electrostatic Potential (ESP) Calculation** on the optimized structure (Gaussian) **->** **RESP Charge Fitting** (AmberTools)

  The drug molecule was first geometry-optimized at the quantum mechanical level, after which the electrostatic potential surrounding the optimized structure was computed, both in Gaussian. The resulting ESP was then used in AmberTools to fit RESP atomic partial charges that reproduce the molecule's electrostatic potential at its molecular surface.
  
- **Topology Generation**: Deriving bonded parameters using the General Amber Force Field 2 (GAFF2). Antechamber successfully assigned parameters for all bond terms; however, two dihedral angles were flagged with a high penalty score and required validation:

  **QM Torsional Scan** (0° -> 360°, 5° step, remaining structure minimized at each point) **->** **QM Torsional Energy Profile** (reference) **->** **Compare against MD-Sampled Dihedral Distribution** (6 replicates x 300 ns)
  For each flagged dihedral, the angle was rotated from 0° to 360° in 5° increments, minimizing the rest of the structure at each step to generate a set of conformers. The energy of each conformer was computed to construct a QM torsional energy profile, which served as the reference. The derived partial charges and GAFF2 bond parameters were then used to define the drug molecule in standard MD simulations (6 replicates, 300 ns each), and the probability distribution of the sampled dihedral angle was compared against the reference QM torsional energy profile to confirm reasonable conformational sampling.

The two videos below present the outcome of the dihedral angle validation process described above. Each video below is organized into three panels:
- **Top row**: the dihedral rotation, generating each conformer along the scan.
- **Bottom-left**: the QM relative energy corresponding to each rotameric conformer displayed above.
- **Bottom-right**: a static graph showing the probability distribution of the inspected dihedral angle, obtained from standard molecular dynamics simulations.
<br><br><br>
<p align="center"><strong> First Dihedral Validation </strong></p>
<div align="center">
	<video src="https://github.com/user-attachments/assets/3ad3d8c6-0684-4b04-9948-1d836c0730b5" controls width="600"></video>
</div>
<br><br><br>
<p align="center"><strong> Second Dihedral Validation </strong></p>
<div align="center">
	<video src="https://github.com/user-attachments/assets/184d6711-5ed7-45d3-8e08-ce8fa22dbfb9" controls width="600"></video>
</div>
<br>

Overall, the GAFF2 bonded parameters are considered reliable. For the first dihedral, the QM global energy minimum corresponds to the most populated angle in the MD simulation, indicating good agreement between the QM and MD-derived profiles.

For the second dihedral , however, the most populated angle in the MD simulation outcome corresponds to the starting structure's dihedral value rather than the QM energy minimum. This result is not surprised as the simulated system need time to get out of the local energy minima to discover lower energy state.

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


Passive permeation through lipid bilayers membrane is a rare event, and the main rate-limiting step is simulated molecules being flip-flop inside the membrane and their entry/exit at the headgroup interface. Unbiased MD simulation (see below), might not sample sufficient permeation of simulated molecules within practical timescales.

<p align="center"><strong> Unbiased MD simulation of drug permeating membrane </strong></p>
<div align="center">
	<video src="https://github.com/user-attachments/assets/df0e4e8b-847c-4395-8639-0f8251c5023f" controls width="200"></video>
</div>
<br>

Unbiased ATMD simulations of drug molecule permeating lipid bilayer membrane. In total, ~ 4000ns of unbiased ATMD simulation (10 replications and 400 ns each) was conducted. Still, no permeation event could be observed. Note: sudden jump of drug molecule's z-position across the membrane from one end of the simulation box
to the other end is due to implementation of periodic boundary condition for approximating a large (infinite) system in MD simulation code. The underlying idea of the algorithm is that after an object passes through one side of the cell, it reappears on the opposite side with the same velocity. The lipid bilayer spans from 4-8 nm in z-direction​. Original video is compressed to 10Mb to ful-fill GitHub requirement, check [OneDrive link](https://1drv.ms/v/c/7232c20154625746/IQDwyWRYEu4SQL7LPFkHT72xAXyQa-e_Zetd6ra1fwxpBXc?e=mdobNg) for original video.

<br>

Enhanced sampling method Accelerated Weight Histogram (AWH) algorithm was used to flatten underlying free energy barriers due to the membrane and enable multiple permeation events within a single simulation.

<p align="center"><strong> AWH MD simulation of drug permeating membrane </strong></p>
<div align="center">
	<video src="https://github.com/user-attachments/assets/937e0cee-af87-47ad-8b3d-2f279a1e27e6" controls width="200"></video>
</div>
<br>

(A) Simulation trajectory of an AWH walker, 
(B) the recorded z-positions of drug’s center of mass over the simulation shown in (A), 
(C) evolution of position-dependent potential of mean force (PMF) profile,
(D) Evolution of position-dependent diffusivity profile.​

The lipid bilayer spans from 4-8 nm in z-direction​. The permeability coefficient of simulated drug can be obtained using the Inhomogeneous solubility-diffusion (ISD) model, which integrates the PMF profile in (C) and position dependent diffusivity profile in (D). Original video is compressed to 10Mb to ful-fill GitHub requirement, check [OneDrive link](https://1drv.ms/v/c/7232c20154625746/IQBWu0Kex-KBS4UXXdRDqFVwAaBBvI2Y6QyfafIDgkjLoHE?e=U0blV0) for original video.


<br>
Sampled drug conformations when its non-bonded interactions are fully ON are shown at different positions relative to the membrane: bulk water, membrane–water interface (adsorption), and inside membrane. Note that the simulation box is ~12 nm along z-direction, with the membrane spanning between z = 4–8 nm. Structures shown are the dominant conformations from RMSD clustering pooled across 20 AWH walkers.



https://github.com/user-attachments/assets/6f6e666c-0341-4989-889e-655da50b9ea5


