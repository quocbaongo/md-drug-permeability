# Drug Permeability Estimation via Molecular Dynamics Simulations

This repository provides a reproducible computational workflow for theoretically estimating the **permeability coefficient** of small drug-like molecules across a lipid bilayer membrane structure using **molecular dynamics (MD) simulations**.

The workflow follows a pipeline of **Parameterization → Validation → Production**. The Inhomogeneous Solubility-Diffusion (ISD) model was used to compute permeation resistance across the lipid bilayer, using the drug molecule **YT0** (RCSB Ligand ID: [YT0](https://www.rcsb.org/ligand/YT0)) as the primary case study.






## Workflow Overview

1. **Drug Molecule Parameterization**  
   Generate simulation-ready topology and coordinate files using quantum mechanics for partial charges and the General Amber Force Field (GAFF) for bonded parameters.

2. **Parameter Validation – log Pₒw Calculation**  
   Compute solvation free energies in water and octanol → calculate log Pₒw as a quality check of the force field parameters.

3. **Permeability Coefficient Calculation**  
   Perform umbrella sampling across a lipid bilayer → reconstruct PMF (ΔG(z)) → compute local diffusion coefficients D(z) → integrate to obtain permeation resistance (1/P).

## Software & Dependencies

| Software / Library       | Version       | 
|--------------------------|---------------|
| GROMACS                  | 2022.4        | 
| AmberTools               | 2024          | 
| Gaussian                 | 16            |
| Python                   | 3.9+          |
| NumPy                    | 1.24.3        |
| Scipy                    | 1.5.3         |
| Sklearn                  | 1.3.0         |
| Statsmodels              | 0.14.1        |
| Uncertainties            | 3.2.3         |

## Repository Structure (suggested / typical layout)

├── README.md
├── Parameterization-workflow.sh           # Main script: YT0 parameterization
├── logP_calculation/
│   ├── water/
│   ├── octanol/
│   └── analyze_logP.py
├── permeability/
│   ├── bilayer_setup/
│   ├── umbrella_windows/
│   ├── analysis/
│   └── permeability_integration.py
├── data/
│   ├── YT0_Gaussian.log
│   ├── topologies/
│   └── example_profiles/
└── figures/
└── example_free_energy_profile.png

## Detailed Workflow

### 1. Drug Molecule Parameterization

See [`Parameterization-workflow.sh`](Parameterization-workflow.sh)

This script guides through:

- 3D structure preparation (RDKit or manual)
- Geometry optimization + HF/6-31G* ESP calculation in Gaussian 16
- RESP charge fitting
- GAFF parameter assignment using `antechamber` and `parmchk2`
- Conversion to GROMACS format (via acpype or parmed)

### 2. Water–Octanol Partition Coefficient (Validation)

Solvation free energy difference is used to compute:

$$
\log_{10} P_{ow} = \frac{\Delta G^\circ_\text{water} - \Delta G^\circ_\text{octanol}}{RT \ln(10)}
$$

where  
- $R = 8.31446261815 \times 10^{-3}$ kJ mol⁻¹ K⁻¹  
- $T = 298$ K (standard)  
- $\ln(10) \approx 2.302585$

Typical methods: thermodynamic integration (TI) or Bennett Acceptance Ratio (BAR) in explicit solvent.

### 3. Membrane Permeability Coefficient

Follows the position-dependent permeability model (see reference below).

The inverse permeability (resistance) is calculated as:

$$
R = \frac{1}{P} = (30 \pm 6) \times \int_{z_1}^{z_2} \frac{e^{\beta \Delta G_\text{rel, water}(z)}}{D(z)}\ \text{d}z
$$

where  
- $\beta = 1/(k_B T)$  
- $\Delta G_\text{rel, water}(z)$ = potential of mean force relative to bulk water  
- $D(z)$ = local diffusion coefficient along the membrane normal  
- The factor $(30 \pm 6)$ is an empirical prefactor (in s/m or equivalent units after calibration; commonly used in recent literature)

Reference methodology:  
Lundborg et al., *Biophysical Journal* (2022)  
https://doi.org/10.1016/j.bpj.2022.07.016

## Getting Started

1. Clone the repository  
   ```bash
   git clone https://github.com/yourusername/drug-permeability-md-workflow.git
   cd drug-permeability-md-workflow
