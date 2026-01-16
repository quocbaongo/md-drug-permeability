. This working folder aims to describe workflow to compute logP value of research drug molecule using molecular dynamics simulation. As log P, specifically refers to the ratio of a compound in its neutral (unionized) state when in equilibrium between an organic solvent (typically n-octanol) and water. Thus the drug molecule in neutral form is chosen for this computation.

. As described here https://pmc.ncbi.nlm.nih.gov/articles/PMC10101867: 

	logP = (ΔG of water hydration free energy - ΔG of octanol solvent free energy) / RTln10
	
hence, we need to obtain hydration free energy and octanol solvation free energy of researched molecule, which is feasible using molecular dynamics simulation.

. The two sub-folders named "Water_hydration-in-parallel-abcg2_gaff2" and "Octanol_hydration-in-parallel-abcg2_gaff2" present workflow for computing water and octanol like hydration free energy. Command line details can be found in the file named "Simulation_workflow.sh" in each folder. The folder "Solvation_free_energy_results" within each sub-folder presents raw output generated from the solvation free energy computation procedure.

. It is important to note that the the output values of the computations are named water and octanol like hydration free energy as they do not truly represent water and octanol hydration free energy. But the resulting numbers still serve the purpose of calculating logP of researched drug molecule. More detailed explanation for this is available in the powerpoint slide "logP-thermodynamics-cycle"

. In general the computational workflow to computing either the water and octanol like hydration free energy can be described in three steps:

	Step1: Establishing simulation system and running short MD simulation for system equilibration purpose. The step is pretty much like a standard MD simulation
	Step2: In this step, we aimed firstly to gradually turn off the electrostatics interactions of drug molecule (both intramolecular and extramolecular electrostatics interactions), which was followed by gradually turning off van der Waals interactions (both intramolecular and extramolecular vdw interactions). After we turned off both the intramolecular and extramolecular vdw + electrostatics interactions, the drug molecule behaves like a dummy molecule, where only its bonded parameters keep the molecule from falling apart. We also had two turning off trajectories after turning off all the non-bonded interacted, which are "$WorkDir/Starting_conformations_generation/Electrostatics_decoupling/traj_comp.xtc" and "$WorkDir/Starting_conformations_generation/vdw_decoupling/traj_comp.xtc". Along those two trajectories, we then extracted multiple conformations of drug molecule, and each of those conformations represents the structure of drug, when its non-bonded interaction parameters are scaled to a specific degree.
	
		To explain in further detail, we have a concept of lambda parameter. When vdw_lambdas and electrostatics_lambdas equal to 0, both the vdw and electrostatics interactions are fully on. When the two lambdas equal to 1.0, both the interactions are off. Similarly, if vdw_lambdas = 0.15 and electrostatics_lambdas = 1.00, it means both the intramolecular + extramolecular electrostatics interaction of drug molecule are totally off, while its vdw interaction strength are scaled 15% down from the original values. Below are the extracted drug conformations with different scaled vdw+electrostatics interactions strength. Each of these conformation would serve as starting structure for more thoroughly sampling simulations (described in step 3).
		
		
		
					lambda index	vdw_lambdas	electrostatics_lambdas
		Conformation 1		0		0.00		0.00				# All the non-bonded interactions are on
		Conformation 2		1		0.00		0.10
		Conformation 3		2		0.00		0.20
		Conformation 4		3		0.00		0.30
		Conformation 5		4		0.00		0.40
		Conformation 6		5		0.00		0.50
		Conformation 7		6		0.00		0.60
		Conformation 8		7		0.00		0.70
		Conformation 9		8		0.00		0.80			
		Conformation 10		9		0.00		0.90
		Conformation 11		10		0.00		1.00				# Totally turning off electrostatics interactions
		Conformation 12		11		0.05		1.00
		Conformation 13		12		0.10		1.00
		Conformation 14		13		0.15		1.00
		Conformation 15		14		0.20		1.00
		Conformation 16		15		0.25		1.00
		Conformation 17		16		0.30		1.00
		Conformation 18		17		0.35		1.00
		Conformation 19		18		0.40		1.00
		Conformation 20		19		0.45		1.00
		Conformation 21		20		0.50		1.00
		Conformation 22		21		0.55		1.00
		Conformation 23		22		0.60		1.00
		Conformation 24		23		0.65		1.00
		Conformation 25		24		0.70		1.00
		Conformation 26		25		0.75		1.00
		Conformation 27		26		0.80		1.00
		Conformation 28		27		0.85		1.00
		Conformation 29		28		0.90		1.00
		Conformation 30		29		0.95		1.00
		Conformation 31		30		1.00		1.00				# All the non-boned interactions are off


	Step3: Each of the conformations above will be used as a starting structure for a 100ns standard MD simulation. The parameter "init-lambda-state" in .mdp file indicates the index of lambda state (base 0)


























