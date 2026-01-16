. The sub-folder "Drug_molecule_permeation" and "Drug_molecule_hydration" contains structure+topology file and computational workflow for performing simulation of drug's permeation process through the DMPC membrane and hydration free energy for water calibration purpose, respectively (will be explained in more detail below). The enhanced sampling simulation technique named "Accelerated weight histogram" (AWH) was used for both of the simulations.

	
	. Within "Drug_molecule_permeation" folder: 
	
		. The folder "System_coordinates" contains multiple starting coordinates (files with names starting with "conf") for the permeation AWH simulation, where the drug molecule is placed randomly along the reaction coordinate Z-axis. Note that Z-axis is perpendicular to the surface of DMPC membrane => the movement of drug molecule along the axis corresponds to its movement from one water bulk phase to the another one through the membrane.
	
		. The folder "System_topology" contains the topology files for all the molecules (drug+DMPC molecule+water) in the simulated system.
		
		. The folder "MDSimulation" contains file "Simulation_workflow.sh" describing step-by-step how to perform the AWH simulation for drug molecule permeation process and all the python scripts within "Analysis_scripts" folder were used for post-simulation analysis. 


	. Within "Drug_molecule_water_solvation" folder:

		. File named "PreEquilibrated-water-box_100ns.gro" is the pre-equilibrated water box, in which drug molecule will be placed randomly inside the box.
		. File "Simulation_workflow.sh" describes the step-by-step how to compute water solvation free energy using AWH simulation technique.



. AWH strategy illustrations for computing drug permeation procedure and water solvation free energy. It is important to note that in the scope of this study, we assumed that only the neutral state of researched drug molecule can pass through the DMPC membrane. Thus, we did not perform the permeation simulation and water solvation free energy calculation of drug molecule in charged state (even though it was predicted to be the dominant state in physiological ph, ph=7.5). 

. The assumption was made due to an observation in https://pubs.acs.org/doi/10.1021/acs.jcim.4c00722 that the permeability coeffcients computed based on the assumption that only neutral states can pass through the membrane is similar to those based on the assumption that both the neutral and charged states can pass through.

		######################################## Drug permeation through DMPC membrane ################################################ 

		The AWH drug permeation through DMPC membrane procedure was performed with 20 walkers. At the start of each AWH walker simulation, the drug molecule was randomly inserted along Z-axis (check file name starting with "conf" in "Drug_molecule_permeation/System_coordinates" folder). This means that the simulations were run in parallel, and the sampling was communicated between them when the bias and potential of mean force (PMF) are updated. That means that they share the same bias and effectively share the sampling of the free energy landscape. This sharing allows the sampling to be more effective, enhances the convergence of PMF profile and altogether reduces the required simulation time for a converged PMF profile.
		
		The drug molecule was inserted along the z-axis at a random position with all interactions with its environment turned off. The PMF profile through DMPC membrane structure was calculated using a two-dimensional AWH setup, using a harmonic potential to steer the drug across the membrane (Z-axis), and an alchemical free energy reaction coordinate. The alchemical free energy calculations imply that interactions of the drug with its surroundings can be gradually decoupled (switched on or off) => estimating the free energy of insertion (same way as calculating the solvation free energy).
		
		A major advantage of two dimension AWH is that this allows sampling the free energy along the permeation direction and also the relative insertion free energy of the drug from the vacuum => enables a calibration to the hydration free energy since the vacuum state is the same in both cases. After calibration, each point in the PMF corresponds to the free energy of transfer from the water vehicle to that point of the DMPC membrane structure. The thermodynamic cycle, which is used as a base theory for the water calibration is shown in the powerpoint slide "solvent-to-membrane-thermodynamics-cycle.pptx".
		
		Another advatage of this two dimension AWH approach is that this also improves the  sampling in regions with very slow diffusion (long correlation times), since turning off the interactions with the surroundings will let the drug leave that region to sample other parts of the landscape, and later return in a different configuration.
		


		############################################## Water solvation free energy ####################################################
		
		The solvation free energy like computation of drug molecule starts with a water box large enough to solvate the molecule. The drug molecule was then inserted into the system at a random position with its interactions to the surroundings turned off, as if in vacuum. The AWH algorithm was used to sample the alchemical free energy λ states to calculate the solvation free energy. Here, λ = 1.0 corresponds to the state when all the non-bonded interactions of the drug are fully on, and λ = 0.0 corresponds to the state when all the non-bonded interactions of the drug (both intramolecular and extramolecular) are fully off.
	
	
	
	. TO better understand the AWH algorithm, please refer to this two Gromacs tutorials: https://tutorials.gromacs.org/awh-tutorial.html for using harmonic potential to pull along defined reaction coordinate and "https://tutorials.gromacs.org/awh-free-energy-of-solvation.html" for using AWH to compute free energy of solvation)
		
		
		
		
		
		

		 
















