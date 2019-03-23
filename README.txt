Code to obtain the results from
"Non-monotonic auto-regulation in single gene circuits"
Lana Descheemaeker, Eveline Peeters and Sophie de Buyl

Fig 1: Conceptual representation of the two hypotheses
	- File: figures_concepts.py
	- Function: concept_figure(project=Project.PAPER)

Fig 2: Theory, dynamics and implementation of non-monotonic response curves
	- File: figures_concepts.py
	- Function: bistability_oscillation_timeseries_figure(project=Project.PAPER)

Fig 3: Classification of (non-)monotonic curves
	- File: figures_concepts.py
	- Function: nonmonotonicities()

Fig 4: Oscillatory solutions for the different toy models
	1. Random scan for oscillatory solutions
		- File: scan_parameter_space.py
		- Function: random_scan(file_sgs, sgs, system, concept=Concept.OSCILLATION)
						with 	file_sgs							sgs 	system
								'oscillation/MDS/variables.csv'		MDS 	System.RANDOM
								'oscillation/2DS/variables.csv'		DDS 	System.RANDOM
								'oscillation/3DS/variables.csv'		DDDS 	System.RANDOM
								'oscillation/SsLrpB/variables.csv'	DDDS 	System.SSLRPB
	2. Perform timeseries on oscillators to determine period and amplitudes, add data to files
		- File: scan_parameter_space.py
		- Function: add_oscillation_data(file_sgs, sgs)
						with	file_sgs							sgs 
								'oscillation/MDS/variables.csv'		MDS 
								'oscillation/2DS/variables.csv'		DDS 
								'oscillation/3DS/variables.csv'		DDDS 
								'oscillation/SsLrpB/variables.csv'	DDDS
	3. Do bifurcation analysis
		- File: bifurcation.py
		- Function: bifurcation_analysis(sgs, system, concept=Concept.OSCILLATION)
						with 	sgs 	system
								MDS 	System.RANDOM
								DDS 	System.RANDOM
								DDDS 	System.RANDOM
								DDDS 	System.SSLRPB
	4. Plot mean value and mean width of the parameters 
		- File: figures_oscillation.py
		- Function: parameter_ranges_figure(concept=Concept.OSCILLATION)

Fig 5: Amplitudes of oscillation in the monomer (m) - dimer (d) plane for the diffent toy models
	1. Random scan for oscillatory solutions 
			-> See Fig 4. 1
	2. Perform timeseries on oscillators to determine period and amplitudes, add data to files
 			-> See Fig 4. 2
 	3. Plot the amplitudes
 		- File: figures_oscillation.py
 		- Function: amplitudes_figure()

Fig 6: Distribution of induction times for the different toy models
	1. Systematic scan over parameter space to look for bistable solutions 
		and calculate the induction times
		- File: scan_parameter_space.py
		- Function: bistability_scan(sgs, Hs=[500,600])
						with sgs = MDS, DDS and DDDS
					bistability_scan_SsLrpB()
	2. Plot the distribution of the induction times
		- File: figures_bistability.py
		- Function: induction_time_histograms_figure()

Fig 7: Schematic of the Ss-LrpB system

Fig 8: Stochastic time series.
	1. Perform some stochastic timeseries
		- File: plot_SGS.py
		- Set the sgs, give the filename and index to find the parameters
		- Determine the DNAcopynumber
		- Set stochsim to True
		- Add savetxt comment after stochsim to save stochastic time lapse (comment on line 148)
		- Run code
	2. Plot the stochastic timetraces
		- File: figure_stochastic_timelapses.py
		- Set paths to stochastic timeseries
		- Run code

Fig 9: Bistable region with induction times in the f123-f13 plane
	1. Systematic scan over parameter space to look for bistable solutions 
		and calculate the induction times
			-> See Fig 6. 1
	2. Look for random oscillatory solutions
			-> See Fig 4.1 in the case of 	file_sgs							sgs 	system
											'oscillation/SsLrpB/variables.csv'	DDDS 	System.SSLRPB
	3. Calculate percentage of compatible solutions
		- File: figure_SsLrpB_analysis.py
		- Function: percentage_SsLrpB_compatible_solutions(fname='bistability/SsLrpB/compatible_solutions.csv')
	3. Plot all data in one figure
		- File: figure_SsLrpB_analysis.py
		- Function: plot_bistable_compatible_regions('bistability/SsLrpB/grid/', 'oscillation/SsLrpB/variables.csv')

Fig 10: Bursting behavior of the 3DS systems
	

Table 2: Summary of the search for oscillations and the scan for bistability
	1. Random scan for oscillatory solutions
		-> See Fig 4. 1
	2. Systematic scan over parameter space to look for bistable solutions 
		-> See Fig 6. 1
	3. Count oscillators for every monotonicity type
		- File: scan_parameter_space.py
		- Function: count_solutions_oscillations(file)
			with file = 'oscillation/MDS/variables.csv'
						'oscillation/2DS/variables.csv'
						'oscillation/3DS/variables.csv'
						'oscillation/SsLrpB/variables.csv'
	4. Count bistable systems for every monotonicity type
		- File: scan_parameter_space.py
		- Function: count_solutions_bistability()


