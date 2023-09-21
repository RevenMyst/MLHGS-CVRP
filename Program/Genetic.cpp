#include "Genetic.h"

void Genetic::run()
{	
	/* INITIAL POPULATION */
	population.generatePopulation();

	int nbIter;
	int nbIterNonProd = 1;
	if (params.verbose) std::cout << "----- STARTING GENETIC ALGORITHM" << std::endl;
	for (nbIter = 0 ; nbIterNonProd <= params.ap.nbIter && (params.ap.timeLimit == 0 || (double)(clock()-params.startTime)/(double)CLOCKS_PER_SEC < params.ap.timeLimit) && (params.ap.nbIterTotal == -1 || params.ap.nbIterTotal > nbIter); nbIter++)
	{	
		/* SELECTION AND CROSSOVER */
		crossoverOX(offspring, population.getBinaryTournament(),population.getBinaryTournament());

		/* LOCAL SEARCH */
		localSearch.run(offspring, params.penaltyCapacity, params.penaltyDuration);
		bool isNewBest = population.addIndividual(offspring,true, nbIter);
		if (!offspring.eval.isFeasible && params.ran()%2 == 0) // Repair half of the solutions in case of infeasibility
		{
			localSearch.run(offspring, params.penaltyCapacity*10., params.penaltyDuration*10.);
			if (offspring.eval.isFeasible) isNewBest = (population.addIndividual(offspring,false, nbIter) || isNewBest);
		}

		/* TRACKING THE NUMBER OF ITERATIONS SINCE LAST SOLUTION IMPROVEMENT */
		if (isNewBest) nbIterNonProd = 1;
		else nbIterNonProd ++ ;

		/* DIVERSIFICATION, PENALTY MANAGEMENT AND TRACES */
		if (nbIter % params.ap.nbIterPenaltyManagement == 0) population.managePenalties();
		if (nbIter % params.ap.nbIterTraces == 0){ 
			population.printState(nbIter, nbIterNonProd);
			population.fullSearchProgress.push_back({ clock() - params.startTime, population.getBestFound()->eval.penalizedCost, nbIter, population.getDiversity(population.getFeasibleSubPop()), population.getDiversity(population.getInfeasibleSubPop())});
			}

		/* FOR TESTS INVOLVING SUCCESSIVE RUNS UNTIL A TIME LIMIT: WE RESET THE ALGORITHM/POPULATION EACH TIME maxIterNonProd IS ATTAINED*/
		if (params.ap.timeLimit != 0 && nbIterNonProd == params.ap.nbIter)
		{
			population.restart();
			nbIterNonProd = 1;
		}
	}
	if (params.verbose) std::cout << "----- GENETIC ALGORITHM FINISHED AFTER " << nbIter << " ITERATIONS. TIME SPENT: " << (double)(clock() - params.startTime) / (double)CLOCKS_PER_SEC << std::endl;
}

void Genetic::run2()
{	
	/* INITIAL POPULATION */
	population.generatePopulation();

	int nbIter;
	int nbIterNonProd = 1;
	if (params.verbose) std::cout << "----- STARTING GENETIC ALGORITHM" << std::endl;
	for (nbIter = 0 ; nbIterNonProd <= params.ap.nbIter && (params.ap.timeLimit == 0 || (double)(clock()-params.startTime)/(double)CLOCKS_PER_SEC < params.ap.timeLimit) && (params.ap.nbIterTotal == -1 || params.ap.nbIterTotal > nbIter); nbIter++)
	{	
		/* SELECTION AND CROSSOVER */
		std::vector<Individual> offsprings = exhaustiveCrossoverOX(population.getRanked(),population.getRanked());
		double best_score = -1;
		std::vector<double> diversities = std::vector<double>(offsprings.size(),0.0);
		std::vector<double> costs = std::vector<double>(offsprings.size(),0.0);
		double min_cost = -1;
		double max_cost = -1;
		//std::cout<<"start"<<std::endl;
		for(int i = 0; i < offsprings.size(); i++)
		{
			offspring = offsprings[i];
			localSearch.run(offspring, params.penaltyCapacity, params.penaltyDuration);
			offsprings[i] = offspring;
			diversities[i] = population.averageBrokenPairsDistanceAll(offsprings[i]);
			costs[i] = offspring.eval.penalizedCost;
			if(min_cost > offspring.eval.penalizedCost || min_cost == -1) min_cost = offspring.eval.penalizedCost;
			if(max_cost < offspring.eval.penalizedCost || max_cost == -1) max_cost = offspring.eval.penalizedCost;
				
		}
		for(int i = 0; i < offsprings.size(); i++)
		{
			// set cost between 0 and 1 (1 is the best)
			costs[i] = (max_cost - costs[i]) / (max_cost - min_cost);
			if( costs[i]*diversities[i] > best_score || best_score == -1)
			{
				best_score = costs[i]*diversities[i];
				offspring = offsprings[i];
			}
		}

		//std::cout << offspring.eval.penalizedCost << " " << min_cost << " " << best_score << std::endl;


		bool isNewBest = population.addIndividual(offspring,true, nbIter);
		if (!offspring.eval.isFeasible && params.ran()%2 == 0) // Repair half of the solutions in case of infeasibility
		{
			localSearch.run(offspring, params.penaltyCapacity*10., params.penaltyDuration*10.);
			if (offspring.eval.isFeasible) isNewBest = (population.addIndividual(offspring,false, nbIter) || isNewBest);
		}

		/* TRACKING THE NUMBER OF ITERATIONS SINCE LAST SOLUTION IMPROVEMENT */
		if (isNewBest) nbIterNonProd = 1;
		else nbIterNonProd ++ ;

		/* DIVERSIFICATION, PENALTY MANAGEMENT AND TRACES */
		if (nbIter % params.ap.nbIterPenaltyManagement == 0) population.managePenalties();
		if (nbIter % params.ap.nbIterTraces == 0) {
			population.printState(nbIter, nbIterNonProd);
			population.fullSearchProgress.push_back({ clock() - params.startTime, population.getBestFound()->eval.penalizedCost, nbIter, population.getDiversity(population.getFeasibleSubPop()), population.getDiversity(population.getInfeasibleSubPop())});
		}

		/* FOR TESTS INVOLVING SUCCESSIVE RUNS UNTIL A TIME LIMIT: WE RESET THE ALGORITHM/POPULATION EACH TIME maxIterNonProd IS ATTAINED*/
		if (params.ap.timeLimit != 0 && nbIterNonProd == params.ap.nbIter)
		{
			population.restart();
			nbIterNonProd = 1;
		}
	}
	if (params.verbose) std::cout << "----- GENETIC ALGORITHM FINISHED AFTER " << nbIter << " ITERATIONS. TIME SPENT: " << (double)(clock() - params.startTime) / (double)CLOCKS_PER_SEC << std::endl;
}

void Genetic::run3()
{	
	/* INITIAL POPULATION */
	population.generatePopulation();

	int nbIter;
	int nbIterNonProd = 1;
	if (params.verbose) std::cout << "----- STARTING GENETIC ALGORITHM" << std::endl;
	for (nbIter = 0 ; nbIterNonProd <= params.ap.nbIter && (params.ap.timeLimit == 0 || (double)(clock()-params.startTime)/(double)CLOCKS_PER_SEC < params.ap.timeLimit) && (params.ap.nbIterTotal == -1 || params.ap.nbIterTotal > nbIter); nbIter++)
	{	
		/* SELECTION AND CROSSOVER */
		std::vector<Individual> offsprings = exhaustiveCrossoverOX(population.getRanked(),population.getRanked());
		double best_score = -1;
		bool select_div = false;
		double current = 0;
		Individual & bestIndiv = offspring ;
		if(params.ran()%2==0){
			select_div = true;
		}
		for(int i = 0; i < offsprings.size(); i++)
		{
			offspring = offsprings[i];
			localSearch.run(offspring, params.penaltyCapacity, params.penaltyDuration);
			offsprings[i] = offspring;
			if(select_div){
				current =  1-population.averageBrokenPairsDistanceAll(offsprings[i]);
			}else{
				current = offspring.eval.penalizedCost;
			}
			if(best_score == -1 || current < best_score) {
				best_score = current;
				bestIndiv = offspring;
			};
				
		}

		offspring = bestIndiv;
		//std::cout << offspring.eval.penalizedCost << " " << min_cost << " " << best_score << std::endl;


		best_cuts.push_back({parent1.chromT, parent2.chromT, offspring.cut1, offspring.cut2});

		bool isNewBest = population.addIndividual(offspring,true, nbIter);
		if (!offspring.eval.isFeasible && params.ran()%2 == 0) // Repair half of the solutions in case of infeasibility
		{
			localSearch.run(offspring, params.penaltyCapacity*10., params.penaltyDuration*10.);
			if (offspring.eval.isFeasible) isNewBest = (population.addIndividual(offspring,false, nbIter) || isNewBest);
		}

		/* TRACKING THE NUMBER OF ITERATIONS SINCE LAST SOLUTION IMPROVEMENT */
		if (isNewBest) nbIterNonProd = 1;
		else nbIterNonProd ++ ;

		/* DIVERSIFICATION, PENALTY MANAGEMENT AND TRACES */
		if (nbIter % params.ap.nbIterPenaltyManagement == 0) population.managePenalties();
		if (nbIter % params.ap.nbIterTraces == 0) {
			population.printState(nbIter, nbIterNonProd);
			population.fullSearchProgress.push_back({ clock() - params.startTime, population.getBestFound()->eval.penalizedCost, nbIter, population.getDiversity(population.getFeasibleSubPop()), population.getDiversity(population.getInfeasibleSubPop())});
		}

		/* FOR TESTS INVOLVING SUCCESSIVE RUNS UNTIL A TIME LIMIT: WE RESET THE ALGORITHM/POPULATION EACH TIME maxIterNonProd IS ATTAINED*/
		if (params.ap.timeLimit != 0 && nbIterNonProd == params.ap.nbIter)
		{
			population.restart();
			nbIterNonProd = 1;
		}
	}
	if (params.verbose) std::cout << "----- GENETIC ALGORITHM FINISHED AFTER " << nbIter << " ITERATIONS. TIME SPENT: " << (double)(clock() - params.startTime) / (double)CLOCKS_PER_SEC << std::endl;
}

void Genetic::crossoverOX(Individual & result, const Individual & parent1, const Individual & parent2)
{
	// Frequency table to track the customers which have been already inserted
	std::vector <bool> freqClient = std::vector <bool> (params.nbClients + 1, false);

	// Picking the beginning and end of the crossover zone
	std::uniform_int_distribution<> distr(0, params.nbClients-1);
	int start = distr(params.ran);
	int end = distr(params.ran);

	// Avoid that start and end coincide by accident
	while (end == start) end = distr(params.ran);

	// Copy from start to end
	int j = start;
	while (j % params.nbClients != (end + 1) % params.nbClients)
	{
		result.chromT[j % params.nbClients] = parent1.chromT[j % params.nbClients];
		freqClient[result.chromT[j % params.nbClients]] = true;
		j++;
	}

	// Fill the remaining elements in the order given by the second parent
	for (int i = 1; i <= params.nbClients; i++)
	{
		int temp = parent2.chromT[(end + i) % params.nbClients];
		if (freqClient[temp] == false)
		{
			result.chromT[j % params.nbClients] = temp;
			j++;
		}
	}

	// Complete the individual with the Split algorithm
	split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::specificCrossoverOX(Individual& result, const Individual & parent1, const Individual & parent2, int start, int end)
{
	// Frequency table to track the customers which have been already inserted
	std::vector <bool> freqClient = std::vector <bool> (params.nbClients + 1, false);

	// Copy from start to end
	int j = start;
	while (j % params.nbClients != (end + 1) % params.nbClients)
	{
		result.chromT[j % params.nbClients] = parent1.chromT[j % params.nbClients];
		freqClient[result.chromT[j % params.nbClients]] = true;
		j++;
	}

	// Fill the remaining elements in the order given by the second parent
	for (int i = 1; i <= params.nbClients; i++)
	{
		int temp = parent2.chromT[(end + i) % params.nbClients];
		if (freqClient[temp] == false)
		{
			result.chromT[j % params.nbClients] = temp;
			j++;
		}
	}

	// Complete the individual with the Split algorithm
<<<<<<< HEAD
	split.generalSplit(result, parent1.eval.nbRoutes);
=======
	split.generalSplit(*result, parent1.eval.nbRoutes);
	result->cut1 = start;
	result->cut2 = end;
>>>>>>> bc09749d17410f94438bebf58afe163a9754ff59
}

std::vector<Individual> Genetic::exhaustiveCrossoverOX(const Individual & parent1, const Individual & parent2)
{
	std::vector<Individual> result;
	int count = 0;
	for (int i = 0; i < params.nbClients; i++)
	{
		for (int j = i + 1; j < params.nbClients; j++)
		{
			result.push_back(Individual(params));
			specificCrossoverOX(result[count++], parent1, parent2, i, j);
		}
	}
	return result;
}

Genetic::Genetic(Params & params) : 
	params(params), 
	split(params),
	localSearch(params),
	population(params,this->split,this->localSearch),
	offspring(params){}

