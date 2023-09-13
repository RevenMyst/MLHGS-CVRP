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
		if (nbIter % params.ap.nbIterTraces == 0) population.printState(nbIter, nbIterNonProd);

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
		std::vector<Individual*> offsprings = exhaustiveCrossoverOX(population.getBinaryTournament(),population.getBinaryTournament());
		Individual* mostDiverse = nullptr;
		double maxDiversity = -1.;
		for(Individual* ind : offsprings)
		{
			offspring = *ind;
			/* LOCAL SEARCH */
			localSearch.run(offspring, params.penaltyCapacity, params.penaltyDuration);
			double diversity = population.averageBrokenPairsDistanceAll(offspring);
			if(diversity > maxDiversity)
			{
				maxDiversity = diversity;
				mostDiverse = ind;
			}
		}
		offspring = *mostDiverse;

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
		if (nbIter % params.ap.nbIterTraces == 0) population.printState(nbIter, nbIterNonProd);

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

void Genetic::specificCrossoverOX(Individual* result, const Individual & parent1, const Individual & parent2, int start, int end)
{
	// Frequency table to track the customers which have been already inserted
	std::vector <bool> freqClient = std::vector <bool> (params.nbClients + 1, false);

	// Copy from start to end
	int j = start;
	while (j % params.nbClients != (end + 1) % params.nbClients)
	{
		result->chromT[j % params.nbClients] = parent1.chromT[j % params.nbClients];
		freqClient[result->chromT[j % params.nbClients]] = true;
		j++;
	}

	// Fill the remaining elements in the order given by the second parent
	for (int i = 1; i <= params.nbClients; i++)
	{
		int temp = parent2.chromT[(end + i) % params.nbClients];
		if (freqClient[temp] == false)
		{
			result->chromT[j % params.nbClients] = temp;
			j++;
		}
	}

	// Complete the individual with the Split algorithm
	split.generalSplit(*result, parent1.eval.nbRoutes);
}

std::vector<Individual*> Genetic::exhaustiveCrossoverOX(const Individual & parent1, const Individual & parent2)
{
	std::vector<Individual*> result;
	for (int i = 0; i < params.nbClients; i++)
	{
		for (int j = i + 1; j < params.nbClients; j++)
		{
			Individual* temp = new Individual(params);
			specificCrossoverOX(temp, parent1, parent2, i, j);
			result.push_back(temp);
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

