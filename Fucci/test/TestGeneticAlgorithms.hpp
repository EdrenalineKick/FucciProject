#ifndef PROJECTS_FUCCI_TEST_TESTGENETICALGORITHMS_HPP_
#define PROJECTS_FUCCI_TEST_TESTGENETICALGORITHMS_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"

#include "AttilaHybrid2006OdeSystem.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "G2TysonNovak2001OdeSystem.hpp"

#include <ga.h>

#include <GARealGenome.c>
#include <GARealGenome.h>
#include <GASimpleGA.c>
#include <GASimpleGA.h>
#include "RungeKutta4IvpOdeSolver.hpp"
#include "TysonNovak2001OdeSystem.hpp"
#include "FakePetscSetup.hpp"

float Objective1(GAGenome&);

class TestGeneticAlgorithms : public AbstractCellBasedTestSuite
{
public:
    void TestSolvingOdes()
    {
        // int popsize = 30;
        // int ngen = 400;
        // double pmut = 0.01;
        // double pcross = 0.6;

        GARealAlleleSet alleles(0, 1);

        int param = 4;
        GARealGenome genome(param, alleles, Objective1);

        GASimpleGA ga(genome);

        // GAParameterList params;
        // GASimpleGA::registerDefaultParameters(params);
        // params.set(gaNnGenerations, ngen);
        // params.set(gaNpopulationSize, popsize);
        // params.set(gaNpMutation,pmut);
        // params.set(gaNpCrossover,pcross);
        // params.set(gaNscoreFrequency, 10);
        // params.set(gaNflushFrequency, 50);
        // params.set(gaNselectScores, (int)GAStatistics::AllScores);

        std::cout << "evolving...";
        while (!ga.done())
        {
            ga.step();
            if (ga.generation() % 50 == 0)
            {
                std::cout << ".";
            }
        }
        std::cout << "\n\n";

        std::cout << "the ga generated:\n"
                  << ga.statistics().bestIndividual() << "\n";
    }

    static float Objective1(GAGenome& g)
    {
        GARealGenome& genome = (GARealGenome&)g;
        AttilaHybrid2006OdeSystem my_ode;

        RungeKutta4IvpOdeSolver euler_solver;
        std::vector<double> initial_conditions = my_ode.GetInitialConditions();

        my_ode.SetmV(genome.gene(0) * 0.35);
        my_ode.SetmVRand(genome.gene(1) * 0.1);
        my_ode.SetmK4b(genome.gene(2) * 3);
        my_ode.SetmK4bRand(genome.gene(3) * 0.5);

        double g1Times[30];
        double stopTimes[30];

        for (unsigned j = 0; j < 30; j++)
        {
            OdeSolution solutions = euler_solver.Solve(&my_ode, initial_conditions, 0, 30, 0.0001, 0.1);
            //OdeSolution solutions = euler_solver.Solve(&my_ode,)

            g1Times[j] = my_ode.GetG1EventTime();
            stopTimes[j] = my_ode.GetStopTime();
        }

        double comparedG1[7] = { 5.0, 6.0, 4.0, 4.5, 4.8, 5.4, 5.6 };
        double comparedStop[7] = { 15.3, 16.4, 14.5, 16.7, 15.8, 14.6, 15.1 };

        double g1combined[37][2];

        int i = 0, j = 0, k = 0;
        int n1 = 30;
        int n2 = 7;

        std::sort(g1Times, g1Times + 30);
        std::sort(stopTimes, stopTimes + 30);
        std::sort(comparedG1, comparedG1 + 7);
        std::sort(comparedStop, comparedStop + 7);

        // Traverse both array
        while (i < n1 && j < n2)
        {
            // Check if current element of first
            // array is smaller than current element
            // of second array. If yes, store first
            // array element and increment first array
            // index. Otherwise do same with second array
            if (g1Times[i] < comparedG1[j])
            {
                double index = (double)k;
                g1combined[k++][0] = g1Times[i++];
                g1combined[k][1] = index;
            }
            else
            {
                double index = (double)k;
                g1combined[k++][0] = comparedG1[j++];
                g1combined[k][1] = -index;
            }
        }

        while (i < n1)
        {
            double index = (double)k;
            g1combined[k++][0] = g1Times[i++];
            g1combined[k][1] = index;
        }

        // Store remaining elements of second array
        while (j < n2)
        {
            double index = (double)k;
            g1combined[k++][0] = comparedG1[j++];
            g1combined[k][1] = -index;
        }

        double score = g1combined[1][0];

        return score;
    }
};

#endif /* PROJECTS_FUCCI_TEST_TESTGENETICALGORITHMS_HPP_ */