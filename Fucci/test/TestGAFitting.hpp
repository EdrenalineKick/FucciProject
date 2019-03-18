#ifndef PROJECTS_FUCCI_TEST_TESTGAFITTING_HPP_
#define PROJECTS_FUCCI_TEST_TESTGAFITTING_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>

#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"

#include "AttilaHybrid2006OdeSystem.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "G2TysonNovak2001OdeSystem.hpp"

#include <fstream>
#include "RungeKutta4IvpOdeSolver.hpp"
#include "TysonNovak2001OdeSystem.hpp"
#include "FakePetscSetup.hpp"

//----------------------------------------------------------------------------------------
//
//  define a data structure which will define a chromosome
//
//----------------------------------------------------------------------------------------

int PARAM_SIZE = 4;
double pmut = 0.1;
int popsize = 16;
double pcross = 0.6;

struct chromo_typ
{
    std::vector<double> params;

    double fitness;

    chromo_typ() : params(4, 0), fitness(0.0){};
    chromo_typ(std::vector<double> parms, double ftns) : params(parms), fitness(ftns) {}
};

double Objective(std::vector<double> parms);
std::vector<double> GetRandomParams(int size);
std::vector<double> Roulette(double total_fitness, chromo_typ* Population);
void Mutate(std::vector<double>& params);
void Crossover(std::vector<double>& offspring1, std::vector<double>& offspring2);

class TestGaFitting : public AbstractCellBasedTestSuite
{
public:
    void TestSolvingOdes()
    {
        int ngen = 10;

        chromo_typ Population[16];

        for (int i = 0; i < popsize; i++)
        {
            Population[i].params = GetRandomParams(PARAM_SIZE);
            Population[i].fitness = 0.0;
        }

        //enter the main GA loop
        for (int i = 0; i < ngen; i++)
        {
            //this is used during roulette wheel sampling
            double TotalFitness = 0.0;

            // test and update the fitness of every chromosome in the
            // population
            for (int i = 0; i < popsize; i++)
            {
                Population[i].fitness = Objective(Population[i].params);
                std::cout << "Fitness: " << Population[i].fitness << "\n";
                TotalFitness += Population[i].fitness;
            }

            chromo_typ temp[16];

            int cPop = 0;

            //loop until we have created popsize new chromosomes
            while (cPop < popsize)
            {
                // we are going to create the new population by grabbing members of the old population
                // two at a time via roulette wheel selection.
                std::vector<double> offspring1 = Roulette(TotalFitness, Population);
                std::vector<double> offspring2 = Roulette(TotalFitness, Population);

                //add crossover dependent on the crossover rate
                Crossover(offspring1, offspring2);

                //now mutate dependent on the mutation rate
                Mutate(offspring1);
                Mutate(offspring2);

                //add these offspring to the new population. (assigning zero as their
                //fitness scores)
                temp[cPop++] = chromo_typ(offspring1, 0.0);
                temp[cPop++] = chromo_typ(offspring2, 0.0);

            } //end loop

            for (int i = 0; i < popsize; i++)
            {
                Population[i] = temp[i];
            }
        }

        //One last assessment of fitness to find the best individual.
        for (int i = 0; i < popsize; i++)
        {
            Population[i].fitness = Objective(Population[i].params);
        }

        int index = 0;
        double current_max = 0.0;
        for (int i = 0; i < popsize; i++)
        {
            if (Population[i].fitness > current_max)
            {
                index = i;
                current_max = Population[i].fitness;
                std::cout << "currentMax" << current_max << "\n";
            }
        }

        for (int i = 0; i < PARAM_SIZE; i++)
        {
            std::cout << Population[index].params[i] << "\n";
        }

        //     std::ofstream myfile("Params.txt");
        //     if (myfile.is_open())
        //     {
        //         for (int j = 0; j < popsize; j++)
        //         {
        //             for (int i = 0; i < PARAM_SIZE; i++)
        //             {
        //                 myfile << i << " " << Population[index].params[i] << " ";
        //             }
        //             myfile << Population[index].fitness << "\n";
        //         }
        //         myfile.close();
        //     }
        //     else
        //     {
        //         std::cout << "Unable to open file";
        //     }
    }

    void Mutate(std::vector<double>& params)
    {
        for (unsigned i = 0; i < params.size(); i++)
        {
            if (RandomNumberGenerator::Instance()->ranf() < pmut)
            {
                std::cout << "Pre-Mut: " << params[i] << "\n";
                params[i] = params[i] + RandomNumberGenerator::Instance()->NormalRandomDeviate(0, 0.1);
                std::cout << "Post-Mut: " << params[i] << "\n";

                if (params[i] > 1 || params[i] < 0)
                {
                    params[i] = RandomNumberGenerator::Instance()->ranf();
                }
            }
        }

        return;
    }

    void Crossover(std::vector<double>& offspring1, std::vector<double>& offspring2)
    {
        //dependent on the crossover rate
        if (RandomNumberGenerator::Instance()->ranf() < pcross)
        {
            //create a random crossover point
            int crossover = (int)(RandomNumberGenerator::Instance()->ranf() * PARAM_SIZE);

            std::vector<double> t1;
            std::vector<double> t2;

            for (int i = 0; i < PARAM_SIZE; i++)
            {
                if (i < crossover)
                {
                    t1.push_back(offspring1[i]);
                    t2.push_back(offspring2[i]);
                }
                else
                {
                    t1.push_back(offspring2[i]);
                    t2.push_back(offspring1[i]);
                }
            }
            std::cout << "Cross: " << crossover << "\n"
                      << "offspring1: ";
            for (int i = 0; i < PARAM_SIZE; i++)
            {
                std::cout << offspring1[i] << " ";
            }
            std::cout << "\n"
                      << "offspring2: ";
            for (int i = 0; i < PARAM_SIZE; i++)
            {
                std::cout << offspring2[i] << " ";
            }
            std::cout << "\n"
                      << "t1: ";
            for (int i = 0; i < PARAM_SIZE; i++)
            {
                std::cout << t1[i] << " ";
            }
            std::cout << "\n"
                      << "t2: ";
            for (int i = 0; i < PARAM_SIZE; i++)
            {
                std::cout << t2[i] << " ";
            }
            std::cout << "\n";

            offspring1 = t1;
            offspring2 = t2;
        }
    }

    //--------------------------------Roulette-------------------------------------------
    //
    //  selects a chromosome from the population via roulette wheel selection
    //------------------------------------------------------------------------------------
    std::vector<double> Roulette(double total_fitness, chromo_typ* Population)
    {
        //generate a random number between 0 & total fitness count
        double Slice = (double)(RandomNumberGenerator::Instance()->ranf() * total_fitness);

        //go through the chromosones adding up the fitness so far
        double FitnessSoFar = 0.0;

        for (int i = 0; i < popsize; i++)
        {
            if (Population[i].fitness >= 0)
            {
                FitnessSoFar += Population[i].fitness;
            }
            else
            {
                FitnessSoFar -= Population[i].fitness;
            }

            std::cout << "current fitness: " << FitnessSoFar << "\n";
            //if the fitness so far > random number return the chromo at this point
            if (FitnessSoFar >= Slice)
            {
                std::cout << "Selected Fitness: " << FitnessSoFar << "\n";
                return Population[i].params;
            }
        }
        std::cout << "oops";
        return Population[0].params;
    }

    double Objective(std::vector<double> parms)
    {
        RungeKutta4IvpOdeSolver euler_solver;

        double g1Times[16];
        double stopTimes[16];

        for (unsigned j = 0; j < 16; j++)
        {
            AttilaHybrid2006OdeSystem my_ode;
            std::vector<double> initial_conditions = my_ode.GetInitialConditions();

            double mVRand = parms[1] * RandomNumberGenerator::Instance()->ranf();
            my_ode.SetmV(parms[0] * 0.35 + mVRand);
            double mK4bRand = parms[3] * 2* RandomNumberGenerator::Instance()->ranf();
            my_ode.SetmK4b(parms[2] * 3 + mK4bRand);

            OdeSolution solutions = euler_solver.Solve(&my_ode, initial_conditions, 0, 30, 0.0001, 0.1);
            //OdeSolution solutions = euler_solver.Solve(&my_ode,)

            g1Times[j] = my_ode.GetG1EventTime();
            stopTimes[j] = my_ode.GetStopTime();

            std::cout << "G1 Time: " << g1Times[j] << "\n"
                      << "Stop Time: " << stopTimes[j] << "\n";
        }

        // double comparedG1[7] = { 5.0, 6.0, 4.0, 4.5, 4.8, 5.4, 5.6 };
        // double comparedStop[7] = { 15.3, 16.4, 14.5, 16.7, 15.8, 14.6, 15.1 };

        double comparedG1[19] = { 6.5, 6, 6, 3.5, 4, 4.5, 6, 5, 4.5, 4, 4.5, 7, 6.5, 5, 4.5, 4, 4, 4, 3.5 };
        double comparedStop[19] = { 19, 19, 24, 25, 25, 22, 20.5, 18, 21, 20, 18.5, 23.4, 18.5, 22, 20.5, 16.5, 17.5, 23.5, 22 };

        double g1combined[35][2];
        double stopcombined[35][2];

        int i = 0, j = 0, k = 0;
        int n1 = 16;
        int n2 = 19;

        std::sort(g1Times, g1Times + 16);
        std::sort(stopTimes, stopTimes + 16);
        std::sort(comparedG1, comparedG1 + 19);
        std::sort(comparedStop, comparedStop + 19);

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
                g1combined[k][0] = g1Times[i++];
                g1combined[k][1] = index;
                k++;
            }
            else
            {
                double index = (double)k;
                g1combined[k][0] = comparedG1[j++];
                g1combined[k][1] = -index;
                k++;
            }
        }

        while (i < n1)
        {
            double index = (double)k;
            g1combined[k][0] = g1Times[i++];
            g1combined[k][1] = index;
            k++;
        }
        // Store remaining elements of second array
        while (j < n2)
        {
            double index = (double)k;
            g1combined[k][0] = comparedG1[j++];
            g1combined[k][1] = -index;
            k++;
        }

        double score = 0;
        for (int i = 0; i < 35; i++)
        {
            score += g1combined[i][1];
            if (g1combined[i][0] > 50)
            {
                score += 40;
                std::cout << "OOF"
                          << "\n";
            }
        }

        if (score < 0)
        {
            score = -score;
        }

        i = 0, j = 0, k = 0;
        n1 = 16;
        n2 = 19;

        // Traverse both array
        while (i < n1 && j < n2)
        {
            // Check if current element of first
            // array is smaller than current element
            // of second array. If yes, store first
            // array element and increment first array
            // index. Otherwise do same with second array
            if (stopTimes[i] < comparedStop[j])
            {
                double index = (double)k;
                stopcombined[k][0] = stopTimes[i++];
                stopcombined[k][1] = index;
                k++;
            }
            else
            {
                double index = (double)k;
                stopcombined[k][0] = comparedStop[j++];
                stopcombined[k][1] = -index;
                k++;
            }
        }

        while (i < n1)
        {
            double index = (double)k;
            stopcombined[k][0] = stopTimes[i++];
            stopcombined[k][1] = index;
            k++;
        }
        // Store remaining elements of second array
        while (j < n2)
        {
            double index = (double)k;
            stopcombined[k][0] = comparedStop[j++];
            stopcombined[k][1] = -index;
            k++;
        }

        for (int i = 0; i < 35; i++)
        {
            score += stopcombined[i][1]*4;
            if (stopcombined[i][0] > 50)
            {
                score += 50;
                std::cout << "Ouch" << "\n";
            }
        }

        if (score < 0)
        {
            score = -score;
        }

        if (score == 0)
        {
            score = 100;
        }
        std::cout << score << "\n";
        std::ofstream myfile;

        myfile.open ("scores_genetic_algorithms.csv", std::ios::app);
        myfile << score << "," << parms[0] << "," << parms[1] << "," << parms[2] << "," << parms[3] << ",\n";

        std::ofstream finalfile;
        finalfile.open ("genetic_algorithms_final_score.csv", std::ios::app);
        finalfile << parms[0] * 0.35 << "," << parms[1] << "," << parms[2] * 3 << "," << parms[3] * 2 << ",\n";
        return 1 / score;
    }

    std::vector<double> GetRandomParams(int size)
    {
        std::vector<double> params(size);
        for (int i = 0; i < size; i++)
        {
            params[i] = RandomNumberGenerator::Instance()->ranf();
        }
        return params;
    }
};

#endif /* PROJECTS_FUCCI_TEST_TESTGAFITTING_HPP_ */