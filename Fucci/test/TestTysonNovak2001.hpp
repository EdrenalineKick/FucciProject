#ifndef PROJECTS_FUCCI_TEST_TESTTYSONNOVAK2001_HPP_
#define PROJECTS_FUCCI_TEST_TESTTYSONNOVAK2001_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include "G2TysonNovak2001OdeSystem.hpp"

#include "FakePetscSetup.hpp"

class TestTysonNovakOde: public CxxTest::TestSuite
{
public:
    void TestOde()
    {
        G2TysonNovak2001OdeSystem my_ode;
        EulerIvpOdeSolver euler_solver;
        //std::vector<double> initial_condition;
        //initial_condition.push_back(1.0);
        std::vector<double> initial_condition = my_ode.GetInitialConditions();

        OdeSolution solutions = euler_solver.Solve(&my_ode, initial_condition, 0, 2, 0.0001, 0.1);
        for (unsigned i=0; i<solutions.rGetTimes().size(); i++)
        {
            std::cout << solutions.rGetTimes()[i] << " " << solutions.rGetSolutions()[i][0] << " " << solutions.rGetSolutions()[i][1] 
            << " " << solutions.rGetSolutions()[i][2] << " " << solutions.rGetSolutions()[i][3]<< "\n";
        }

        solutions.WriteToFile("TestG2ODE", "my_ode_solution", "sec");
        
    }
};


#endif /* PROJECTS_FUCCI_TEST_TESTTYSONNOVAK2001_HPP_ */
