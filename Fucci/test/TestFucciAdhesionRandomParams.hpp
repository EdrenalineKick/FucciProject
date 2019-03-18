/*
 * TestFucciAdhesionRandomParams.hpp
 *
 *  Created on: 2 Nov 2017
 *      Author: edward
 */

#ifndef PROJECTS_FUCCI_TEST_TESTFUCCIADHESIONRANDOMPARAMS_HPP_
#define PROJECTS_FUCCI_TEST_TESTFUCCIADHESIONRANDOMPARAMS_HPP_

#include <cxxtest/TestSuite.h>
#include "GammaG1CellCycleModel.hpp"
#include "RandomMotionForce.hpp"
#include "CheckpointArchiveTypes.hpp"
#include <time.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "StemCellProliferativeType.hpp"
#include "OffLatticeSimulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "RandomCellKiller.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "FakePetscSetup.hpp"

class TestFucciAdhesionRandomParams : public AbstractCellBasedTestSuite
{
public:
	void TestAdhesion()
	{
		//Take Command line arguments to be utilised as the random seed number
		CommandLineArguments* p_args = CommandLineArguments::Instance();
        unsigned argc = *(p_args->p_argc); // has the number of arguments, and
        char **argv = *(p_args->p_argv); // is a char** of them.
        std::cout << "#" << argc-1 << " arguments supplied.\n" << std::flush;

        if (argc != 2)
        {
            std::cerr << "TestAdhesion::Please input one argument\n"
                      << std::flush;
            return;
        }

        unsigned randomSeed = atof(argv[1]);

        std::cout << "Random Seed: " << randomSeed << "\n" << std::flush;

        //generate mesh
        HoneycombMeshGenerator generator(60, 60, 2);
        MutableMesh<2, 2>* p_mesh = generator.GetMesh();

        //Get location indices of the nodes on the mesh
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        std::vector<unsigned> real_indices;

        //Set seeding probability
        double seeding_probability = 0.016;

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

        //Iterate over the nodes, random cells are pushed into real_indices according to the seeding probability
        for (unsigned i = 0; i < location_indices.size(); i++)
        {
            unsigned cell_index = location_indices[i];
            if (p_gen->ranf() < seeding_probability)
            {
                real_indices.push_back(cell_index);
            }
		}

//		for (unsigned i = 0 ; i < location_indices.size(); i++)
//		{
//			unsigned cell_index = location_indices[i];
//			if (i == 1 || i == 7 || i == 16)
//			{
//				real_indices.push_back(cell_index);
//			}
//		}

		//Generate Cell Pointer vector and pointer to stem cell type and wild type mutation state
		std::vector<CellPtr> cells;
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		boost::shared_ptr<AbstractCellProperty> p_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();

		//Iterate over real_indices, assigning cell cycle and proliferative type to the cells
		for (unsigned i = 0; i<real_indices.size(); i++)
		{
			//Set cell cycle
			GammaG1CellCycleModel* p_cycle_model = new GammaG1CellCycleModel();
			p_cycle_model->SetDimension(2);

			//To avoid a 'pulsing' behaviour with birth events, we set each cell's initial age to be
			// ~U(-12, 0) in the past, as each cell cycle duration is U(11, 13).
			double birth_time = 12.0*RandomNumberGenerator::Instance()->ranf();
			p_cycle_model->SetBirthTime(-birth_time);
			p_cycle_model->SetShape(12);
			p_cycle_model->SetScale(0.4);

			CellPtr p_cell(new Cell(p_state, p_cycle_model));
			p_cell->SetCellProliferativeType(p_stem_type); //Set the cell to be stem cell type
			p_cell->InitialiseCellCycleModel();

			cells.push_back(p_cell);
		}

		//Link mesh and cells together into a cell population
		MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, real_indices);

		cell_population.AddPopulationWriter<VoronoiDataWriter>();
		cell_population.AddCellWriter<CellProliferativePhasesWriter>();




		p_gen->Reseed(randomSeed);
		

		double movementRand = p_gen->ranf();
		double deathRand = p_gen->ranf();
		double springRand = p_gen->ranf();

		double movementParameter = movementRand*1.5;
		double deathParameter = deathRand*0.05;
		double springParameter = springRand*60;		

		OffLatticeSimulation<2> simulator(cell_population);

		MAKE_PTR(RandomMotionForce<2>, p_random_force);
		p_random_force->SetMovementParameter(movementParameter);
		simulator.AddForce(p_random_force);
		
		MAKE_PTR_ARGS(RandomCellKiller<2>,p_killer , (&cell_population, deathParameter));
		simulator.AddCellKiller(p_killer);

		MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
		p_linear_force->SetCutOffLength(2.0);
		p_linear_force->SetMeinekeSpringStiffness(springParameter);
		simulator.AddForce(p_linear_force);

		std::string movementString = std::to_string(movementParameter);
		std::string deathString = std::to_string(deathParameter);
		std::string springString = std::to_string(springParameter);
		std::string seedString = std::to_string(randomSeed);

		std::string outputString = seedString + "Movement" + movementString + "Death" + deathString + "Spring" + springString;


		simulator.SetOutputDirectory(outputString);
		simulator.SetEndTime(10);
		simulator.SetSamplingTimestepMultiple(60);

		c_vector<double,2> point = zero_vector<double>(2);
		c_vector<double,2> normal = zero_vector<double>(2);

		point(0) = 5.0;
		normal(0) = -1.0;
		MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
		simulator.AddCellPopulationBoundaryCondition(p_bc1);

		point(0) = 50.0;
		normal(0) = 1.0;
		MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal));
		simulator.AddCellPopulationBoundaryCondition(p_bc2);

		point(0) = 0.0;
		point(1) = 5.0;
		normal(0) = 0.0;
		normal(1) = -1.0;
		MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal));
		simulator.AddCellPopulationBoundaryCondition(p_bc3);

		point(1) = 50.0;
		normal(1) = 1.0;
		MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal));
		simulator.AddCellPopulationBoundaryCondition(p_bc4);


		simulator.Solve();
	}
};

#endif /* PROJECTS_FUCCI_TEST_TESTFUCCIADHESIONRANDOMPARAMS_HPP_ */
