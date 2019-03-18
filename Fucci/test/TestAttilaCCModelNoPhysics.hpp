/*
 * TESTATTILACCMODELNOPHYSICS.hpp
 *
 *  Created on: 2 Nov 2017
 *      Author: edward
 */

#ifndef PROJECTS_FUCCI_TEST_TESTATTILACCMODELNOPHYSICS_HPP_
#define PROJECTS_FUCCI_TEST_TESTATTILACCMODELNOPHYSICS_HPP_

#include <stdlib.h>

#include <cxxtest/TestSuite.h>
#include "Attila2006CellCycleModel.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "DirectedMotionForce.hpp"
#include "G2TysonNovakCellCycleModel.hpp"
#include "RandomMotionForce.hpp"

#include "CellVolumesWriter.hpp"

#include "ColonyLinearSpringForce.hpp"

#include "ColonyIDState.hpp"
#include "ColonyIDWriter.hpp"

#include "DensityWriter.hpp"
#include "LocalNeighbourhoodDensity.hpp"

#include "AggregationForce.hpp"
#include "GammaG1CellCycleModel.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellsGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "ColonyBoundaryDetector.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "GhostNodeMeshBasedPopulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "OffLatticeSimulation.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "RandomCellKiller.hpp"
#include "SmartPointers.hpp"
#include "StemCellProliferativeType.hpp"
#include "VoronoiDataWriter.hpp"
#include "FakePetscSetup.hpp"

class TestAttilaCCModelNoPhysics : public AbstractCellBasedTestSuite
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
        
        //generate mesh
        HoneycombMeshGenerator generator(60, 60, 2);
        MutableMesh<2, 2>* p_mesh = generator.GetMesh();

        //Get location indices of the nodes on the mesh
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        std::vector<unsigned> real_indices;

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

        double seeding_probability = 0.016;

        p_gen->Reseed(randomSeed);
        //Iterate over the nodes, random cells are pushed into real_indices according to the seeding probability
        for (unsigned i = 0; i < location_indices.size(); i++)
        {
            unsigned cell_index = location_indices[i];
            if (p_gen->ranf() < seeding_probability)
            {
                real_indices.push_back(cell_index);
            }
        }

        //Generate Cell Pointer vector and pointer to stem cell type and wild type mutation state
        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        boost::shared_ptr<AbstractCellProperty> p_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();

        //Iterate over real_indices, assigning cell cycle and proliferative type to the cells
        for (unsigned i = 0; i < real_indices.size(); i++)
        {
            //Set cell cycle
            Attila2006CellCycleModel* p_cycle_model = new Attila2006CellCycleModel();
            p_cycle_model->SetDimension(2);

            //To avoid a 'pulsing' behaviour with birth events, we set each cell's initial age to be
            // ~U(-12, 0) in the past, as each cell cycle duration is U(11, 13).
            double birth_time = 16.0 * p_gen->ranf();
            p_cycle_model->SetBirthTime(-birth_time);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_stem_type); //Set the cell to be stem cell type

            cells.push_back(p_cell);
        }

        //Link mesh and cells together into a cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, real_indices, false, 30);

        cell_population.AddCellWriter<CellProliferativePhasesWriter>();

        OffLatticeSimulation<2> simulator(cell_population);


        std::string fileName = "TestAttilaCCModelNoPhysics";
        std::string seedNumber = std::to_string(randomSeed);
        std::string outputString = fileName + seedNumber;
        simulator.SetOutputDirectory(outputString);
        simulator.SetEndTime(60);
        simulator.SetSamplingTimestepMultiple(60);
        // MAKE_PTR(DirectedMotionForce<2>, p_random_force);
        // p_random_force->SetMovementParameter(0.6);
        // simulator.AddForce(p_random_force);


        // MAKE_PTR_ARGS(RandomCellKiller<2>, p_killer, (&cell_population, 0.01));
        // simulator.AddCellKiller(p_killer);

        simulator.Solve();
    }
};

#endif /* PROJECTS_FUCCI_TEST_TESTATTILACCMODELNOPHYSICS_HPP_ */
