/*
 * TESTCOLONYBOUNDARY.hpp
 *
 *  Created on: 2 Nov 2017
 *      Author: edward
 */

#ifndef PROJECTS_FUCCI_TEST_TESTCOLONYBOUNDARY_HPP_
#define PROJECTS_FUCCI_TEST_TESTCOLONYBOUNDARY_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "DirectedMotionForce.hpp"
#include "GammaG1CellCycleModel.hpp"
#include "RandomMotionForce.hpp"
#include "UniformCellCycleModel.hpp"

#include "CellVolumesWriter.hpp"

#include "ColonyLinearSpringForce.hpp"

#include "ColonyIDWriter.hpp"
#include "ColonyIDState.hpp"

#include "DensityWriter.hpp"
#include "LocalNeighbourhoodDensity.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellsGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "ColonyBoundaryDetector.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "OffLatticeSimulation.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "RandomCellKiller.hpp"
#include "SmartPointers.hpp"
#include "StemCellProliferativeType.hpp"
#include "VoronoiDataWriter.hpp"
#include "FakePetscSetup.hpp"

class TestColonyBoundary : public AbstractCellBasedTestSuite
{
public:
    void TestAdhesion()
    {

        //generate mesh
        HoneycombMeshGenerator generator(40, 40, 2);
        MutableMesh<2, 2>* p_mesh = generator.GetMesh();

        //Get location indices of the nodes on the mesh
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        std::vector<unsigned> real_indices;

        //Set seeding probability
        double seeding_probability = 0.016;

        //Iterate over the nodes, random cells are pushed into real_indices according to the seeding probability
        for (unsigned i = 0; i < location_indices.size(); i++)
        {
            unsigned cell_index = location_indices[i];
            if (RandomNumberGenerator::Instance()->ranf() < seeding_probability)
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
            GammaG1CellCycleModel* p_cycle_model = new GammaG1CellCycleModel();
            p_cycle_model->SetDimension(2);

            //To avoid a 'pulsing' behaviour with birth events, we set each cell's initial age to be
            // ~U(-12, 0) in the past, as each cell cycle duration is U(11, 13).
            double birth_time = 20.0 * RandomNumberGenerator::Instance()->ranf();
            p_cycle_model->SetBirthTime(-birth_time);
            p_cycle_model->SetShape(2);
            p_cycle_model->SetScale(5);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_stem_type); //Set the cell to be stem cell type
            p_cell->InitialiseCellCycleModel();

            cells.push_back(p_cell);
        }

        //Link mesh and cells together into a cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, real_indices);

        cell_population.AddPopulationWriter<VoronoiDataWriter>();
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        cell_population.AddPopulationWriter<ColonyBoundaryDetector>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<ColonyIDWriter>();
        cell_population.AddCellWriter<DensityWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("DensityTest2");
        simulator.SetEndTime(60.0);
        simulator.SetSamplingTimestepMultiple(60);

        
        MAKE_PTR(ColonyIDState<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(LocalNeighbourhoodDensity<2>, p_density);
        simulator.AddSimulationModifier(p_density);

        // MAKE_PTR(RandomMotionForce<2>, p_random_force);
        // p_random_force->SetMovementParameter(0.4);
        // simulator.AddForce(p_random_force);

        MAKE_PTR(DirectedMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(0.6);
        simulator.AddForce(p_random_force);

        MAKE_PTR(ColonyLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(2.0);
        p_linear_force->SetCentralRestingSpringLength(0.35);
        p_linear_force->SetMeinekeSpringStiffness(45);
        simulator.AddForce(p_linear_force);

        // MAKE_PTR_ARGS(RandomCellKiller<2>, p_killer, (&cell_population, 0.01));
        // simulator.AddCellKiller(p_killer);

        c_vector<double, 2> point = zero_vector<double>(2);
        c_vector<double, 2> normal = zero_vector<double>(2);

        point(0) = 7.0;
        normal(0) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc1);

        point(0) = 48.0;
        normal(0) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        point(0) = 0.0;
        point(1) = 7.0;
        normal(0) = 0.0;
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc3);

        point(1) = 48.0;
        normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc4);

        simulator.Solve();
    }
};

#endif /* PROJECTS_FUCCI_TEST_TESTCOLONYBOUNDARY_HPP_ */
