#ifndef PROJECTS_FUCCI_TEST_TESTMODELG2PHASECELLCYCLE_HPP_
#define PROJECTS_FUCCI_TEST_TESTMODELG2PHASECELLCYCLE_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellsGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "G2TysonNovakCellCycleModel.hpp"
#include "FucciProportionalCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "VoronoiDataWriter.hpp"
#include "FakePetscSetup.hpp"

class TestG2PhaseCellCycle : public AbstractCellBasedTestSuite
{
public:
    void TestG2Phase()
    {
        // HoneycombMeshGenerator generator(2, 2);    // Parameters are: cells across, cells up
        // MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // std::vector<CellPtr> cells;
        // CellsGenerator<G2TysonNovakCellCycleModel, 2> cells_generator;
        // cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        // MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // cell_population.SetWriteVtkAsPoints(false);
        // cell_population.AddPopulationWriter<VoronoiDataWriter>();

        // OffLatticeSimulation<2> simulator(cell_population);
        // simulator.SetOutputDirectory("G2Testy");
        // simulator.SetEndTime(5.0);

        // simulator.SetSamplingTimestepMultiple(2);

        // simulator.Solve();

        //generate mesh
        HoneycombMeshGenerator generator(2, 2, 4);
        MutableMesh<2, 2>* p_mesh = generator.GetMesh();

        //Get location indices of the nodes on the mesh
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        std::vector<unsigned> real_indices;

        for (unsigned i = 0; i < location_indices.size(); i++)
        {
            unsigned cell_index = location_indices[i];
            if (i == 1)
            {
                real_indices.push_back(cell_index);
            }
        }

		std::vector<CellPtr> cells;

        //Iterate over real_indices, assigning cell cycle and proliferative type to the cells
        for (unsigned i = 0; i < real_indices.size(); i++)
        {
            G2TysonNovakCellCycleModel* p_cell_cycle_model = new G2TysonNovakCellCycleModel;
            p_cell_cycle_model->SetDimension(2);

            boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

            double birth_time = -0.5 * RandomNumberGenerator::Instance()->ranf();

            if (p_cell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
            {
                birth_time = -0.5 * RandomNumberGenerator::Instance()->ranf();
            }

            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        //Link mesh and cells together into a cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, real_indices);

        cell_population.AddPopulationWriter<VoronoiDataWriter>();
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();

        OffLatticeSimulation<2> simulator(cell_population);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(2.0);
        p_linear_force->SetMeinekeSpringStiffness(45);
        simulator.AddForce(p_linear_force);

        simulator.SetOutputDirectory("G2Print");
        simulator.SetEndTime(15);
        simulator.SetSamplingTimestepMultiple(20);

        simulator.Solve();
    }
};

#endif /* PROJECTS_FUCCI_TEST_TESTMODELSCENARIOSODES_HPP_ */