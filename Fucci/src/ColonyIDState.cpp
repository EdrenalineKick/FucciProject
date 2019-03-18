/*

Copyright (c) 2005-2017, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "ColonyIDState.hpp"
#include "Exception.hpp"

template <unsigned DIM>
ColonyIDState<DIM>::ColonyIDState()
        : AbstractCellBasedSimulationModifier<DIM>()
{
}

template <unsigned DIM>
ColonyIDState<DIM>::~ColonyIDState()
{
}

template <unsigned DIM>
void ColonyIDState<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template <unsigned DIM>
void ColonyIDState<DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get the location index corresponding to this cell
        unsigned index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

        // MeshBasedCellPopulation* pCellPopulation =
        CellPtr pCell = *cell_iter;

        pCell->GetCellData()->SetItem("Colony ID", index);
    }
    UpdateCellData(rCellPopulation);
}

template <unsigned DIM>
void ColonyIDState<DIM>::UpdateCellData(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    rCellPopulation.Update();

    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get the location index corresponding to this cell
        unsigned index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

        // MeshBasedCellPopulation* pCellPopulation =
        CellPtr pCell = *cell_iter;

        // Get the set of neighbouring location indices
        std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringNodeIndices(index);

        MeshBasedCellPopulation<DIM, DIM>* meshCellPopulation = static_cast<MeshBasedCellPopulation<DIM, DIM>*>(&rCellPopulation);
        unsigned colonyId = pCell->GetCellData()->GetItem("Colony ID");
        bool hasNeighbour = false;
        // Iterate over these neighbours
        for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
             neighbour_iter != neighbour_indices.end();
             ++neighbour_iter)
        {
            unsigned neighbour_index = *neighbour_iter;
            if (!meshCellPopulation->IsGhostNode(neighbour_index))
            {
                CellPtr pNeighbour = rCellPopulation.GetCellUsingLocationIndex(neighbour_index);
     
                pNeighbour->GetCellData()->SetItem("Colony ID", colonyId);

                hasNeighbour = true;
            }
        }

        if (hasNeighbour == false)
        {
            pCell->GetCellData()->SetItem("Colony ID", index);
        }
    }
}

template <unsigned DIM>
void ColonyIDState<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ColonyIDState)