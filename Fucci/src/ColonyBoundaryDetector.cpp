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

#include "ColonyBoundaryDetector.hpp"

#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"

#include "CellLabel.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ColonyBoundaryDetector<ELEMENT_DIM, SPACE_DIM>::ColonyBoundaryDetector()
        : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("colonyboundary.dat")
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ColonyBoundaryDetector<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
         cell_iter != pCellPopulation->End();
         ++cell_iter)
    {
        // Get the location index corresponding to this cell
        unsigned index = pCellPopulation->GetLocationIndexUsingCell(*cell_iter);
        CellPtr pCell = *cell_iter;
        // Get the set of neighbouring location indices
        std::set<unsigned> neighbour_indices = pCellPopulation->GetNeighbouringNodeIndices(index);

        boost::shared_ptr<AbstractCellProperty> p_label = pCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<CellLabel>();

        pCell->RemoveCellProperty<CellLabel>();
        // Iterate over these neighbours
        for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
             neighbour_iter != neighbour_indices.end();
             ++neighbour_iter)
        {
            unsigned neighbour_index = *neighbour_iter;

            if (pCellPopulation->IsGhostNode(neighbour_index))
            {
                pCell->AddCellProperty(p_label);
            }
        }

        // Store whether this cell is labelled
        bool cell_is_labelled = pCell->template HasCellProperty<CellLabel>();
        *this->mpOutStream << cell_is_labelled << " ";
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ColonyBoundaryDetector<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ColonyBoundaryDetector<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ColonyBoundaryDetector<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ColonyBoundaryDetector<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{}

    // Explicit instantiation
    template class ColonyBoundaryDetector<1, 1>;
    template class ColonyBoundaryDetector<1, 2>;
    template class ColonyBoundaryDetector<2, 2>;
    template class ColonyBoundaryDetector<1, 3>;
    template class ColonyBoundaryDetector<2, 3>;
    template class ColonyBoundaryDetector<3, 3>;

#include "SerializationExportWrapperForCpp.hpp"
    // Declare identifier for the serializer
    EXPORT_TEMPLATE_CLASS_ALL_DIMS(ColonyBoundaryDetector)
