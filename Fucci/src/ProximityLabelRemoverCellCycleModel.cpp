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

#include "ProximityLabelRemoverCellCycleModel.hpp"

#include "CellLabel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"

ProximityLabelRemoverCellCycleModel::ProximityLabelRemoverCellCycleModel()
    : AbstractSimplePhaseBasedCellCycleModel(),
      mCompressedVolumeFraction(DOUBLE_UNSET),
      mEquilibriumVolume(DOUBLE_UNSET)
{
}

ProximityLabelRemoverCellCycleModel::ProximityLabelRemoverCellCycleModel(const ProximityLabelRemoverCellCycleModel& rModel)
    : AbstractSimplePhaseBasedCellCycleModel(rModel),
      mCompressedVolumeFraction(rModel.mCompressedVolumeFraction),
      mEquilibriumVolume(rModel.mEquilibriumVolume)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractPhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mG1Duration is (re)set as soon as InitialiseDaughterCell()
     * is called on the new cell-cycle model.
     */
}

void ProximityLabelRemoverCellCycleModel::UpdateCellCyclePhase()
{
    if ((mCompressedVolumeFraction == DOUBLE_UNSET) || (mEquilibriumVolume == DOUBLE_UNSET))
    {
        EXCEPTION("The member variables mCompressedVolumeFraction and mEquilibriumVolume have not yet been set.");
    }

    // Get cell volume
//    double cell_volume = mpCell->GetCellData()->GetItem("volume");


    // Update G1 duration based on cell volume
//    double compressed_volume = mEquilibriumVolume * mCompressedVolumeFraction;

//    if (cell_volume < compressed_volume)
//    {
//    	/*
//    	 * This method is usually called within a CellBasedSimulation, after the CellPopulation
//    	 * has called CellPropertyRegistry::TakeOwnership(). This means that were we to call
//    	 * CellPropertyRegistry::Instance() here when adding the CellLabel, we would be creating
//    	 * a new CellPropertyRegistry. In this case the CellLabel's cell count would be incorrect.
//    	 * We must therefore access the CellLabel via the cell's CellPropertyCollection.
//    	 */
//    	boost::shared_ptr<AbstractCellProperty> p_label =
//    			mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<CellLabel>();
//    	mpCell->AddCellProperty(p_label);
//    }

    double time_since_birth = GetAge();
    assert(time_since_birth >= 0);

    if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mCurrentCellCyclePhase = G_ZERO_PHASE;
    }
    else if (time_since_birth < GetMDuration())
    {
        mCurrentCellCyclePhase = M_PHASE;
    }
    else if (time_since_birth < GetMDuration() + mG1Duration)
    {
        mCurrentCellCyclePhase = G_ONE_PHASE;
    }
    else if (time_since_birth < GetMDuration() + mG1Duration + GetSDuration())
    {
        mCurrentCellCyclePhase = S_PHASE;
    }
    else if (time_since_birth < GetMDuration() + mG1Duration + GetSDuration() + GetG2Duration())
    {
        mCurrentCellCyclePhase = G_TWO_PHASE;
    }

    // Removes the cell label
    mpCell->RemoveCellProperty<CellLabel>();
}

AbstractCellCycleModel* ProximityLabelRemoverCellCycleModel::CreateCellCycleModel()
{
    return new ProximityLabelRemoverCellCycleModel(*this);
}

void ProximityLabelRemoverCellCycleModel::SetCompressedVolumeFraction(double compressedVolumeFraction)
{
    mCompressedVolumeFraction = compressedVolumeFraction;
}

double ProximityLabelRemoverCellCycleModel::GetCompressedVolumeFraction() const
{
    return mCompressedVolumeFraction;
}

void ProximityLabelRemoverCellCycleModel::SetEquilibriumVolume(double equilibriumVolume)
{
    mEquilibriumVolume = equilibriumVolume;
}

double ProximityLabelRemoverCellCycleModel::GetEquilibriumVolume() const
{
    return mEquilibriumVolume;
}

void ProximityLabelRemoverCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CompressedVolumeFraction>" << mCompressedVolumeFraction << "</CompressedVolumeFraction>\n";
    *rParamsFile << "\t\t\t<EquilibriumVolume>" << mEquilibriumVolume << "</EquilibriumVolume>\n";

    // Call method on direct parent class
    AbstractSimplePhaseBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ProximityLabelRemoverCellCycleModel)
