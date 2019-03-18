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

#include "NormalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

NormalCellCycleModel::NormalCellCycleModel()
    : AbstractSimplePhaseBasedCellCycleModelModified(),
      mMinCellCycleDuration(13.0), // Hours
      mMaxCellCycleDuration(40.0)  // Hours
{
}

NormalCellCycleModel::NormalCellCycleModel(const NormalCellCycleModel& rModel)
   : AbstractSimplePhaseBasedCellCycleModelModified(rModel),
     mMinCellCycleDuration(rModel.mMinCellCycleDuration),
     mMaxCellCycleDuration(rModel.mMaxCellCycleDuration)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variable mCellCycleDuration is initialized in the
     * AbstractSimpleCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mCellCycleDuration is (re)set as soon as
     * InitialiseDaughterCell() is called on the new cell-cycle model.
     */
}

AbstractCellCycleModel* NormalCellCycleModel::CreateCellCycleModel()
{
    SetCellCycleDuration();
    SetG1Duration();
    SetG2Duration();
    return new NormalCellCycleModel(*this);
}

void NormalCellCycleModel::Initialise()
{
    SetCellCycleDuration();
}

void NormalCellCycleModel::ResetForDivision()
{
    AbstractPhaseBasedCellCycleModel::ResetForDivision();
    SetCellCycleDuration();
}


void NormalCellCycleModel::InitialiseDaughterCell()
{
    AbstractCellCycleModel::InitialiseDaughterCell();
    SetCellCycleDuration();
}

void NormalCellCycleModel::SetCellCycleDuration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
   
    while(mCellCycleDuration < mMinCellCycleDuration || mCellCycleDuration > mMaxCellCycleDuration)
    {
        mSDuration = 0.1;
        //Sampled
        mCellCycleDuration = p_gen->NormalRandomDeviate(23.6618,6.31478);
        //Simulated
        // mCellCycleDuration = p_gen->NormalRandomDeviate(20.85,2.48);
        // std::cout << "CC Length: " << mCellCycleDuration << "\n";
    }
    SetG1Duration();
    SetG2Duration();
}

double NormalCellCycleModel::GetMinCellCycleDuration()
{
    return mMinCellCycleDuration;
}

void NormalCellCycleModel::SetMinCellCycleDuration(double minCellCycleDuration)
{
    mMinCellCycleDuration = minCellCycleDuration;
}

double NormalCellCycleModel::GetMaxCellCycleDuration()
{
    return mMaxCellCycleDuration;
}

void NormalCellCycleModel::SetG1Duration()
{
    mG1Duration = 0.3* mCellCycleDuration;
    // std::cout << "G1: " << mG1Duration << "\n";
}

void NormalCellCycleModel::SetG2Duration()
{
    mG2Duration = 0.7 * mCellCycleDuration;
    // std::cout << "G2: " << mG2Duration << "\n";
}

double NormalCellCycleModel::GetAverageTransitCellCycleTime()
{
    return(23.66);
}

double NormalCellCycleModel::GetAverageStemCellCycleTime()
{
    return(23.66);
}

void NormalCellCycleModel::SetMaxCellCycleDuration(double maxCellCycleDuration)
{
    mMaxCellCycleDuration = maxCellCycleDuration;
}

void NormalCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MinCellCycleDuration>" << mMinCellCycleDuration << "</MinCellCycleDuration>\n";
    *rParamsFile << "\t\t\t<MaxCellCycleDuration>" << mMaxCellCycleDuration << "</MaxCellCycleDuration>\n";

    // Call method on direct parent class
    AbstractSimplePhaseBasedCellCycleModelModified::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(NormalCellCycleModel)
