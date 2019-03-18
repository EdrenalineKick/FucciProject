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

#include "FucciProportionalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "Exception.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"

FucciProportionalCellCycleModel::FucciProportionalCellCycleModel()
        : AbstractSimplePhaseBasedCellCycleModel(),
          mG1Rate(1.0 / mTransitCellG1Duration),
          mSRate(1.0 / mSDuration),
          mG2Rate(1.0 / mG2Duration)
{
}

FucciProportionalCellCycleModel::FucciProportionalCellCycleModel(const FucciProportionalCellCycleModel& rModel)
        : AbstractSimplePhaseBasedCellCycleModel(rModel),
          mG1Rate(rModel.mG1Rate),
          mSRate(rModel.mSRate),
          mG2Rate(rModel.mG2Rate)
{
    /*
     * Initialize only the member variable defined in this class.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractSimplePhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mG1Duration is (re)set as soon as InitialiseDaughterCell()
     * is called on the new cell-cycle model.
     */
}

AbstractCellCycleModel* FucciProportionalCellCycleModel::CreateCellCycleModel()
{
    return new FucciProportionalCellCycleModel(*this);
}

void FucciProportionalCellCycleModel::Initialise()
{
    SetG1Duration();
    SetSDuration();
    SetG2Duration();
}

void FucciProportionalCellCycleModel::InitialiseDaughterCell()
{
    AbstractSimplePhaseBasedCellCycleModel::InitialiseDaughterCell();
    SetSDuration();
    SetG2Duration();
}

void FucciProportionalCellCycleModel::ResetForDivision()
{
    AbstractSimplePhaseBasedCellCycleModel::ResetForDivision();
    SetSDuration();
    SetG2Duration();    
}

void FucciProportionalCellCycleModel::SetG1Duration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>()
        || mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
    {
        // Generate an exponential random number with mScale
        mG1Duration = p_gen->ExponentialRandomDeviate(mG1Rate);
        std::cout << "G1 Duration: " << mG1Duration << "\n";
    }
    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }
    else
    {
        NEVER_REACHED;
    }
}

void FucciProportionalCellCycleModel::SetSDuration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    // Generate an exponential random number with mScale
    mSDuration = p_gen->ExponentialRandomDeviate(mSRate);
    std::cout << "S Duration: " << mSDuration << "\n";

}

void FucciProportionalCellCycleModel::SetG2Duration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    // Generate an exponential random number with mScale
    mG2Duration = p_gen->ExponentialRandomDeviate(mG2Rate);
    std::cout << "G2 Duration: " << mG2Duration << "\n";

}

double FucciProportionalCellCycleModel::GetG1Rate()
{
    return mG1Rate;
}

double FucciProportionalCellCycleModel::GetSRate()
{
    return mSRate;
}

double FucciProportionalCellCycleModel::GetG2Rate()
{
    return mG2Rate;
}

void FucciProportionalCellCycleModel::SetG1Rate(double rate)
{
    mG1Rate = rate;
}

void FucciProportionalCellCycleModel::SetSRate(double rate)
{
    mSRate = rate;
}

void FucciProportionalCellCycleModel::SetG2Rate(double rate)
{
    mG2Rate = rate;
}


void FucciProportionalCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<G1Rate>" << mG1Rate << "</G1Rate>\n";
    *rParamsFile << "\t\t\t<SRate>" << mSRate << "</SRate>\n";
    *rParamsFile << "\t\t\t<G2Rate>" << mG2Rate << "</G2Rate>\n";

    // Call method on direct parent class
    AbstractSimplePhaseBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(FucciProportionalCellCycleModel)
