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
#include "Attila2006CellCycleModel.hpp"

Attila2006CellCycleModel::Attila2006CellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
        : AbstractImprovedOdeBasedPhaseBasedCellCycleModel(SimulationTime::Instance()->GetTime(), pOdeSolver)
{
    if (!mpOdeSolver)
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<Attila2006CellCycleModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        // Chaste solvers always check for stopping events, CVODE needs to be instructed to do so
        mpOdeSolver->CheckForStoppingEvents();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<Attila2006CellCycleModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
#endif //CHASTE_CVODE
    }
}

Attila2006CellCycleModel::~Attila2006CellCycleModel()
{
}

Attila2006CellCycleModel::Attila2006CellCycleModel(const Attila2006CellCycleModel& rModel)
        : AbstractImprovedOdeBasedPhaseBasedCellCycleModel(rModel)
{
    /*
     * Initialize only those member variables defined in this class.
     * Create the new cell-cycle model's ODE system and use the current
     * values of the state variables in mpOdeSystem as an initial condition.
     *
     * The member variables mDivideTime and mG2PhaseStartTime are
     * initialized in the AbstractOdeBasedPhaseBasedCellCycleModel
     * constructor.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractPhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     */
    assert(rModel.GetOdeSystem());
    SetOdeSystem(new AttilaHybrid2006OdeSystem);
    SetStateVariables(rModel.GetOdeSystem()->rGetStateVariables());
    mG1Time = DOUBLE_UNSET;
}

AbstractCellCycleModel* Attila2006CellCycleModel::CreateCellCycleModel()
{
    // Create a new cell-cycle model
    return new Attila2006CellCycleModel(*this);
}

void Attila2006CellCycleModel::Initialise()
{
    assert(mpOdeSystem == nullptr);
    assert(mpCell != nullptr);

    mpOdeSystem = new AttilaHybrid2006OdeSystem;
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());

    AbstractImprovedOdeBasedPhaseBasedCellCycleModel::Initialise();
    SetSDuration(0.1);
}

void Attila2006CellCycleModel::ResetForDivision()
{
    AbstractImprovedOdeBasedPhaseBasedCellCycleModel::ResetForDivision();

    assert(mpOdeSystem != nullptr);

    /**
     * This model needs the protein concentrations and phase resetting to G0/G1.
     *
     * In theory, the solution to the Tyson-Novak equations should exhibit stable
     * oscillations, and we only need to halve the mass of the cell each period.
     *
     * However, the backward Euler solver used to solve the equations
     * currently returns a solution that diverges after long times, so
     * we must reset the initial conditions each period.
     *
     * When running with CVODE however we can use the halving the mass of the cell method.
     */
#ifdef CHASTE_CVODE
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
#else
    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
#endif //CHASTE_CVODE
    mG1Time = DOUBLE_UNSET;
}

void Attila2006CellCycleModel::InitialiseDaughterCell()
{
    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        /*
         * This method is usually called within a CellBasedSimulation, after the CellPopulation
         * has called CellPropertyRegistry::TakeOwnership(). This means that were we to call
         * CellPropertyRegistry::Instance() here when setting the CellProliferativeType, we
         * would be creating a new CellPropertyRegistry. In this case the cell proliferative
         * type counts, as returned by AbstractCellPopulation::GetCellProliferativeTypeCount(),
         * would be incorrect. We must therefore access the CellProliferativeType via the cell's
         * CellPropertyCollection.
         */
        boost::shared_ptr<AbstractCellProperty> p_stem_type = mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<StemCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_stem_type);
    }
    SetSDuration(0.1);
}

void Attila2006CellCycleModel::UpdateCellCyclePhase()
{
    assert(mpOdeSystem != nullptr);

    double current_time = SimulationTime::Instance()->GetTime();
    bool mOdeFail = false;

    // Update the phase from M to G1 when necessary
    if (mCurrentCellCyclePhase == M_PHASE)
    {
        double m_duration = GetMDuration();
        if (GetAge() >= m_duration)
        {
            static_cast<AttilaHybrid2006OdeSystem*>(mpOdeSystem)->SetG1EventTime(DBL_MAX);
            mCurrentCellCyclePhase = G_ONE_PHASE;
            mLastTime = m_duration + mBirthTime;
        }
        else
        {
            // Still dividing; don't run ODEs
            return;
        }
    }

    if (current_time > mLastTime)
    {
        if (!this->mFinishedRunningOdes)
        {
            // Update whether a stopping event has occurred
            this->mFinishedRunningOdes = SolveOdeToTime(current_time);
            UpdateG1Time(*this);

            // Check no concentrations have gone negative
            for (unsigned i = 0; i < mpOdeSystem->GetNumberOfStateVariables(); i++)
            {
                if (mpOdeSystem->rGetStateVariables()[i] < -DBL_EPSILON)
                {
                    // LCOV_EXCL_START
                    EXCEPTION("A protein concentration " << i << " has gone negative (" << mpOdeSystem->rGetStateVariables()[i] << ")\n"
                                                         << "Chaste predicts that the CellCycleModel numerical method is probably unstable.");
                    // LCOV_EXCL_STOP
                }
            }

            if (!this->mFinishedRunningOdes)
            {
                //Constantly update G1 ending time during G1
                if (mCurrentCellCyclePhase == G_ONE_PHASE)
                {
                    mG1Time = static_cast<AttilaHybrid2006OdeSystem*>(mpOdeSystem)->GetG1EventTime();
                }
                // std::cout << mG1Duration << "\n";
                // std::cout << "Age: " << GetAge() << "\n";
                // std::cout << "G1 Time: " << mG1Time << "\n";

                //If current time is greater than the G1 ending time, update the cell cycle phase to the new phase.
                if (current_time >= mG1Time && mCurrentCellCyclePhase == G_ONE_PHASE)
                {
                    mCurrentCellCyclePhase = S_PHASE;
                }

                double s_duration = GetSDuration();
                if (current_time >= (mG1Time + s_duration) && mCurrentCellCyclePhase == S_PHASE)
                {
                    mCurrentCellCyclePhase = G_TWO_PHASE;
                }
            }

            if (this->mFinishedRunningOdes && mOdeFail == false)
            {
                // Update durations of each phase
                mG1Duration = static_cast<AttilaHybrid2006OdeSystem*>(mpOdeSystem)->GetG1EventTime() - mBirthTime - GetMDuration();
                mG2PhaseStartTime = static_cast<AttilaHybrid2006OdeSystem*>(mpOdeSystem)->GetG1EventTime() + GetSDuration();
                mDivideTime = GetOdeStopTime();
                mG2Duration = mDivideTime - mG2PhaseStartTime;
                double ccLength = GetOdeStopTime() - mBirthTime;
                std::cout << "Birth Time: " << mBirthTime << "\n"
                          << "G1 Time: " << static_cast<AttilaHybrid2006OdeSystem*>(mpOdeSystem)->GetG1EventTime() << "\n"
                          << "G1: " << mG1Duration << "\n"
                          << "G2 Start: " << mG2PhaseStartTime << "\n"
                          << "Divide Time" << mDivideTime << "\n"
                          << "G2 Length: " << mG2Duration << "\n"
                          << "CC Length: " << ccLength << "\n";

                std::ofstream myfile;

                myfile.open ("cell_cycle_times.csv", std::ios::app);
                myfile << mG1Duration << ",," << mG2Duration << ",\n";

            }
        }
        else
        {
            // ODE model finished, just increasing time until division...
            if (current_time >= mG2PhaseStartTime)
            {
                // std::cout << "ODE finished \n";
            }
        }
    }
}

void Attila2006CellCycleModel::UpdateG1Time(Attila2006CellCycleModel& rModel)
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    //Consult the ODE system to determine whether the event ending G1 has occured.
    //G1 Trigger records that the event has occured so it only updates a G1 time once.

    bool g1_true = static_cast<AttilaHybrid2006OdeSystem*>(mpOdeSystem)->CalculateG1Event(SimulationTime::Instance()->GetTime(), mpOdeSystem->rGetStateVariables());

    if (g1_true == true)
    {
        rModel.G1Trigger();
        rModel.GetG1Time();
    }
}

void Attila2006CellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractImprovedOdeBasedPhaseBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Attila2006CellCycleModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Attila2006CellCycleModel)
