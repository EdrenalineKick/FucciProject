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

#include "AttilaHybrid2006OdeSystem.hpp"
#include "MathsCustomFunctions.hpp"
#include "OdeSystemInformation.hpp"

AttilaHybrid2006OdeSystem::AttilaHybrid2006OdeSystem(std::vector<double> stateVariables)
        : AbstractOdeSystem(7)
{
    mpSystemInfo = OdeSystemInformation<AttilaHybrid2006OdeSystem>::Instance();

    Init();

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

AttilaHybrid2006OdeSystem::~AttilaHybrid2006OdeSystem()
{
    // Do nothing
}

void AttilaHybrid2006OdeSystem::Init()
{
    //     mV = 0.14 + RandomNumberGenerator::Instance()->ranf()*0.05;

    //     mK4b = 1.5 + RandomNumberGenerator::Instance()->ranf() * 1.5;

    //Currently used parameters
    // mV = 0.16 + RandomNumberGenerator::Instance()->ranf() * 0.05;

    // mK4b = 2.0 + RandomNumberGenerator::Instance()->ranf() * 1.0;


    //New Fitted parameters
    // mV = 0.0926 + RandomNumberGenerator::Instance()->ranf() * 0.194;
    // mK4b = 2.433 + RandomNumberGenerator::Instance()->ranf() * 0.974;
    // mK4b = 2.433 + RandomNumberGenerator::Instance()->ranf() * 1.2;

    // Newer with weight towards stopTime
    mV = 0.092595 + 0.186332 * RandomNumberGenerator::Instance()->ranf();
    mK4b = 2.3227 + 1.20598 *RandomNumberGenerator::Instance()->ranf();

    // std::cout << "mVRand: " << mVRand << "\n" << "mK4bRand" << mK4bRand << "\n";
    // mV = 0.437587 * 0.35;
    // mK4b = 0.8548 *3;
    std::cout << mV << " " << mK4b << "\n";
    mKwr = 0.2;
    mJwr = 0.1;
    mKw = 0.5;
    mJw = 0.1;
    mK3a = 0.1;
    mK3b = 3;
    mJ3 = 0.01;
    // mK4a = 0;
    mK4a = 0.1;
    // mK4b = 2;
    mJ4 = 0.05;
    mK25 = 0.5;
    mJ25 = 0.1;
    mK25r = 0.2;
    mJ25r = 0.1;
    mVawee = 0.01;
    mVbwee = 0.8;
    mK2a = 0.05;
    // mK2b = 1;
    mK2b = 0.8;
    mVac25 = 0.02;
    mVbc25 = 0.5;
    mKs20p = 0.001;
    mJ20 = 1;
    mKs20pp = 0.5;
    mN = 4;
    mKa20 = 0.5;
    mJa20 = 0.1;
    mK1 = 0.15;
    mKi20 = 0.25;
    mJi20 = 0.05;
    mKCdc20i = 0.5;
    mG1 = false;
    mG1EventTime = DBL_MAX;
    mOdeStopTime = DBL_MAX;
}

void AttilaHybrid2006OdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    double x1 = rY[0];
    double x2 = rY[1];
    double x3 = rY[2];
    double x4 = rY[3];
    double x5 = rY[4];
    double x6 = rY[5];
    double x7 = rY[6];

    double dx1 = 0.0;
    double dx2 = 0.0;
    double dx3 = 0.0;
    double dx4 = 0.0;
    double dx5 = 0.0;
    double dx6 = 0.0;
    double dx7 = 0.0;

    double temp1 = 0.0;
    double temp2 = 0.0;
    double temp3 = 0.0;
    double temp4 = 0.0;
    /**
     * 1. [CDKa]
     * 2. [CDKp]
     * 3. [Wee1]
     * 4. [Cdc25]
     * 5. [Cdc20i]
     * 6. [Cdc20a]
     * 7. [Fzr1]
     */

    // std::cout << mV << " " << mK4b << "\n";

    temp1 = (mVac25 * (1 - x4) + mVbc25 * x4) * x2;
    temp2 = (mVawee * (1 - x3) + mVbwee * x3) * x1;
    temp3 = (mK2a * (1 - x7) + mK2b * x7) * x1;

    dx1 = mV + temp1 - temp2 - temp3;

    temp3 = (mK2a * (1 - x7) + mK2b * x7) * x2;

    dx2 = temp2 - temp1 - temp3;

    temp1 = (mKwr * (1 - x3)) / (mJwr + 1 - x3);
    temp2 = (mKw * x1 * x3) / (mJw + x3);
    dx3 = temp1 - temp2;

    temp1 = (mK25 * x1 * (1 - x4)) / (mJ25 + 1 - x4);
    temp2 = (mK25r * x4) / (mJ25r + x4);
    dx4 = temp1 - temp2;

    temp1 = (mKs20p + mKs20pp * SmallPow(x1, mN)) / (SmallPow(mJ20, mN) + SmallPow(x1, mN));
    temp2 = (mKi20 * x6) / (mJi20 + x6);
    temp3 = (mKa20 * x1 * x5) / (mJa20 + x5);
    temp4 = mKCdc20i * x5;
    dx5 = temp1 + temp2 - temp3 - temp4;

    temp1 = mK1 * x6;
    dx6 = temp3 - temp1 - temp2;

    temp1 = ((mK3a + mK3b * x6) * (1 - x7)) / (mJ3 + 1 - x7);
    temp2 = ((mK4a + mK4b * x1) * x7) / (mJ4 + x7);
    dx7 = temp1 - temp2;

    // std::cout << time << " "
    // << rY[0] << " "
    // << rY[1] << " "
    // << rY[2] << " "
    // << rY[3] << " "
    // << rY[4] << " "
    // << rY[5] << " "
    // << rY[6] << "\n";

    // std::cout << rDY[6] << " ";
    // Multiply by 60 beacuase the Tyson and Novak 2001 paper has time in minutes, not hours
    // double scalingFactor = 1.5;
    double scalingFactor = 1.2;
    rDY[0] = dx1 * scalingFactor;
    rDY[1] = dx2 * scalingFactor;
    rDY[2] = dx3 * scalingFactor;
    rDY[3] = dx4 * scalingFactor;
    rDY[4] = dx5 * scalingFactor;
    rDY[5] = dx6 * scalingFactor;
    rDY[6] = dx7 * scalingFactor;
}

bool AttilaHybrid2006OdeSystem::CalculateG1Event(double time, const std::vector<double>& rY)
{
    std::vector<double> dy(rY.size());
    EvaluateYDerivatives(time, rY, dy);

    // if ((rY[3] > 0.2) && (dy[3] > 0.0) && (rY[1] > 0.3) && (rY[2] > 0.3))
    if ((dy[1] > 0) && (rY[6] < 0.3) && (dy[2] < 0.0))
    {
        if (mG1 == false)
        {
            mG1EventTime = time;
        }
        mG1 = true;
    }

    // return ((rY[3] > 0.2) && (dy[3] > 0.0) && (rY[1] > 0.3) && (rY[2] > 0.3));
    // return ((dy[1] > 0) && (rY[6] < 0.3) && (rY[2] > 0.3));
    return ((dy[1] > 0) && (rY[6] < 0.3) && (dy[2] < 0.0));
}

bool AttilaHybrid2006OdeSystem::CalculateStoppingEvent(double time, const std::vector<double>& rY)
{
    std::vector<double> dy(rY.size());
    EvaluateYDerivatives(time, rY, dy);

    CalculateG1Event(time, rY);
    // Only call this a stopping condition if the mass of the cell is over 0.6
    // (normally cycles from 0.5-1.0 ish!)

    if (mG1 == true)
    {
        if (((rY[0] < 0.26) && (dy[0] < 0.0)) && (rY[6] > 0.8))
        {
            mG1 = false;
            mOdeStopTime = time;
            std::cout << "mV: " << mV << "\n";
            std::cout << "mK4b: " << mK4b << "\n";
        }
        return (((rY[0] < 0.26) && (dy[0] < 0.0)) && (rY[6] > 0.8));
    }
    else
    {
        return 0;
    }
}

double AttilaHybrid2006OdeSystem::CalculateRootFunction(double time, const std::vector<double>& rY)
{
    std::vector<double> dy(rY.size());
    EvaluateYDerivatives(time, rY, dy);

    CalculateG1Event(time, rY);

    if (mG1 == false)
    {
        return 1.0;
    }

    if (dy[0] >= 0.0)
    {
        return 1.0;
    }

    if (rY[6] <= 0.8)
    {
        return 1.0;
    }

    if (rY[0] < 0.26)
    {
        mOdeStopTime = time;
        mG1 = false;
        std::cout << "Root Calculated" << "\n";
    }

    return rY[0] - 0.26;
}

double AttilaHybrid2006OdeSystem::GetG1EventTime()
{
    return mG1EventTime;
}

double AttilaHybrid2006OdeSystem::GetStopTime()
{
    return mOdeStopTime;
}

void AttilaHybrid2006OdeSystem::SetmV(double mv)
{
    mV = mv;
}

void AttilaHybrid2006OdeSystem::SetmVRand(double mvrand)
{
    mVRand = mvrand;
}

void AttilaHybrid2006OdeSystem::SetG1EventTime(double g1eventtime)
{
    mG1EventTime = g1eventtime;
}

void AttilaHybrid2006OdeSystem::SetmK4b(double mk4b)
{
    mK4b = mk4b;
}

void AttilaHybrid2006OdeSystem::SetmK4bRand(double mk4brand)
{
    mK4bRand = mk4brand;
}

template <>
void OdeSystemInformation<AttilaHybrid2006OdeSystem>::Initialise()
{
    this->mVariableNames.push_back("CDKa");
    this->mVariableUnits.push_back("nM");
    //    this->mInitialConditions.push_back(0.1);
    this->mInitialConditions.push_back(0.26);

    this->mVariableNames.push_back("CDKp");
    this->mVariableUnits.push_back("nM");
    //    this->mInitialConditions.push_back(9.8770e-01);
    this->mInitialConditions.push_back(0.0303);

    this->mVariableNames.push_back("Wee1");
    this->mVariableUnits.push_back("nM");
    //    this->mInitialConditions.push_back(1.5011e+00);
    this->mInitialConditions.push_back(0.1825);

    this->mVariableNames.push_back("Cdc25");
    this->mVariableUnits.push_back("nM");
    //    this->mInitialConditions.push_back(1.2924e+00);
    this->mInitialConditions.push_back(0.8175);

    this->mVariableNames.push_back("Cdc20i");
    this->mVariableUnits.push_back("nM");
    //    this->mInitialConditions.push_back(6.5405e-01);
    this->mInitialConditions.push_back(0.2379);

    this->mVariableNames.push_back("Cdc20a");
    this->mVariableUnits.push_back("nM");
    //    this->mInitialConditions.push_back(4.7039e-01);
    this->mInitialConditions.push_back(0.6136);

    this->mVariableNames.push_back("Fzr1");
    this->mVariableUnits.push_back("nM");
    //    this->mInitialConditions.push_back(4.7039e-01);
    this->mInitialConditions.push_back(1.0);

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(AttilaHybrid2006OdeSystem)
