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

#include "G2TysonNovak2001OdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include "MathsCustomFunctions.hpp"

G2TysonNovak2001OdeSystem::G2TysonNovak2001OdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(6)
{
    mpSystemInfo = OdeSystemInformation<G2TysonNovak2001OdeSystem>::Instance();

    Init();

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

G2TysonNovak2001OdeSystem::~G2TysonNovak2001OdeSystem()
{
    // Do nothing
}

void G2TysonNovak2001OdeSystem::Init()
{
    // Initialise model parameter values
    mK1 = 0.04;
    mK2d = 0.04;
    mK2dd = 1.0;
    mK2ddd = 1.0;
    mCycB_threshold = 0.1;
    mK3d = 1.0;
    mK3dd = 10.0;
    mK4d = 2.0;
    mK4 = 35;
    mJ3 = 0.04;
    mJ4 = 0.04;
    mK5d = 0.005;
    mK5dd = 0.2;
    mK6 = 0.1;
    mJ5 = 0.3;
    mN = 4u;
    mK7 = 1.0;
    mK8 = 0.5;
    mJ7 = 1e-3;
    mJ8 = 1e-3;
    mMad = 1.0;
    mK9 = 0.1;
    mK10 = 0.02;
    mK11 = 1.0;
    mK12d = 0.2;
    mK12dd = 50.0;
    mK12ddd = 100.0;
    mKeq = 1e3;
    mK13 = 1.0;
    mK14 = 1.0;
    mK15d = 1.5;
    mK15dd = 0.05;
    mK16d = 1.0;
    mK16dd = 3.0;
    mJ15 = 0.01;
    mJ16 = 0.01;
    mMu = 0.01;
    mMstar = 10.0;
}

void G2TysonNovak2001OdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    double x1 = rY[0];
    double x2 = rY[1];
    double x3 = rY[2];
    double x4 = rY[3];
    double x5 = rY[4];
    double x6 = rY[5];

    double dx1 = 0.0;
    double dx2 = 0.0;
    double dx3 = 0.0;
    double dx4 = 0.0;
    double dx5 = 0.0;
    double dx6 = 0.0;

    double temp1 = 0.0;
    double temp2 = 0.0;
    double temp3 = 0.0;

    /**
     * 1. [CycB]
     * 2. [Cdh1]
     * 3. [Cdc20T]
     * 4. [Cdc20A]
     * 5. [IEP]
     * 6. m - mass of the cell
     */

    dx1 = mK1-(mK2d+mK2dd*x2)*x1;

    // The commented line below models the start transition, no cycling, without Cdc20A
//    temp1 = ((mK3d)*(1.0-x2))/(mJ3+1.0-x2);

    temp1 = ((mK3d+mK3dd*x4)*(1.0-x2))/(mJ3+1.0-x2);
    temp2 = (mK4*x6*x1*x2)/(mJ4+x2);
    dx2 = temp1-temp2;

    temp1 = mK5dd*(SmallPow(x1*x6/mJ5,mN)/(1+SmallPow(x1*x6/mJ5,mN)));
    temp2 = mK6*x3;
    dx3 = mK5d + temp1 - temp2;

    temp1 = (mK7*x5*(x3-x4))/(mJ7+x3-x4);
    temp2 = (mK8*mMad*x4)/(mJ8+x4);
    temp3 = mK6*x4;
    dx4 = temp1 - temp2 - temp3;

    dx5 = mK9*x6*x1*(1.0-x5) - mK10*x5;

    dx6 = mMu*x6*(1.0-x6/mMstar);

    // std::cout << time << " "
    // << rY[0] << " "
    // << rY[1] << " "
    // << rY[2] << " "
    // << rY[3] << " "
    // << rY[4] << " "
    // << rY[5] << "\n";

    // std::cout << rDY[6] << " ";
    // Multiply by 60 beacuase the Tyson and Novak 2001 paper has time in minutes, not hours
    rDY[0] = dx1*10.0;
    rDY[1] = dx2*10.0;
    rDY[2] = dx3*10.0;
    rDY[3] = dx4*10.0;
    rDY[4] = dx5*10.0;
    rDY[5] = dx6*10.0;
}

bool G2TysonNovak2001OdeSystem::CalculateG1Event(double time, const std::vector<double>& rY)
{
    std::vector<double> dy(rY.size());
    EvaluateYDerivatives(time, rY, dy);
    return( (rY[0] > 0.3 ) && (rY[3] < 0.25) && dy[0] > 0.0 );
}

bool G2TysonNovak2001OdeSystem::CalculateStoppingEvent(double time, const std::vector<double>& rY)
{
    std::vector<double> dy(rY.size());
    EvaluateYDerivatives(time, rY, dy);
    
    CalculateG1Event(time, rY);
    // Only call this a stopping condition if the mass of the cell is over 0.6
    // (normally cycles from 0.5-1.0 ish!)

    return ( (rY[5] > 0.6 ) && (rY[0] < mCycB_threshold) && dy[0] < 0.0 );
}

double G2TysonNovak2001OdeSystem::CalculateRootFunction(double time, const std::vector<double>& rY)
{
    std::vector<double> dy(rY.size());
    EvaluateYDerivatives(time, rY, dy);

    // Only call this a stopping condition if the mass of the cell is over 0.6
    // (normally cycles from 0.5-1.0 ish!)
    if (rY[5]<0.6)
    {
        return 1.0;
    }

    if (dy[0] >= 0.0)
    {
        return 1.0;
    }
    return rY[0]-mCycB_threshold;
}

template<>
void OdeSystemInformation<G2TysonNovak2001OdeSystem>::Initialise()
{
    /*
     * Initialise state variables.
     *
     * These initial conditions are the approximate steady state
     * solution values while the commented out conditions are taken
     * from the Tyson and Novak 2001 paper.
     */
    this->mVariableNames.push_back("CycB");
    this->mVariableUnits.push_back("nM");
//    this->mInitialConditions.push_back(0.1);
    this->mInitialConditions.push_back(0.099999999999977);

    this->mVariableNames.push_back("Cdh1");
    this->mVariableUnits.push_back("nM");
//    this->mInitialConditions.push_back(9.8770e-01);
    this->mInitialConditions.push_back(0.989026454281841);

    this->mVariableNames.push_back("Cdc20T");
    this->mVariableUnits.push_back("nM");
//    this->mInitialConditions.push_back(1.5011e+00);
    this->mInitialConditions.push_back(1.547942029285891);

    this->mVariableNames.push_back("Cdc20A");
    this->mVariableUnits.push_back("nM");
//    this->mInitialConditions.push_back(1.2924e+00);
    this->mInitialConditions.push_back(1.421110920155839);

    this->mVariableNames.push_back("IEP");
    this->mVariableUnits.push_back("nM");
//    this->mInitialConditions.push_back(6.5405e-01);
    this->mInitialConditions.push_back(0.672838844290094);

    this->mVariableNames.push_back("mass");
    this->mVariableUnits.push_back("");
//    this->mInitialConditions.push_back(4.7039e-01);
    this->mInitialConditions.push_back(0.970831277863956 / 2);

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(G2TysonNovak2001OdeSystem)
