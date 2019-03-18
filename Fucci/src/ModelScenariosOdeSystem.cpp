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

#include "ModelScenariosOdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include "MathsCustomFunctions.hpp"

ModelScenariosOdeSystem::ModelScenariosOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(5)
{
    mpSystemInfo = OdeSystemInformation<ModelScenariosOdeSystem>::Instance();

    Init();

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

ModelScenariosOdeSystem::~ModelScenariosOdeSystem()
{
    // Do nothing
}

void ModelScenariosOdeSystem::Init()
{
    // Initialise model parameter values
    u = 0.005;
    k1 = 0.05;
    apcThreshold = 0.2;
    k2a = 0.05;
    k2b = 1;
    k3a = 0.1;
    k3b = 3;
    k4a = 0;
    k4b = 2;
    kas = 0.05;
    kada = 0.005;
    kadb = 1;
    kaa = 1;
    j3 = 0.05;
    j4 = 0.05;
}

void ModelScenariosOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    double x1 = rY[0];
    double x2 = rY[1];
    double x3 = rY[2];
    double x4 = rY[3];
    double x5 = rY[4];

    double dx1 = 0.0;
    double dx2 = 0.0;
    double dx3 = 0.0;
    double dx4 = 0.0;
    double dx5 = 0.0;
 
    double temp1 = 0.0;
    double temp2 = 0.0;
    double temp3 = 0.0;

    /**
     * 1. [CDK]
     * 2. [APC]
     * 3. [ACTt]
     * 4. [ACT]
     * 5. m - mass of the cell
     */

    dx1 = k1*x5 - (k2a*(1-x2)+k2b*x2*x1);

    // The commented line below models the start transition, no cycling, without Cdc20A
//    temp1 = ((mK3d)*(1.0-x2))/(mJ3+1.0-x2);

    temp1 = (((k3a+k3b*x4)*(1-x2))/(j3+1-x2));
    temp2 = ((k4a+k4b*x1)*x2)/(j4+x2);
    dx2 = temp1-temp2;

    temp1 = kas;
    temp2 = (kada*(1-x2)+kadb*x2)*x3;
    dx3 = temp1 - temp2;

    if (rY[1] <= 0.2 && dx2 < 0.0 && rY[0] < 0.8)
    {
        kais = 3.0;
    }
    else
    {
        kais = 0.0;
    }

    if (rY[0] >= 1.0 && dx1 > 0.0)
    {
        kaim = 7.0;
    }
    else
    {
        kaim = 0.0;
    }

    kai = kais + kaim;

    temp1 = kaa*(x3-x4);
    temp2 = kai*x4;
    temp3 = (kada*(1.0-x2)+kadb*x2)*x4;
    dx4 = temp1 - temp2 - temp3;

    dx5 = u*x5;

    // Multiply by 60 beacuase the Tyson and Novak 2001 paper has time in minutes, not hours
    rDY[0] = dx1*60;
    rDY[1] = dx2*60;
    rDY[2] = dx3*60;
    rDY[3] = dx4*60;
    
    
    if ((rY[4] > 1.6 ) && (rY[1] < apcThreshold) && rDY[0] < 0.0 )
    {
        dx5 = -x5/(2);
    }

    rDY[4] = dx5*60;
}


bool ModelScenariosOdeSystem::CalculateStoppingEvent(double time, const std::vector<double>& rY)
{
    std::vector<double> dy(rY.size());
    EvaluateYDerivatives(time, rY, dy);

    // Only call this a stopping condition if the mass of the cell is over 0.6
    // (normally cycles from 0.5-1.0 ish!)
    return ( (rY[4] > 1.6 ) && (rY[1] < apcThreshold) && dy[0] < -0.5 );
}

double ModelScenariosOdeSystem::CalculateRootFunction(double time, const std::vector<double>& rY)
{
    std::vector<double> dy(rY.size());
    EvaluateYDerivatives(time, rY, dy);

    // Only call this a stopping condition if the mass of the cell is over 1.6
    // (normally cycles from 1.0-2.0 ish!)
    if (rY[4]<1.6)
    {
        return 1.0;
    }

    if (dy[0] >= 0.0)
    {
        return 1.0;
    }
    return rY[1]-apcThreshold;
}

template<>
void OdeSystemInformation<ModelScenariosOdeSystem>::Initialise()
{

    this->mVariableNames.push_back("CDK");
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(0.05);

    this->mVariableNames.push_back("APC");
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(0.5);

    this->mVariableNames.push_back("ACTt");
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("ACTa");
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(0.5);

    this->mVariableNames.push_back("mass");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(1.0);

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ModelScenariosOdeSystem)
