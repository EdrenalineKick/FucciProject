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

#ifndef ATTILAHYBRID2006ODESYSTEM_HPP_
#define ATTILAHYBRID2006ODESYSTEM_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"

#include <cmath>
#include <iostream>
#include "RandomNumberGenerator.hpp"

#include "AbstractOdeSystem.hpp"

/**
 * Represents the Tyson & Novak (2001) system of ODEs.
 * [doi:10.1006/jtbi.2001.2293]
 */
class AttilaHybrid2006OdeSystem : public AbstractOdeSystem
{
private:
    double mV;

    double mKwr;

    double mJwr;

    double mKw;

    double mJw;

    double mK3a;

    double mK3b;

    double mJ3;

    double mK4a;

    double mK4b;

    double mJ4;

    double mK25;

    double mJ25;

    double mK25r;

    double mJ25r;

    double mVawee;

    double mVbwee;

    double mK2a;

    double mK2b;

    double mVac25;

    double mVbc25;

    double mKs20p;

    double mJ20;

    double mKs20pp;

    unsigned mN;

    double mKa20;

    double mJa20;

    double mK1;

    double mKi20;

    double mJi20;

    double mKCdc20i;

    bool mG1;

    double mG1EventTime;

    double mOdeStopTime;

    double mVRand;

    double mK4bRand;

    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

public:
    /**
     * Constructor.
     *
     * @param stateVariables optional initial conditions for state variables (only used in archiving)
     */
    AttilaHybrid2006OdeSystem(std::vector<double> stateVariables = std::vector<double>());

    /**
     * Destructor.
     */
    ~AttilaHybrid2006OdeSystem();

    /**
     * Initialise parameter values.
     */
    void Init();

    /**
     * Compute the RHS of the Alarcon et al. (2004) system of ODEs.
     *
     * Returns a vector representing the RHS of the ODEs at each time step, y' = [y1' ... yn'].
     * An ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *
     * @param time used to evaluate the RHS.
     * @param rY value of the solution vector used to evaluate the RHS.
     * @param rDY filled in with the resulting derivatives (using Alarcons et al. (2004) system of equations).
     */
    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);

    bool CalculateG1Event(double time, const std::vector<double>& rY);

    /**
     * Calculate whether the conditions for the cell cycle to finish have been met.
     * (Used by Chaste solvers to find whether or not to stop solving)
     *
     * @param time at which to calculate whether the stopping event has occurred
     * @param rY value of the solution vector used to evaluate the RHS.
     *
     * @return whether or not stopping conditions have been met
     */
    bool CalculateStoppingEvent(double time, const std::vector<double>& rY);

    /**
     * Calculate whether the conditions for the cell cycle to finish have been met.
     * (Used by CVODE solver to find exact stopping position)
     *
     * @param time at which to calculate whether the stopping event has occurred
     * @param rY value of the solution vector used to evaluate the RHS.
     *
     * @return How close we are to the root of the stopping condition
     */
    double CalculateRootFunction(double time, const std::vector<double>& rY);

    double GetG1EventTime();

    void SetG1EventTime(double g1eventtime);

    double GetStopTime();

    void SetmV(double mv);
    void SetmVRand(double mvrand);
    void SetmK4b(double mk4b);
    void SetmK4bRand(double mk4brand);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(AttilaHybrid2006OdeSystem)

namespace boost
{
namespace serialization
{
    /**
 * Serialize information required to construct a AttilaHybrid2006OdeSystem.
 */
    template <class Archive>
    inline void save_construct_data(
        Archive& ar, const AttilaHybrid2006OdeSystem* t, const unsigned int file_version)
    {
        // Save data required to construct instance
        const std::vector<double>& state_variables = t->rGetConstStateVariables();
        ar& state_variables;
    }

    /**
 * De-serialize constructor parameters and initialise a AttilaHybrid2006OdeSystem.
 */
    template <class Archive>
    inline void load_construct_data(
        Archive& ar, AttilaHybrid2006OdeSystem* t, const unsigned int file_version)
    {
        // Retrieve data from archive required to construct new instance
        std::vector<double> state_variables;
        ar& state_variables;

        // Invoke inplace constructor to initialise instance
        ::new (t) AttilaHybrid2006OdeSystem(state_variables);
    }
} // namespace serialization
} // namespace boost

#endif /*ATTILAHYBRID2006ODESYSTEM_HPP_*/
