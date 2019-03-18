#include "DirectedMotionForce.hpp"

template <unsigned DIM>
DirectedMotionForce<DIM>::DirectedMotionForce()
        : AbstractForce<DIM>(),
          mMovementParameter(0.01),
          mPreviousDirection(0.0)
{
}

template <unsigned DIM>
DirectedMotionForce<DIM>::~DirectedMotionForce()
{
}

template <unsigned DIM>
void DirectedMotionForce<DIM>::SetMovementParameter(double movementParameter)
{
    assert(movementParameter > 0.0);
    mMovementParameter = movementParameter;
}

template <unsigned DIM>
double DirectedMotionForce<DIM>::GetMovementParameter()
{
    return mMovementParameter;
}

template <unsigned DIM>
void DirectedMotionForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    double dt = SimulationTime::Instance()->GetTimeStep();
    // double PI = 3.14159;

    // Iterate over the nodes
    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        c_vector<double, DIM> force_contribution;
        
        // double force_angle = mPreviousDirection + (RandomNumberGenerator::Instance()->StandardNormalRandomDeviate())* * PI;
        double force_angle = mPreviousDirection + (RandomNumberGenerator::Instance()->StandardNormalRandomDeviate())* 0.001;

        double xi = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
        double strengthX = (sqrt(2.0 * mMovementParameter * dt) / dt) * cos(force_angle) * xi;

        xi = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
        double strengthY = (sqrt(2.0 * mMovementParameter * dt) / dt) * sin(force_angle) * xi;
        force_contribution[0] = strengthX;
        force_contribution[1] = strengthY;

        node_iter->AddAppliedForceContribution(force_contribution);
        mPreviousDirection = force_angle;
    }
}

template <unsigned DIM>
void DirectedMotionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MovementParameter>" << mMovementParameter << "</MovementParameter> \n";

    // Call direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class DirectedMotionForce<1>;
template class DirectedMotionForce<2>;
template class DirectedMotionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DirectedMotionForce)