#include "AggregationForce.hpp"

template <unsigned DIM>
AggregationForce<DIM>::AggregationForce()
        : AbstractForce<DIM>(),
          mMovementLikelihood(0.4),
          mScalingFactor(0.6)
{
}

template <unsigned DIM>
AggregationForce<DIM>::~AggregationForce()
{
}

template <unsigned DIM>
void AggregationForce<DIM>::SetMovementLikelihood(double movementLikelihood)
{
    assert(movementLikelihood > 0.0);
    mMovementLikelihood = movementLikelihood;
}

template <unsigned DIM>
double AggregationForce<DIM>::GetMovementLikelihood()
{
    return mMovementLikelihood;
}

template <unsigned DIM>
double AggregationForce<DIM>::GetScalingFactor()
{
    return mScalingFactor;
}

template <unsigned DIM>
void AggregationForce<DIM>::SetScalingFactor(double scalingFactor)
{
    assert(scalingFactor > 0.0);
    mScalingFactor = scalingFactor;
}

template <unsigned DIM>
unsigned AggregationForce<DIM>::GetClosestNeighbour(AbstractCellPopulation<DIM>& rCellPopulation, unsigned index)
{
    double closestDistance = DBL_MAX;
    unsigned closestNeighbour;
    unsigned index_a = index;

    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned index_b = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

        Node<DIM>* p_node_a = rCellPopulation.GetNode(index_a);
        Node<DIM>* p_node_b = rCellPopulation.GetNode(index_b);

        CellPtr p_cell_a = rCellPopulation.GetCellUsingLocationIndex(index_a);
        CellPtr p_cell_b = rCellPopulation.GetCellUsingLocationIndex(index_b);
        unsigned cell_colony = p_cell_a->GetCellData()->GetItem("Colony ID");
        unsigned neighbour_colony = p_cell_b->GetCellData()->GetItem("Colony ID");

        if (p_node_a != p_node_b && cell_colony != neighbour_colony)
        {
            const c_vector<double, DIM>& r_node_a_location = p_node_a->rGetLocation();
            const c_vector<double, DIM>& r_node_b_location = p_node_b->rGetLocation();

            c_vector<double, DIM> difference;
            difference = rCellPopulation.rGetMesh().GetVectorFromAtoB(r_node_a_location, r_node_b_location);

            double distance_between_nodes = norm_2(difference);
            if (distance_between_nodes > 0)
            {
                assert(!std::isnan(distance_between_nodes));

                if (distance_between_nodes < closestDistance)
                {
                    closestDistance = distance_between_nodes;
                    closestNeighbour = index_b;
                }
            }
        }
    }

    if (closestDistance > 3.0)
    {
        return index;
    }
    else
    {
        return closestNeighbour;
    }
}

template <unsigned DIM>
void AggregationForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        double randomChance = RandomNumberGenerator::Instance()->ranf();
        if (randomChance < mMovementLikelihood)
        {
            c_vector<double, DIM> force_contribution;
            CellPtr p_cell_a = *cell_iter;

            // if (p_cell_a->template HasCellProperty<CellLabel>())
            // {
            unsigned cell_index = rCellPopulation.GetLocationIndexUsingCell(p_cell_a);
            unsigned neighbour_index = GetClosestNeighbour(rCellPopulation, cell_index);

            if (cell_index != neighbour_index)
            {
                Node<DIM>* p_node_a = rCellPopulation.GetNode(cell_index);
                Node<DIM>* p_node_b = rCellPopulation.GetNode(neighbour_index);

                const c_vector<double, DIM>& r_node_a_location = p_node_a->rGetLocation();
                const c_vector<double, DIM>& r_node_b_location = p_node_b->rGetLocation();

                c_vector<double, DIM> unit_difference;
                unit_difference = rCellPopulation.rGetMesh().GetVectorFromAtoB(r_node_a_location, r_node_b_location);

                double distance_between_nodes = norm_2(unit_difference);
                assert(distance_between_nodes > 0);
                assert(!std::isnan(distance_between_nodes));

                unit_difference /= distance_between_nodes;
                force_contribution = unit_difference * mScalingFactor;
                rCellPopulation.GetNode(cell_index)->AddAppliedForceContribution(force_contribution);
            }
            // }
        }
    }
}

template <unsigned DIM>
void AggregationForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MovementLikelihood>" << mMovementLikelihood << "</MovementLikelihood> \n";
    *rParamsFile << "\t\t\t<ScalingFactor>" << mScalingFactor << "</ScalingFactor> \n";

    // Call direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class AggregationForce<1>;
template class AggregationForce<2>;
template class AggregationForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AggregationForce)