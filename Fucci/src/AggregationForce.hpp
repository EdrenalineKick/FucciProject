#ifndef AGGREGATIONFORCE_HPP_
#define AGGREGATIONFORCE_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"

#include "AbstractForce.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"

/**
 * A force class to model random cell movement.
 */
template <unsigned DIM>
class AggregationForce : public AbstractForce<DIM>
{
private :
    /**
     * Archiving.
     */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive& mMovementLikelihood;
        archive& mScalingFactor;
    }

    /**
     * Random Movement Parameter.
     */
    double mMovementLikelihood;
    double mScalingFactor;

public :
    /**
     * Constructor.
     */
    AggregationForce();

    /**
     * Destructor.
     */
    ~AggregationForce();

    /**
     * Set the diffusion constant for the cells.
     *
     * @param MovementLikelihood the movement parameter to use
     */
    void SetMovementLikelihood(double movementLikelihood);

    /**
     * Get the random motion coefficient.
     *
     * @return mMovementLikelihood
     */
    double GetMovementLikelihood();

    void SetScalingFactor(double scalingFactor);

    double GetScalingFactor();

    unsigned GetClosestNeighbour(AbstractCellPopulation<DIM>& rCellPopulation, unsigned index);

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the tissue
     *
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AggregationForce)

#endif /*AGGREGATIONFORCE_HPP_*/