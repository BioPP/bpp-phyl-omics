// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _MAXIMUMLIKELIHOODDISTANCEESTIMATIONMAFITERATOR_H_
#define _MAXIMUMLIKELIHOODDISTANCEESTIMATIONMAFITERATOR_H_

#include "AbstractDistanceEstimationMafIterator.h"

// From bpp-phyl:
#include <Bpp/Phyl/OptimizationTools.h>

namespace bpp
{
/**
 * @brief Compute A simple distance using observed counts.
 */
class MaximumLikelihoodDistanceEstimationMafIterator :
  public AbstractDistanceEstimationMafIterator
{
private:
  std::unique_ptr<DistanceEstimation> distEst_;
  double propGapsToKeep_; // Exclude sites with too many gaps
  bool gapsAsUnresolved_;  // For most models, should be yes as they do not allow for gap characters
  std::string paramOpt_;

public:
  /**
   * @brief Build a new distance estimation maf iterator, based on the DistanceEstimation class.
   *
   * @see DistanceEstimation
   * @param iterator The input iterator.
   * @param distEst A DistanceEstimation object, initialized with an appropriate substitution model.
   * @param propGapsToKeep The maximum gapfrequency in a site to include it in the analysis.
   * @param gapsAsUnresolved Tell if gap characters should be considered as unresolved states. In ost cases it should be set to true, as very few substitution models consider gaps as genuine states.
   * @param paramOpt Tell if substitution model parameters should be optimized in a pairwise manner or not. See OptimizationTools::estimateDistanceMatrix for more details.
   * @param addCoordinatesInSequenceNames Should full sequence coordinates be included in the sequence name in the output matrix?
   * @param verbose Tell if some information should be output in the default message stream.
   */
  MaximumLikelihoodDistanceEstimationMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      std::unique_ptr<DistanceEstimation> distEst,
      double propGapsToKeep = 0,
      bool gapsAsUnresolved = true,
      const string& paramOpt = OptimizationTools::DISTANCEMETHOD_INIT,
      bool addCoordinatesInSequenceNames = true,
      bool verbose = true) :
    AbstractDistanceEstimationMafIterator(iterator, addCoordinatesInSequenceNames),
    distEst_(std::move(distEst)),
    propGapsToKeep_(propGapsToKeep),
    gapsAsUnresolved_(gapsAsUnresolved),
    paramOpt_(paramOpt)
  {
    setVerbose(verbose);
    distEst_->setVerbose(verbose ? 3 : 0);
  }

private:
  MaximumLikelihoodDistanceEstimationMafIterator(const MaximumLikelihoodDistanceEstimationMafIterator& iterator) :
    AbstractDistanceEstimationMafIterator(0, true),
    distEst_(),
    propGapsToKeep_(iterator.propGapsToKeep_),
    gapsAsUnresolved_(iterator.gapsAsUnresolved_),
    paramOpt_(iterator.paramOpt_)
  {}

  MaximumLikelihoodDistanceEstimationMafIterator& operator=(const MaximumLikelihoodDistanceEstimationMafIterator& iterator)
  {
    distEst_.reset();
    propGapsToKeep_ = iterator.propGapsToKeep_;
    gapsAsUnresolved_ = iterator.gapsAsUnresolved_;
    paramOpt_ = iterator.paramOpt_;
    return *this;
  }

public:
  std::string getPropertyName() const { return "MLDistance"; }
  std::unique_ptr<DistanceMatrix> estimateDistanceMatrixForBlock(const MafBlock& block);
};
} // end of namespace bpp.

#endif // _MAXIMUMLIKELIHOODDISTANCEESTIMATIONMAFITERATOR_H_
