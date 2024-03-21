// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _COUNTDISTANCEESTIMATIONMAFITERATOR_H_
#define _COUNTDISTANCEESTIMATIONMAFITERATOR_H_

#include "AbstractDistanceEstimationMafIterator.h"

namespace bpp
{
/**
 * @brief Compute A simple distance using observed counts.
 */
class CountDistanceEstimationMafIterator :
  public AbstractDistanceEstimationMafIterator
{
private:
  std::string gapOption_;
  bool unresolvedAsGap_;

public:
  /**
   * @brief Build a new distance estimation maf iterator, based on the SiteContainerTools::computeSimilarityMatrix method.
   *
   * @see SiteContainerTools
   * @param iterator The input iterator.
   * @param gapOption How to deal with gaps. Option forawarded to the computeSimilarityMatrix method.
   * @param unresolvedAsGap Tell if unresolved characters should be considered as gaps. Option forawarded to the computeSimilarityMatrix method.
   * @param addCoordinatesInSequenceNames Should full sequence coordinates be included in the sequence name in the output matrix?
   */
  CountDistanceEstimationMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      const std::string& gapOption,
      bool unresolvedAsGap,
      bool addCoordinatesInSequenceNames = true) :
    AbstractDistanceEstimationMafIterator(iterator, addCoordinatesInSequenceNames),
    gapOption_(gapOption),
    unresolvedAsGap_(unresolvedAsGap)
  {}

public:
  std::string getPropertyName() const { return "CountDistance"; }
  std::unique_ptr<DistanceMatrix> estimateDistanceMatrixForBlock(const MafBlock& block);
};
} // end of namespace bpp.

#endif // _COUNTDISTANCEESTIMATIONMAFITERATOR_H_
