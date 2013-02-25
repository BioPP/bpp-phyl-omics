//
// File: CountDistanceEstimationMafIterators.h
// Created by: Julien Dutheil
// Created on: Jul 24 2012
//

/*
Copyright or © or Copr. Bio++ Development Team

This software is a computer program whose purpose is to test the
homogeneity of the substitution process of a given alignment.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef _COUNTDISTANCEESTIMATIONMAFITERATOR_H_
#define _COUNTDISTANCEESTIMATIONMAFITERATOR_H_

#include "AbstractDistanceEstimationMafIterator.h"

namespace bpp {

/**
 * @brief Compute A simple distance using observed counts.
 */
class CountDistanceEstimationMafIterator:
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
     */
    CountDistanceEstimationMafIterator(MafIterator* iterator, const std::string& gapOption, bool unresolvedAsGap):
      AbstractDistanceEstimationMafIterator(iterator),
      gapOption_(gapOption), unresolvedAsGap_(unresolvedAsGap)
    {}
    
  public:
    std::string getPropertyName() const { return "CountDistance"; }
    DistanceMatrix* estimateDistanceMatrixForBlock(const MafBlock& block);

};

} //end of namespace bpp.

#endif //_COUNTDISTANCEESTIMATIONMAFITERATOR_H_

