//
// File: MaximumLikelihoodDistanceEstimationMafIterators.h
// Created by: Julien Dutheil
// Created on: Nov 13 2012
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

#ifndef _MAXIMUMLIKELIHOODDISTANCEESTIMATIONMAFITERATOR_H_
#define _MAXIMUMLIKELIHOODDISTANCEESTIMATIONMAFITERATOR_H_

#include "AbstractDistanceEstimationMafIterator.h"

//From bpp-phyl:
#include <Bpp/Phyl/OptimizationTools.h>

namespace bpp {

/**
 * @brief Compute A simple distance using observed counts.
 */
class MaximumLikelihoodDistanceEstimationMafIterator:
  public AbstractDistanceEstimationMafIterator
{
  private:
    std::unique_ptr<DistanceEstimation> distEst_;
    double propGapsToKeep_; //Exclude sites with too many gaps
    bool gapsAsUnresolved_;  //For most models, should be yes as they do not allow for gap characters
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
    MaximumLikelihoodDistanceEstimationMafIterator(MafIterator* iterator,
        DistanceEstimation* distEst,
        double propGapsToKeep = 0,
        bool gapsAsUnresolved = true,
        const string& paramOpt = OptimizationTools::DISTANCEMETHOD_INIT,
        bool addCoordinatesInSequenceNames = true,
        bool verbose = true):
      AbstractDistanceEstimationMafIterator(iterator, addCoordinatesInSequenceNames),
      distEst_(distEst), propGapsToKeep_(propGapsToKeep), gapsAsUnresolved_(gapsAsUnresolved), paramOpt_(paramOpt)
    {
      setVerbose(verbose);
      distEst_->setVerbose(verbose ? 3 : 0);
    }

  private:
    MaximumLikelihoodDistanceEstimationMafIterator(const MaximumLikelihoodDistanceEstimationMafIterator& iterator):
      AbstractDistanceEstimationMafIterator(0, true),
      distEst_(), propGapsToKeep_(iterator.propGapsToKeep_),
      gapsAsUnresolved_(iterator.gapsAsUnresolved_), paramOpt_(iterator.paramOpt_)
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
    DistanceMatrix* estimateDistanceMatrixForBlock(const MafBlock& block);

};

} //end of namespace bpp.

#endif //_MAXIMUMLIKELIHOODDISTANCEESTIMATIONMAFITERATOR_H_

