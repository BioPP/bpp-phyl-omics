//
// File: MaximumLikelihoodModelFitMafStatistics.cpp
// Created by: Julien Dutheil
// Created on: Mar 25 2014
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

#include "MaximumLikelihoodModelFitMafStatistics.h"

//From bpp-seq:
#include <Bpp/Seq/Container/SiteContainerTools.h>

//From bpp-phyl:
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/OptimizationTools.h>

using namespace bpp;
using namespace std;

const string MaximumLikelihoodModelFitMafStatistics::NO_PROPERTY = "RESERVED_NOPROPERTY";

void MaximumLikelihoodModelFitMafStatistics::compute(const MafBlock& block)
{
  //First we get the alignment:
  unique_ptr<SiteContainer> sites(SiteContainerTools::removeGapSites(block.getAlignment(), propGapsToKeep_));
  //Update names if needed:
  if (tree_.get()) {
    sites->setSequencesNames(block.getSpeciesList(), true);
  }

  if (gapsAsUnresolved_)
    SiteContainerTools::changeGapsToUnknownCharacters(*sites);

  //Second we get the tree:
  const Tree* tree = 0;
  if (!tree_.get()) {
    //No default tree is given, we try to retrieve one from the block:
    if (!block.hasProperty(treePropertyIn_))
      throw Exception("MaximumLikelihoodModelFitMafIterator::fitModelBlock. No property available for " + treePropertyIn_);
    try {
      tree = &(dynamic_cast<const Tree&>(block.getProperty(treePropertyIn_)));
      if (rootFreqs_.get()) {
        if (!tree->isRooted())
          throw Exception("MaximumLikelihoodModelFitMafIterator::fitModelBlock. Tree must be rooted.");
      } else {
        if (tree->isRooted())
          throw Exception("MaximumLikelihoodModelFitMafIterator::fitModelBlock. Tree must be unrooted.");
      }
    } catch (bad_cast& e) {
      throw Exception("MaximumLikelihoodModelFitMafIterator::fitModelBlock. A property was found for '" + treePropertyIn_ + "' but does not appear to contain a phylogenetic tree.");
    }
  } else {
    tree = tree_.get();
  }

  //We build a new TreeLikelihood object:
  unique_ptr<DiscreteRatesAcrossSitesTreeLikelihood> tl;
  
  if (rootFreqs_.get()) {
    //Homogeneous, non-stationary
    modelSet_.reset(SubstitutionModelSetTools::createHomogeneousModelSet(model_->clone(), rootFreqs_->clone(), tree)); 
    tl.reset(new RNonHomogeneousTreeLikelihood(*tree, *sites, modelSet_.get(), rDist_.get(), false, true, false));
    //Initialize:
    if (initParameters_.size() == 0)
      init_(); //so far, even if tree changed, parameter names are supposingly the same. This might not be true in some complex cases...
  } else {
    if (modelSet_.get()) {
      //Non-homogeneous
      tl.reset(new RNonHomogeneousTreeLikelihood(*tree, *sites, modelSet_.get(), rDist_.get(), false, true, false));
    } else {
      //Homogeneous, stationary
      tl.reset(new RHomogeneousTreeLikelihood(*tree, *sites, model_.get(), rDist_.get(), false, false, true));
    }
  }
  tl->initialize();
  tl->setParameters(fixedParameters_);
  
  //We optimize parameters:
  ParameterList initParameters = initParameters_;
  if (reestimateBrLen_)
    initParameters.addParameters(tl->getBranchLengthsParameters());
  unsigned int nbIt = OptimizationTools::optimizeNumericalParameters2(tl.get(), initParameters, 0, 0.000001, 10000, 0, 0, reparametrize_, useClock_, 0);

  //And we save interesting parameter values:
  result_.setValue("NbIterations", static_cast<double>(nbIt));
  for (size_t i = 0;i < parametersOut_.size(); ++i) {
    result_.setValue(parametersOut_[i], tl->getParameterValue(parametersOut_[i]));
  }
}

void MaximumLikelihoodModelFitMafStatistics::init_() {
  if (modelSet_.get()) {
    initParameters_.addParameters(modelSet_->getIndependentParameters());
  } else {
    initParameters_.addParameters(model_->getIndependentParameters());
  }
  initParameters_.addParameters(rDist_->getIndependentParameters());
  //Remove from initParameters the ones to consider fixed:
  initParameters_.deleteParameters(fixedParameters_.getParameterNames());
 }

