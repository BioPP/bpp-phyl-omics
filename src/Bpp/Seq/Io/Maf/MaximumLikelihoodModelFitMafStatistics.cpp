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
#include <Bpp/Phyl/Tree/PhyloTreeTools.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>

using namespace bpp;
using namespace std;

const string MaximumLikelihoodModelFitMafStatistics::NO_PROPERTY = "RESERVED_NOPROPERTY";

void MaximumLikelihoodModelFitMafStatistics::compute(const MafBlock& block)
{
  //First we get the alignment:
  shared_ptr<SiteContainerInterface> sites = block.getAlignment();
  SiteContainerTools::removeGapSites(*sites, propGapsToKeep_);
  //Update names if needed:
  if (tree_.get()) {
    sites->setSequenceNames(block.getSpeciesList(), true);
  }

  if (gapsAsUnresolved_)
    SiteContainerTools::changeGapsToUnknownCharacters(*sites);

  //Second we get the tree:
  const TreeTemplate<Node>* tree = 0;
  if (!tree_.get()) {
    //No default tree is given, we try to retrieve one from the block:
    if (!block.hasProperty(treePropertyIn_))
      throw Exception("MaximumLikelihoodModelFitMafStatistics::fitModelBlock. No property available for " + treePropertyIn_);
    try {
      tree = &(dynamic_cast<const TreeTemplate<Node>&>(block.getProperty(treePropertyIn_)));
      if (process_->hasRootFrequencySet()) {
        if (!tree->isRooted())
          throw Exception("MaximumLikelihoodModelFitMafStatistics::fitModelBlock. Tree must be rooted.");
      } else {
        if (tree->isRooted())
          throw Exception("MaximumLikelihoodModelFitMafStatistics::fitModelBlock. Tree must be unrooted.");
      }
    } catch (bad_cast& e) {
      throw Exception("MaximumLikelihoodModelFitMafStatistics::fitModelBlock. A property was found for '" + treePropertyIn_ + "' but does not appear to contain a phylogenetic tree.");
    }
  } else {
    tree = tree_.get();
  }
  //Convert to PhyloTree:
  shared_ptr<PhyloTree> phyloTree = PhyloTreeTools::buildFromTreeTemplate(*tree);
  process_->setPhyloTree(*phyloTree);
  

  //We build a new TreeLikelihood object:
  Context context;
  auto lik = make_shared<LikelihoodCalculationSingleProcess>(context, sites, process_);
  auto treeLik = make_shared<SingleProcessPhyloLikelihood>(context, lik);
  
  treeLik->setParameters(fixedParameters_);
  
  //We optimize parameters:
  ParameterList initParameters = initParameters_;
  if (reestimateBrLen_)
    initParameters.addParameters(treeLik->getBranchLengthParameters());
  unsigned int nbIt = OptimizationTools::optimizeNumericalParameters2(treeLik, initParameters, 0, 0.000001, 10000, 0, 0, reparametrize_, useClock_, 0);

  //And we save interesting parameter values:
  result_.setValue("NbIterations", static_cast<double>(nbIt));
  for (size_t i = 0;i < parametersOut_.size(); ++i) {
    result_.setValue(parametersOut_[i], treeLik->getParameterValue(parametersOut_[i]));
  }
}

void MaximumLikelihoodModelFitMafStatistics::init_() {
  initParameters_.addParameters(process_->getIndependentParameters());
  //Remove from initParameters the ones to consider fixed:
  initParameters_.deleteParameters(fixedParameters_.getParameterNames());
  ApplicationTools::displayMessage("-- Available parameters:");
  std::vector<std::string> pl = process_->getIndependentParameters().getParameterNames();
  for (auto p: pl) {
    ApplicationTools::displayMessage("    " + p);
  }
}

