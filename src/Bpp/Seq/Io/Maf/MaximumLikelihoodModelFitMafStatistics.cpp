// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
  const TreeTemplate<Node>* tree = nullptr;
  if (!tree_) {
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

