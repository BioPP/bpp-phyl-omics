// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "MaximumLikelihoodDistanceEstimationMafIterator.h"

// From bpp-seq:
#include <Bpp/Seq/Container/SiteContainerTools.h>

using namespace bpp;

unique_ptr<DistanceMatrix> MaximumLikelihoodDistanceEstimationMafIterator::estimateDistanceMatrixForBlock(const MafBlock& block)
{
  // First we get the alignment:
  auto sites = block.getAlignment();
  SiteContainerTools::removeGapSites(*sites, propGapsToKeep_);
  if (gapsAsUnresolved_)
  {
    SiteContainerTools::changeGapsToUnknownCharacters(*sites);
  }

  // Check that there is enough data:
  if (sites->getNumberOfSites() == 0) {
    vector<string> names = sites->getSequenceNames();
    auto mat = std::make_unique<DistanceMatrix>(names);
    return mat;
  } else {
    // Set the data and fit the matrix:
    distEst_->setData(std::move(sites));
    ParameterList p;
    auto mat = OptimizationTools::estimateDistanceMatrix(*distEst_, p, paramOpt_, verbose_);
    return mat;
  }
}
