// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "CountDistanceEstimationMafIterator.h"

#include <Bpp/Seq/Container/SiteContainerTools.h>

using namespace bpp;

std::unique_ptr<DistanceMatrix> CountDistanceEstimationMafIterator::estimateDistanceMatrixForBlock(const MafBlock& block)
{
  auto aln = block.getAlignment();
  auto dist = SiteContainerTools::computeSimilarityMatrix(*aln, true, gapOption_, unresolvedAsGap_);
  return dist;
}
