// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "DistanceBasedPhylogenyReconstructionMafIterator.h"

#include <Bpp/Seq/DistanceMatrix.h>

using namespace bpp;
using namespace std;

std::unique_ptr<Tree> DistanceBasedPhylogenyReconstructionMafIterator::buildTreeForBlock(const MafBlock& block)
{
  // First get the distance matrix for this block:
  if (!block.hasProperty(distanceProperty_))
    throw Exception("DistanceBasedPhylogenyReconstructionMafIterator::buildTreeForBlock. No property available for " + distanceProperty_);
  try
  {
    const DistanceMatrix& dist = dynamic_cast<const DistanceMatrix&>(block.getProperty(distanceProperty_));
    builder_->setDistanceMatrix(dist);
    builder_->computeTree();
    auto& tree = builder_->tree();
    return std::unique_ptr<Tree>(tree.clone());
  }
  catch (bad_cast& e)
  {
    throw Exception("DistanceBasedPhylogenyReconstructionMafIterator::buildTreeForBlock. A property was found for '" + distanceProperty_ + "' but does not appear to contain a distance matrix.");
  }
  catch (Exception& e)
  {
    throw Exception("DistanceBasedPhylogenyReconstructionMafIterator::buildTreeForBlock. Tree reconstruction failed!");
  }
}
