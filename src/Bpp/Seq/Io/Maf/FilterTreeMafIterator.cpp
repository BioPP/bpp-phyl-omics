// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "FilterTreeMafIterator.h"

using namespace bpp;
using namespace std;

unique_ptr<MafBlock> FilterTreeMafIterator::analyseCurrentBlock_()
{
  bool test = false;
  currentBlock_ = iterator_->nextBlock();
  while (!test && currentBlock_) {
    //First get the tree for this block:
    if (!currentBlock_->hasProperty(treeProperty_))
      throw Exception("FilterTreeMafIterator::writeBlock. No property available for " + treeProperty_);
    try {
      TreeTemplate<Node> tree(dynamic_cast<const Tree&>(currentBlock_->getProperty(treeProperty_)));
      vector<double> brlen = tree.getBranchLengths();
      test = true;
      for (double l : brlen) {
        if (l > maxBrLen_) {
          test = false;
          break;
        }     
      }
      if (!test) {
        if (keepTrashedBlocks_)
          trashBuffer_.push_back(std::move(currentBlock_));
        currentBlock_ = iterator_->nextBlock();
      }
    } catch (bad_cast& e) {
      throw Exception("FilterTreeMafIterator::writeBlock. A property was found for '" + treeProperty_ + "' but does not appear to contain a phylogenetic tree.");
    }
  }
  return std::move(currentBlock_);
}


