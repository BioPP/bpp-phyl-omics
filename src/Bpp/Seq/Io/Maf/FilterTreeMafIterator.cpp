//
// File: FilterTreeMafIterator.cpp
// Created by: Julien Dutheil
// Created on: Mar 17 2017
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

#include "FilterTreeMafIterator.h"

using namespace bpp;
using namespace std;

MafBlock* FilterTreeMafIterator::analyseCurrentBlock_() throw (Exception)
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
          trashBuffer_.push_back(currentBlock_);
        currentBlock_ = iterator_->nextBlock();
      }
    } catch (bad_cast& e) {
      throw Exception("FilterTreeMafIterator::writeBlock. A property was found for '" + treeProperty_ + "' but does not appear to contain a phylogenetic tree.");
    }
  }
  return currentBlock_;
}


