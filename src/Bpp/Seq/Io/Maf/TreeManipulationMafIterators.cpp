//
// File: TreeManipulationMafIterators.cpp
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

#include "TreeManipulationMafIterators.h"

using namespace bpp;
using namespace std;

unique_ptr<MafBlock> TreeManipulationMafIterator::analyseCurrentBlock_()
{
  currentBlock_ = iterator_->nextBlock();
  if (currentBlock_) {
    if (!currentBlock_->hasProperty(treePropertyRead_))
      throw Exception("TreeManipulationMafIterator::analyseCurrentBlock_(). No property available for " + treePropertyRead_);
    try {
      auto tree = make_unique<TreeTemplate<Node>>(dynamic_cast<const Tree&>(currentBlock_->getProperty(treePropertyRead_)));
      manipulateTree_(*tree);
      currentBlock_->setProperty(treePropertyWrite_, move(tree));
    } catch (bad_cast& e) {
      throw Exception("TreeManipulationMafIterator::analyseCurrentBlock_(). A property was found for '" + treePropertyRead_ + "' but does not appear to contain a phylogenetic tree.");
    }
  }
  return move(currentBlock_);
}


void NewOutgroupMafIterator::manipulateTree_(TreeTemplate<Node>& tree)
{
  vector<Node*> leaves = tree.getLeaves();
  Node* outgroup = 0;
  bool outgroupFound = false;
  for (size_t i = 0; i < leaves.size() && !outgroupFound; ++i) {
    string species, chr;
    MafSequence::splitNameIntoSpeciesAndChromosome(leaves[i]->getName(), species, chr);
    if (species == outgroupSpecies_) {
      outgroup = leaves[i];
      outgroupFound = true;
    }
  }
  if (!outgroupFound)
    throw Exception("NewOutgroupTreeMafIterator::analyseCurrentBlock_(). No ougroup species was found in the attached tree.");
  tree.newOutGroup(outgroup);
}


void DropSpeciesMafIterator::manipulateTree_(TreeTemplate<Node>& tree)
{
  vector<Node*> leaves = tree.getLeaves();
  for (size_t i = 0; i < leaves.size(); ++i) {
    string species, chr;
    MafSequence::splitNameIntoSpeciesAndChromosome(leaves[i]->getName(), species, chr);
    if (species == species_) {
      TreeTemplateTools::dropSubtree(tree, leaves[i]);
    }
  }
}

