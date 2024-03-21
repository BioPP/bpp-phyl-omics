// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
      currentBlock_->setProperty(treePropertyWrite_, std::move(tree));
    } catch (bad_cast& e) {
      throw Exception("TreeManipulationMafIterator::analyseCurrentBlock_(). A property was found for '" + treePropertyRead_ + "' but does not appear to contain a phylogenetic tree.");
    }
  }
  return std::move(currentBlock_);
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

