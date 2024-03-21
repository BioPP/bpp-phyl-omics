// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "OutputTreeMafIterator.h"

using namespace bpp;
using namespace std;

void OutputTreeMafIterator::writeBlock_(std::ostream& out, const MafBlock& block) const
{
  //First get the tree for this block:
  if (!block.hasProperty(treeProperty_))
    throw Exception("OutputTreeMafIterator::writeBlock. No property available for " + treeProperty_);
  try {
    if (extendedSeqNames_) {
      const Tree& tree = dynamic_cast<const Tree&>(block.getProperty(treeProperty_));
      writer_.writeTree(tree, out);
    } else {
      TreeTemplate<Node> tree(dynamic_cast<const Tree&>(block.getProperty(treeProperty_)));
      stripNames_(*tree.getRootNode());
      writer_.writeTree(tree, out);
    }
  } catch (bad_cast& e) {
    throw Exception("OutputTreeMafIterator::writeBlock. A property was found for '" + treeProperty_ + "' but does not appear to contain a phylogenetic tree.");
  }
}

void OutputTreeMafIterator::stripNames_(Node& node) const
{
  if (node.hasName()) {
    string name = node.getName();
    size_t pos = name.find('.');
    if (pos != string::npos) {
      node.setName(name.substr(0, pos));
    }
  }
  for (size_t i = 0; i < node.getNumberOfSons(); ++i) {
    stripNames_(*node.getSon(i));
  }
}

