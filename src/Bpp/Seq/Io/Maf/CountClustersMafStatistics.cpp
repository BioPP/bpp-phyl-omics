// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "CountClustersMafStatistics.h"

// From bpp-seq:
#include <Bpp/Seq/Container/SiteContainerTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree/TreeTemplate.h>

using namespace bpp;
using namespace std;

unsigned int CountClustersMafStatistics::getNumberOfClusters_(const Node* node, map<const Node*, double>& heights)
{
  unsigned int nClust = 0;
  double h = heights[node];
  if (h < threshold_)
  {
    nClust++;
  }
  else
  {
    for (int i = 0; i < static_cast<int>(node->getNumberOfSons()); ++i)
    {
      nClust += getNumberOfClusters_((*node)[i], heights);
    }
  }
  return nClust;
}

void CountClustersMafStatistics::compute(const MafBlock& block)
{
  if (!block.hasProperty(treeProperty_))
    throw Exception("CountClustersMafStatistics::compute. No property available for " + treeProperty_);
  try
  {
    TreeTemplate<Node> tree(dynamic_cast<const Tree&>(block.getProperty(treeProperty_)));
    if (!tree.isRooted())
      throw Exception("CountClustersMafStatistics::compute. Cluster count only works with a rooted tree.");
    // Compute all tree heights:
    map<const Node*, double> heights;
    TreeTemplateTools::getHeights(*tree.getRootNode(), heights);
    unsigned int nClust = getNumberOfClusters_(tree.getRootNode(), heights);
    result_.setValue("NbClusters", nClust);
  }
  catch (bad_cast& e)
  {
    throw Exception("CountClustersMafStatistics::compute. A property was found for '" + treeProperty_ + "' but does not appear to contain a phylogenetic tree.");
  }
}
