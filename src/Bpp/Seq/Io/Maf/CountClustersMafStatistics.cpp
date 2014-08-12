//
// File: MafIterators.cpp
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

#include "CountClustersMafStatistics.h"

//From bpp-seq:
#include <Bpp/Seq/Container/SiteContainerTools.h>

//From bpp-phyl:
#include <Bpp/Phyl/TreeTemplate.h>

using namespace bpp;
using namespace std;

unsigned int CountClustersMafStatistics::getNumberOfClusters_(const Node* node, map<const Node*, double>& heights)
{
  unsigned int nClust = 0;
  double h = heights[node];
  if (h < threshold_) {
    nClust++;
  } else {
    for (int i = 0; i < static_cast<int>(node->getNumberOfSons()); ++i) {
      nClust += getNumberOfClusters_((*node)[i], heights);
    }
  }
  return nClust;
}

void CountClustersMafStatistics::compute(const MafBlock& block)
{
  if (!block.hasProperty(treeProperty_))
    throw Exception("CountClustersMafStatistics::compute. No property available for " + treeProperty_);
  try {
    TreeTemplate<Node> tree(dynamic_cast<const Tree&>(block.getProperty(treeProperty_)));
    if (!tree.isRooted())
      throw Exception("CountClustersMafStatistics::compute. Cluster count only works with a rooted tree.");
    //Compute all tree heights:
    map<const Node*, double> heights;
    TreeTemplateTools::getHeights(*tree.getRootNode(), heights);
    unsigned int nClust = getNumberOfClusters_(tree.getRootNode(), heights);
    result_.setValue("NbClusters", nClust);
  } catch (bad_cast& e) {
    throw Exception("CountClustersMafStatistics::compute. A property was found for '" + treeProperty_ + "' but does not appear to contain a phylogenetic tree.");
  }
}

