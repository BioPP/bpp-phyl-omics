// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _COUNTCLUSTERSMAFSTATISTICS_H_
#define _COUNTCLUSTERSMAFSTATISTICS_H_

#include <Bpp/Seq/Io/Maf/MafStatistics.h>

//From bpp-phyl:
#include <Bpp/Phyl/Tree/Node.h>

//From the STL:
#include <map>

namespace bpp {

/**
 * @brief Count the number of sequence cluster, given a certain threshold.
 * This requires that a phylogenetic tree was previously computed.
 */
class CountClustersMafStatistics:
  public AbstractMafStatisticsSimple
{
  private:
    std::string treeProperty_;
    double threshold_;
  
  public:
    CountClustersMafStatistics(const std::string& property, double threshold):
      AbstractMafStatisticsSimple("NbClusters"),
      treeProperty_(property), threshold_(threshold)
    {}

  public:
    void setTreeProperty(const std::string& property) { treeProperty_ = property; }
    const std::string& getTreeProperty() const { return treeProperty_; }

    std::string getShortName() const { return "CountClusters(" + TextTools::toString(threshold_) + ")"; }
    std::string getFullName() const { return "Number of sequence clusters with divergence <= " + TextTools::toString(threshold_) + "."; }
    void compute(const MafBlock& block);

    std::string getPropertyName() const { return "CountClusters"; }

  private:
    unsigned int getNumberOfClusters_(const Node* node, std::map<const Node*, double>& heights);
};

} //end of namespace bpp.

#endif //_COUNTCLUSTERSMAFSTATISTICS_H_
