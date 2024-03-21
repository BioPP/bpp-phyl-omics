// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _DISTANCEBASEDPHYLOGENYRECONSTRUCTIONMAFITERATOR_H_
#define _DISTANCEBASEDPHYLOGENYRECONSTRUCTIONMAFITERATOR_H_

#include "AbstractPhylogenyReconstructionMafIterator.h"
#include <Bpp/Phyl/Distance/DistanceMethod.h>

namespace bpp {

/**
 * @brief Implementation for distance-based phylogeny reconstruction iterator.
 */
class DistanceBasedPhylogenyReconstructionMafIterator:
  public AbstractPhylogenyReconstructionMafIterator
{
  private:
    std::string distanceProperty_;
    std::unique_ptr<DistanceMethodInterface> builder_;
  
  public:
    DistanceBasedPhylogenyReconstructionMafIterator(
        std::shared_ptr<MafIteratorInterface> iterator, 
	std::unique_ptr<DistanceMethodInterface> method,
       	const std::string& property):
      AbstractPhylogenyReconstructionMafIterator(iterator),
      distanceProperty_(property),
      builder_(std::move(method))
    {}

  private:
    DistanceBasedPhylogenyReconstructionMafIterator(const DistanceBasedPhylogenyReconstructionMafIterator& it):
      AbstractPhylogenyReconstructionMafIterator(0),
      distanceProperty_(it.distanceProperty_),
      builder_()
    {}

    DistanceBasedPhylogenyReconstructionMafIterator& operator=(const DistanceBasedPhylogenyReconstructionMafIterator& it)
    {
      distanceProperty_ = it.distanceProperty_;
      return *this;
    }

  public:
    void setDistanceProperty(const std::string& property) { distanceProperty_ = property; }
    const std::string& getDistanceProperty() const { return distanceProperty_; }

    std::string getPropertyName() const { return builder_->getName(); }
    std::unique_ptr<Tree> buildTreeForBlock(const MafBlock& block);
};

} //end of namespace bpp.

#endif //_DISTANCEBASEDPHYLOGENYRECONSTRUCTIONMAFITERATOR_H_
