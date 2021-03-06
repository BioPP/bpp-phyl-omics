//
// File: DistanceBasedPhylogenyReconstructionMafIterators.h
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
    std::unique_ptr<DistanceMethod> builder_;
  
  public:
    DistanceBasedPhylogenyReconstructionMafIterator(MafIterator* iterator, DistanceMethod* method, const std::string& property):
      AbstractPhylogenyReconstructionMafIterator(iterator),
      distanceProperty_(property), builder_(method)
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
    Tree* buildTreeForBlock(const MafBlock& block);
};

} //end of namespace bpp.

#endif //_DISTANCEBASEDPHYLOGENYRECONSTRUCTIONMAFITERATOR_H_
