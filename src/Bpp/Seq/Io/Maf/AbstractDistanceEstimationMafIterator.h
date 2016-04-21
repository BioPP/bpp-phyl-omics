//
// File: AbstractDistanceEstimationMafIterator.h
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

#ifndef _ABSTRACTDISTANCEESTIMATIONMAFITERATOR_H_
#define _ABSTRACTDISTANCEESTIMATIONMAFITERATOR_H_

#include <Bpp/Seq/Io/Maf/MafIterator.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Distance/PGMA.h>
#include <Bpp/Phyl/Io/Newick.h>

namespace bpp  {

/**
 * @brief Partial implementation for distance estimation iterator.
 *
 * This iterator calls a distance reconstruction method (to be implemented by the derivated class)
 * and store the resulting distance matrix as an associated block property for the block,
 * before forwarding it.
 */
class AbstractDistanceEstimationMafIterator:
  public AbstractFilterMafIterator
{
  private:
    bool addCoordinatesInSequenceNames_;

  public:
    AbstractDistanceEstimationMafIterator(MafIterator* iterator, bool addCoordinatesInSequenceNames):
      AbstractFilterMafIterator(iterator),
      addCoordinatesInSequenceNames_(addCoordinatesInSequenceNames)
    {}

  private:
    MafBlock* analyseCurrentBlock_() throw (Exception)
    {
      MafBlock* block = iterator_->nextBlock();
      if (!block) return 0;
      unique_ptr<DistanceMatrix> dist(estimateDistanceMatrixForBlock(*block));
      if (!addCoordinatesInSequenceNames_) {
        for (size_t i = 0; i < dist->size(); ++i) {
          std::string name = dist->getName(i);
          size_t pos = name.find('.');
          if (pos != std::string::npos) {
            dist->setName(i, name.substr(0, pos));
          }
        }
      }
      block->setProperty(getPropertyName(), dist.release());
      return block;
    }

  public:
    virtual std::string getPropertyName() const = 0;
    virtual DistanceMatrix* estimateDistanceMatrixForBlock(const MafBlock& block) = 0;

};

} //end of namespace bpp.

#endif // _ABSTRACTDISTANCEESTIMATIONMAFITERATOR_H_
