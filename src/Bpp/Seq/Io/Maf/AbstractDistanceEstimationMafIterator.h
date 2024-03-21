// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _ABSTRACTDISTANCEESTIMATIONMAFITERATOR_H_
#define _ABSTRACTDISTANCEESTIMATIONMAFITERATOR_H_

#include <Bpp/Seq/Io/Maf/AbstractMafIterator.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Distance/PGMA.h>
#include <Bpp/Phyl/Io/Newick.h>

namespace bpp
{
/**
 * @brief Partial implementation for distance estimation iterator.
 *
 * This iterator calls a distance reconstruction method (to be implemented by the derivated class)
 * and store the resulting distance matrix as an associated block property for the block,
 * before forwarding it.
 */
class AbstractDistanceEstimationMafIterator :
  public AbstractFilterMafIterator
{
private:
  bool addCoordinatesInSequenceNames_;

public:
  AbstractDistanceEstimationMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      bool addCoordinatesInSequenceNames) :
    AbstractFilterMafIterator(iterator),
    addCoordinatesInSequenceNames_(addCoordinatesInSequenceNames)
  {}

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_()
  {
    auto block = iterator_->nextBlock();
    if (!block) return nullptr;
    auto dist = estimateDistanceMatrixForBlock(*block);
    if (!addCoordinatesInSequenceNames_)
    {
      for (size_t i = 0; i < dist->size(); ++i)
      {
        std::string name = dist->getName(i);
        size_t pos = name.find('.');
        if (pos != std::string::npos)
        {
          dist->setName(i, name.substr(0, pos));
        }
      }
    }
    block->setProperty(getPropertyName(), std::move(dist));
    return block;
  }

public:
  virtual std::string getPropertyName() const = 0;
  virtual std::unique_ptr<DistanceMatrix> estimateDistanceMatrixForBlock(const MafBlock& block) = 0;
};
} // end of namespace bpp.

#endif // _ABSTRACTDISTANCEESTIMATIONMAFITERATOR_H_
