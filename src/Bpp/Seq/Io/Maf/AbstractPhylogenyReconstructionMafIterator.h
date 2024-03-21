// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _ABSTRACTPHYLOGENYRECONSTRUCTIONMAFITERATOR_H_
#define _ABSTRACTPHYLOGENYRECONSTRUCTIONMAFITERATOR_H_

#include <Bpp/Seq/Io/Maf/AbstractMafIterator.h>
#include <Bpp/Phyl/Tree/Tree.h>

namespace bpp
{
/**
 * @brief Partial implementation for phylogeny reconstruction iterator.
 *
 * This iterator calls a tree reconstruction method (to be implemented by the derivated class)
 * and store the resulting tree as an associated block property for the block,
 * before forwarding it.
 */
class AbstractPhylogenyReconstructionMafIterator :
  public AbstractFilterMafIterator
{
public:
  AbstractPhylogenyReconstructionMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator) :
    AbstractFilterMafIterator(iterator)
  {}

  virtual ~AbstractPhylogenyReconstructionMafIterator() {}

private:
  std::unique_ptr<MafBlock> analyseCurrentBlock_()
  {
    auto block = iterator_->nextBlock();
    if (!block) return nullptr;
    auto tree = buildTreeForBlock(*block);
    block->setProperty(getPropertyName(), std::move(tree));
    return block;
  }

public:
  virtual std::string getPropertyName() const = 0;
  virtual std::unique_ptr<Tree> buildTreeForBlock(const MafBlock& block) = 0;
};
} // end of namespace

#endif // _ABSTRACTPHYLOGENYRECONSTRUCTIONMAFITERATOR_H_
