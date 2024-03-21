// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _FILTERTREEMAFITERATOR_H_
#define _FILTERTREEMAFITERATOR_H_

#include <Bpp/Seq/Io/Maf/AbstractMafIterator.h>

//From bpp-phyl:
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/Io/Newick.h>

namespace bpp {

/**
 * @brief This iterator filter blocks based on the attached tree.
 *
 * Currently, filtering is based on the branch lengths. Block with a tree with one branch length exceeding a threshold are removed.
 */
class FilterTreeMafIterator:
  public AbstractFilterMafIterator,
  public virtual MafTrashIteratorInterface
{
  private:
    std::string treeProperty_;
    double maxBrLen_;
    std::deque<std::unique_ptr<MafBlock>> trashBuffer_;
    bool keepTrashedBlocks_;

  public:
    FilterTreeMafIterator(
	std::shared_ptr<MafIteratorInterface> iterator,
       	const std::string& treeProperty,
       	double maxBrLen,
       	bool keepTrashedBlocks) :
      AbstractFilterMafIterator(iterator),
      treeProperty_(treeProperty),
      maxBrLen_(maxBrLen),
      trashBuffer_(),
      keepTrashedBlocks_(keepTrashedBlocks)
    {}

  private:
    FilterTreeMafIterator(const FilterTreeMafIterator& iterator) :
      AbstractFilterMafIterator(0),
      treeProperty_(iterator.treeProperty_),
      maxBrLen_(iterator.maxBrLen_),
      trashBuffer_(),
      keepTrashedBlocks_(iterator.keepTrashedBlocks_)
    {}
    
    FilterTreeMafIterator& operator=(const FilterTreeMafIterator& iterator)
    {
      treeProperty_      = iterator.treeProperty_;
      maxBrLen_          = iterator.maxBrLen_;
      trashBuffer_.clear();
      keepTrashedBlocks_ = iterator.keepTrashedBlocks_;
      return *this;
    }


  public:
    std::unique_ptr<MafBlock> analyseCurrentBlock_();

    std::unique_ptr<MafBlock> nextRemovedBlock() {
      if (trashBuffer_.size() == 0) return nullptr;
      auto block = std::move(trashBuffer_.front());
      trashBuffer_.pop_front();
      return block;
    }

};

} //end of namespace bpp.

#endif //_FILTERTREEMAFITERATOR_H_

