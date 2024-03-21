// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _OUTPUTTREEMAFITERATOR_H_
#define _OUTPUTTREEMAFITERATOR_H_

#include <Bpp/Seq/Io/Maf/AbstractMafIterator.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/Io/Newick.h>

/**
 * @mainpage
 *
 * @par
 * The bpp-phyl-omics library contains 'omics' phylogenetic tools and classes.
 *
 * As for Bio++ version 2.1.0, these consist of bpp::MafIterator classes allowing to reconstruct phylogenies from a genome alignment.
 * More tools are expected to developped in the incoming versions.
 */

namespace bpp
{
/**
 * @brief This iterator print an attached tree to a newick file.
 */
class OutputTreeMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::shared_ptr<std::ostream> output_;
  std::string treeProperty_;
  Newick writer_;
  bool extendedSeqNames_;

public:
  OutputTreeMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      std::shared_ptr<std::ostream> out,
      const std::string& treeProperty,
      bool extendedSeqNames = true) :
    AbstractFilterMafIterator(iterator),
    output_(out),
    treeProperty_(treeProperty),
    writer_(),
    extendedSeqNames_(extendedSeqNames)
  {}

private:
  OutputTreeMafIterator(const OutputTreeMafIterator& iterator) :
    AbstractFilterMafIterator(0),
    output_(iterator.output_),
    treeProperty_(iterator.treeProperty_),
    writer_(),
    extendedSeqNames_(iterator.extendedSeqNames_)
  {}

  OutputTreeMafIterator& operator=(const OutputTreeMafIterator& iterator)
  {
    output_ = iterator.output_;
    treeProperty_ = iterator.treeProperty_;
    writer_ = iterator.writer_;
    extendedSeqNames_ = iterator.extendedSeqNames_;
    return *this;
  }

public:
  std::unique_ptr<MafBlock> analyseCurrentBlock_()
  {
    currentBlock_ = iterator_->nextBlock();
    if (output_ && currentBlock_)
      writeBlock_(*output_, *currentBlock_);
    return std::move(currentBlock_);
  }

private:
  void writeBlock_(std::ostream& out, const MafBlock& block) const;
  void stripNames_(Node& node) const;
};
} // end of namespace bpp.

#endif // _OUTPUTTREEMAFITERATOR_H_
