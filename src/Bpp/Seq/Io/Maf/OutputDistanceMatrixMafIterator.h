// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _OUTPUTDISTANCEMATRIXMAFITERATOR_H_
#define _OUTPUTDISTANCEMATRIXMAFITERATOR_H_

#include <Bpp/Seq/Io/Maf/AbstractMafIterator.h>

// From bpp-seq:
#include <Bpp/Seq/DistanceMatrix.h>
#include <Bpp/Seq/Io/PhylipDistanceMatrixFormat.h>

namespace bpp
{
/**
 * @brief This iterator print an attached distance matrix in phylip format.
 */
class OutputDistanceMatrixMafIterator :
  public AbstractFilterMafIterator
{
private:
  std::shared_ptr<std::ostream> output_;
  std::string distProperty_;
  PhylipDistanceMatrixFormat writer_;
  bool extendedSeqNames_;

public:
  OutputDistanceMatrixMafIterator(
      std::shared_ptr<MafIteratorInterface> iterator,
      std::shared_ptr<std::ostream> out,
      const std::string& distProperty,
      bool extendedSeqNames = true) :
    AbstractFilterMafIterator(iterator),
    output_(out),
    distProperty_(distProperty),
    writer_(),
    extendedSeqNames_(extendedSeqNames)
  {}

private:
  OutputDistanceMatrixMafIterator(const OutputDistanceMatrixMafIterator& iterator) :
    AbstractFilterMafIterator(0),
    output_(iterator.output_),
    distProperty_(iterator.distProperty_),
    writer_(),
    extendedSeqNames_(iterator.extendedSeqNames_)
  {}

  OutputDistanceMatrixMafIterator& operator=(const OutputDistanceMatrixMafIterator& iterator)
  {
    output_ = iterator.output_;
    distProperty_ = iterator.distProperty_;
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
  void stripNames_(std::vector<std::string>& names) const;
};
} // end of namespace bpp.

#endif // _OUTPUTDISTANCEMATRIXMAFITERATOR_H_
