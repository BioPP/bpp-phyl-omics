//
// File: FilterTreeMafIterator.h
// Created by: Julien Dutheil
// Created on: Mar 17 2017
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

#ifndef _FILTERTREEMAFITERATOR_H_
#define _FILTERTREEMAFITERATOR_H_

#include <Bpp/Seq/Io/Maf/MafIterator.h>

//From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Io/Newick.h>

namespace bpp {

/**
 * @brief This iterator filter blocks based on the attached tree.
 *
 * Currently, filtering is based on the branch lengths. Block with a tree with one branch length exceeding a threshold are removed.
 */
class FilterTreeMafIterator:
  public AbstractFilterMafIterator,
  public virtual MafTrashIterator
{
  private:
    std::string treeProperty_;
    double maxBrLen_;
    std::deque<MafBlock*> trashBuffer_;
    bool keepTrashedBlocks_;

  public:
    FilterTreeMafIterator(MafIterator* iterator, const std::string& treeProperty, double maxBrLen, bool keepTrashedBlocks) :
      AbstractFilterMafIterator(iterator), treeProperty_(treeProperty), maxBrLen_(maxBrLen), trashBuffer_(), keepTrashedBlocks_(keepTrashedBlocks)
    {}

  private:
    FilterTreeMafIterator(const FilterTreeMafIterator& iterator) :
      AbstractFilterMafIterator(0),
      treeProperty_(iterator.treeProperty_),
      maxBrLen_(iterator.maxBrLen_),
      trashBuffer_(iterator.trashBuffer_),
      keepTrashedBlocks_(iterator.keepTrashedBlocks_)
    {}
    
    FilterTreeMafIterator& operator=(const FilterTreeMafIterator& iterator)
    {
      treeProperty_      = iterator.treeProperty_;
      maxBrLen_          = iterator.maxBrLen_;
      trashBuffer_       = iterator.trashBuffer_;
      keepTrashedBlocks_ = iterator.keepTrashedBlocks_;
      return *this;
    }


  public:
    MafBlock* analyseCurrentBlock_();

    MafBlock* nextRemovedBlock() {
      if (trashBuffer_.size() == 0) return 0;
      MafBlock* block = trashBuffer_.front();
      trashBuffer_.pop_front();
      return block;
    }

};

} //end of namespace bpp.

#endif //_FILTERTREEMAFITERATOR_H_

