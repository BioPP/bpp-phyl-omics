//
// File: OutputTreeMafIterator.h
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

#ifndef _OUTPUTTREEMAFITERATOR_H_
#define _OUTPUTTREEMAFITERATOR_H_

#include <Bpp/Seq/Io/Maf/MafIterator.h>

//From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
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

namespace bpp {

/**
 * @brief This iterator print an attached tree to a newick file.
 */
class OutputTreeMafIterator:
  public AbstractFilterMafIterator
{
  private:
    std::ostream* output_;
    std::string treeProperty_;
    Newick writer_;

  public:
    OutputTreeMafIterator(MafIterator* iterator, std::ostream* out, const std::string treeProperty) :
      AbstractFilterMafIterator(iterator), output_(out), treeProperty_(treeProperty), writer_()
    {}

  private:
    OutputTreeMafIterator(const OutputTreeMafIterator& iterator) :
      AbstractFilterMafIterator(0),
      output_(iterator.output_),
      treeProperty_(iterator.treeProperty_),
      writer_()
    {}
    
    OutputTreeMafIterator& operator=(const OutputTreeMafIterator& iterator)
    {
      output_ = iterator.output_;
      treeProperty_ = iterator.treeProperty_;
      writer_ = iterator.writer_;
      return *this;
    }


  public:
    MafBlock* analyseCurrentBlock_() throw (Exception) {
      currentBlock_ = iterator_->nextBlock();
      if (output_ && currentBlock_)
        writeBlock(*output_, *currentBlock_);
      return currentBlock_;
    }

  private:
    void writeBlock(std::ostream& out, const MafBlock& block) const;
};

} //end of namespace bpp.

#endif //_OUTPUTTREEMAFITERATOR_H_

