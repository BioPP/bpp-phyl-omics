//
// File: TreeManipulationMafIterators.h
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

#ifndef _TREEMANIPULATIONMAFITERATOR_H_
#define _TREEMANIPULATIONMAFITERATOR_H_

#include <Bpp/Seq/Io/Maf/MafIterator.h>

//From bpp-phyl:
#include <Bpp/Phyl/TreeTemplate.h>

namespace bpp {

/**
 * @brief This iterator root associated trees according to an outgroup sequence.
 */
class TreeManipulationMafIterator:
  public AbstractFilterMafIterator
{
  private:
    std::string treePropertyRead_;
    std::string treePropertyWrite_;

  public:
    //Write can be the same as read.
    TreeManipulationMafIterator(MafIterator* iterator, const std::string& treePropertyRead, const std::string& treePropertyWrite) :
      AbstractFilterMafIterator(iterator), treePropertyRead_(treePropertyRead), treePropertyWrite_(treePropertyWrite)
    {}

  public:
    MafBlock* analyseCurrentBlock_();

  protected:
    virtual void manipulateTree_(TreeTemplate<Node>* tree) = 0;

};



/**
 * @brief This iterator root associated trees according to an outgroup sequence.
 */
class NewOutgroupMafIterator:
  public TreeManipulationMafIterator
{
  private:
    std::string outgroupSpecies_;

  public:
    //Write can be the same as read.
    NewOutgroupMafIterator(MafIterator* iterator, const std::string& treePropertyRead, const std::string& treePropertyWrite, const std::string& outgroupSpecies) :
      TreeManipulationMafIterator(iterator, treePropertyRead, treePropertyWrite), outgroupSpecies_(outgroupSpecies)
    {}

  private:
    void manipulateTree_(TreeTemplate<Node>* tree);

};



/**
 * @brief This iterator removes leaves of a certain species in an attached tree.
 */
class DropSpeciesMafIterator:
  public TreeManipulationMafIterator
{
  private:
    std::string species_;

  public:
    //Write can be the same as read.
    DropSpeciesMafIterator(MafIterator* iterator, const std::string& treePropertyRead, const std::string& treePropertyWrite, const std::string& species) :
      TreeManipulationMafIterator(iterator, treePropertyRead, treePropertyWrite), species_(species)
    {}

  private:
    void manipulateTree_(TreeTemplate<Node>* tree);

};


} //end of namespace bpp.

#endif //_TREEMANIPULATIONMAFITERATOR_H_
