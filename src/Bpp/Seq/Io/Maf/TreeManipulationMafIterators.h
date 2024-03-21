// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _TREEMANIPULATIONMAFITERATOR_H_
#define _TREEMANIPULATIONMAFITERATOR_H_

#include <Bpp/Seq/Io/Maf/AbstractMafIterator.h>

//From bpp-phyl:
#include <Bpp/Phyl/Tree/TreeTemplate.h>

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
    TreeManipulationMafIterator(
	std::shared_ptr<MafIteratorInterface> iterator,
       	const std::string& treePropertyRead,
       	const std::string& treePropertyWrite) :
      AbstractFilterMafIterator(iterator), 
      treePropertyRead_(treePropertyRead),
      treePropertyWrite_(treePropertyWrite)
    {}

  public:
    std::unique_ptr<MafBlock> analyseCurrentBlock_();

  protected:
    virtual void manipulateTree_(TreeTemplate<Node>& tree) = 0;

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
    NewOutgroupMafIterator(
        std::shared_ptr<MafIteratorInterface> iterator,
       	const std::string& treePropertyRead,
       	const std::string& treePropertyWrite,
       	const std::string& outgroupSpecies) :
      TreeManipulationMafIterator(
          iterator,
	  treePropertyRead, 
	  treePropertyWrite),
      outgroupSpecies_(outgroupSpecies)
    {}

  private:
    void manipulateTree_(TreeTemplate<Node>& tree) override;

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
    DropSpeciesMafIterator(
	std::shared_ptr<MafIteratorInterface> iterator,
       	const std::string& treePropertyRead,
       	const std::string& treePropertyWrite,
       	const std::string& species) :
      TreeManipulationMafIterator(
	  iterator, 
	  treePropertyRead,
	  treePropertyWrite),
       	species_(species)
    {}

  private:
    void manipulateTree_(TreeTemplate<Node>& tree) override;

};


} //end of namespace bpp.

#endif //_TREEMANIPULATIONMAFITERATOR_H_
