//
// File: MaximumLikelihoodModelFitMafStatistics.h
// Created by: Julien Dutheil
// Created on: Mar 25 2014
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

#ifndef _MAXIMUMLIKELIHOODMODELFITMAFSTATISTICS_H_
#define _MAXIMUMLIKELIHOODMODELFITMAFSTATISTICS_H_

#include <Bpp/Seq/Io/Maf/MafStatistics.h>
#include <Bpp/Seq/Container/SiteContainer.h>

//From bpp-phyl:
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/FrequencySet/NucleotideFrequencySet.h>
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Likelihood/DiscreteRatesAcrossSitesTreeLikelihood.h>

//From bpp-core
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

namespace bpp {

/**
 * @brief Fit a substitution model.
 *
 * All nucleotide substitution models and rate distributions are supported.
 * Only time-homogeneous models are allowed though.
 */
class MaximumLikelihoodModelFitMafStatistics:
  public AbstractMafStatistics
{

  private:
    std::unique_ptr<SubstitutionModel> model_;
    std::unique_ptr<SubstitutionModelSet> modelSet_; //Only used in case of non-stationary model.
    std::unique_ptr<DiscreteDistribution> rDist_;
    std::unique_ptr<NucleotideFrequencySet> rootFreqs_;
    std::string treePropertyIn_;
    std::unique_ptr<const Tree> tree_;
    std::vector<std::string> parametersOut_;
    bool reestimateBrLen_;
    double propGapsToKeep_; //Exclude sites with too many gaps
    bool gapsAsUnresolved_;  //For most models, should be yes as they do not allow for gap characters
    bool useClock_;
    bool reparametrize_;
    ParameterList initParameters_;
    ParameterList fixedParameters_;

  public:
    /**
     * @brief Build a model fit maf statistic.
     *
     * A tree must be associated to each block before this analysis can be run.
     *
     * @param model The substitution model.
     * @param rDist The distribution of rates.
     * @param rootFreqs Root frequencies for non-stationary model. If set to 0, then a stationary model is assumed.
     * @param treePropertyIn The name of the property where the input tree is stored for each block.
     * @param parametersOut Parameters to output. 
     * @param fixedParameters Parameter which should not be estimated but fixed to the given value instead.
     * @param reestimateBrLen If the branch length from the tree should be reestimated (otherwise kept as is).
     * @param propGapsToKeep The maximum gapfrequency in a site to include it in the analysis. 
     * @param gapsAsUnresolved Tell if gap characters should be considered as unresolved states. In ost cases it should be set to true, as very few substitution models consider gaps as genuine states.
     * @param useClock Should a ultrametric tree be fitted (global clock)?
     * @param reparametrize Transform parameters to remove constraints.
     */
    MaximumLikelihoodModelFitMafStatistics(
        SubstitutionModel* model,
        DiscreteDistribution* rDist,
        NucleotideFrequencySet* rootFreqs,
        const std::string& treePropertyIn,
        const std::vector<std::string>& parametersOut,
        const ParameterList& fixedParameters,
        bool reestimateBrLen = true,
        double propGapsToKeep = 0,
        bool gapsAsUnresolved = true,
        bool useClock = false,
        bool reparametrize = false):
      AbstractMafStatistics(),
      model_(model), modelSet_(), rDist_(rDist), rootFreqs_(rootFreqs),
      treePropertyIn_(treePropertyIn), tree_(), parametersOut_(parametersOut),
      reestimateBrLen_(reestimateBrLen), propGapsToKeep_(propGapsToKeep), gapsAsUnresolved_(gapsAsUnresolved),
      useClock_(useClock), reparametrize_(reparametrize),
      initParameters_(), fixedParameters_(fixedParameters)
    {
      if (!rootFreqs) {
        init_();
        ApplicationTools::displayMessage("-- Available parameters:");
        std::vector<std::string> pl = model->getIndependentParameters().getParameterNames();
        for (size_t i = 0; i < pl.size(); ++i) {
          ApplicationTools::displayMessage("    " + pl[i]);
        }
        pl = rDist->getIndependentParameters().getParameterNames();
        for (size_t i = 0; i < pl.size(); ++i) {
          ApplicationTools::displayMessage("    " + pl[i]);
        }
      } else {
        //Otherwise we do not initialize parameters as the tree might change for each block.
        //We therefore have to initialize once for each block.
        std::vector<std::string> pl = model->getIndependentParameters().getParameterNames();
        for (size_t i = 0; i < pl.size(); ++i) {
          ApplicationTools::displayMessage("    " + pl[i] + "_1");
        }
        pl = rootFreqs->getIndependentParameters().getParameterNames();
        for (size_t i = 0; i < pl.size(); ++i) {
          ApplicationTools::displayMessage("    " + pl[i]);
        }
        pl = rDist->getIndependentParameters().getParameterNames();
        for (size_t i = 0; i < pl.size(); ++i) {
          ApplicationTools::displayMessage("    " + pl[i]);
        }
      }    
    }

    /**
     * @brief Build a model fit maf statistic.
     *
     * This analysis use the same input tree for all blocks.
     *
     * @param model The substitution model.
     * @param rDist The distribution of rates.
     * @param rootFreqs Root frequencies for non-stationary model. If set to 0, then a stationary model is assumed.
     * @param tree The tree to use for fitting the model.
     * @param parametersOut Parameters to output. 
     * @param fixedParameters Parameter which should not be estimated but fixed to the given value instead.
     * @param reestimateBrLen If the branch length from the tree should be reestimated (otherwise kept as is).
     * @param propGapsToKeep The maximum gapfrequency in a site to include it in the analysis. 
     * @param gapsAsUnresolved Tell if gap characters should be considered as unresolved states. In ost cases it should be set to true, as very few substitution models consider gaps as genuine states.
     * @param useClock Should a ultrametric tree be fitted (global clock)?
     * @param reparametrize Transform parameters to remove constraints.
     */
    MaximumLikelihoodModelFitMafStatistics(
        SubstitutionModel* model,
        DiscreteDistribution* rDist,
        NucleotideFrequencySet* rootFreqs,
        const Tree* tree,
        const std::vector<std::string>& parametersOut,
        const ParameterList& fixedParameters,
        bool reestimateBrLen = true,
        double propGapsToKeep = 0,
        bool gapsAsUnresolved = true,
        bool useClock = false,
        bool reparametrize = false):
      AbstractMafStatistics(),
      model_(model), modelSet_(), rDist_(rDist), rootFreqs_(rootFreqs),
      treePropertyIn_(NO_PROPERTY), tree_(tree), parametersOut_(parametersOut),
      reestimateBrLen_(reestimateBrLen), propGapsToKeep_(propGapsToKeep), gapsAsUnresolved_(gapsAsUnresolved),
      useClock_(useClock), reparametrize_(reparametrize),
      initParameters_(), fixedParameters_(fixedParameters)
    {
      if (!rootFreqs) {
        init_();
        ApplicationTools::displayMessage("-- Available parameters:");
        std::vector<std::string> pl = model->getIndependentParameters().getParameterNames();
        for (size_t i = 0; i < pl.size(); ++i) {
          ApplicationTools::displayMessage("    " + pl[i]);
        }
        pl = rDist->getIndependentParameters().getParameterNames();
        for (size_t i = 0; i < pl.size(); ++i) {
          ApplicationTools::displayMessage("    " + pl[i]);
        }
      } else{
        modelSet_.reset(SubstitutionModelSetTools::createHomogeneousModelSet(model->clone(), rootFreqs->clone(), tree));
        init_();
        ApplicationTools::displayMessage("-- Available parameters:");
        std::vector<std::string> pl = modelSet_->getIndependentParameters().getParameterNames();
        for (size_t i = 0; i < pl.size(); ++i) {
          ApplicationTools::displayMessage("    " + pl[i]);
        }
        pl = rDist->getIndependentParameters().getParameterNames();
        for (size_t i = 0; i < pl.size(); ++i) {
          ApplicationTools::displayMessage("    " + pl[i]);
        }
      }
    }

    /**
     * @brief Build a model fit maf statistic, using non-homogeneous models.
     *
     * This analysis use the same input tree for all blocks.
     *
     * @param modelSet The substitution model set.
     * @param rDist The distribution of rates.
     * @param tree The tree to use for fitting the model.
     * @param parametersOut Parameters to output. 
     * @param fixedParameters Parameter which should not be estimated but fixed to the given value instead.
     * @param reestimateBrLen If the branch length from the tree should be reestimated (otherwise kept as is).
     * @param propGapsToKeep The maximum gapfrequency in a site to include it in the analysis. 
     * @param gapsAsUnresolved Tell if gap characters should be considered as unresolved states. In ost cases it should be set to true, as very few substitution models consider gaps as genuine states.
     * @param useClock Should a ultrametric tree be fitted (global clock)?
     * @param reparametrize Transform parameters to remove constraints.
     */
    MaximumLikelihoodModelFitMafStatistics(
        SubstitutionModelSet* modelSet,
        DiscreteDistribution* rDist,
        const Tree* tree,
        const std::vector<std::string>& parametersOut,
        const ParameterList& fixedParameters,
        bool reestimateBrLen = true,
        double propGapsToKeep = 0,
        bool gapsAsUnresolved = true,
        bool useClock = false,
        bool reparametrize = false):
      AbstractMafStatistics(),
      model_(), modelSet_(modelSet), rDist_(rDist), rootFreqs_(),
      treePropertyIn_(NO_PROPERTY), tree_(tree), parametersOut_(parametersOut),
      reestimateBrLen_(reestimateBrLen), propGapsToKeep_(propGapsToKeep), gapsAsUnresolved_(gapsAsUnresolved),
      useClock_(useClock), reparametrize_(reparametrize),
      initParameters_(), fixedParameters_(fixedParameters)
    {
      init_();
      ApplicationTools::displayMessage("-- Available parameters:");
      std::vector<std::string> pl = modelSet_->getIndependentParameters().getParameterNames();
      for (size_t i = 0; i < pl.size(); ++i) {
        ApplicationTools::displayMessage("    " + pl[i]);
      }
      pl = rDist->getIndependentParameters().getParameterNames();
      for (size_t i = 0; i < pl.size(); ++i) {
        ApplicationTools::displayMessage("    " + pl[i]);
      }
    }

  private:
    MaximumLikelihoodModelFitMafStatistics(const MaximumLikelihoodModelFitMafStatistics& mafstat):
      AbstractMafStatistics(),
      model_(), modelSet_(), rDist_(), rootFreqs_(),
      treePropertyIn_(mafstat.treePropertyIn_), tree_(), parametersOut_(mafstat.parametersOut_),
      reestimateBrLen_(mafstat.reestimateBrLen_), propGapsToKeep_(mafstat.propGapsToKeep_), gapsAsUnresolved_(mafstat.gapsAsUnresolved_),
      useClock_(mafstat.useClock_), reparametrize_(mafstat.reparametrize_),
      initParameters_(mafstat.initParameters_), fixedParameters_(mafstat.fixedParameters_)
    {}
    
    MaximumLikelihoodModelFitMafStatistics& operator=(const MaximumLikelihoodModelFitMafStatistics& mafstat)
    {
      model_.reset();
      modelSet_.reset();
      rDist_.reset();
      rootFreqs_.reset();
      treePropertyIn_ = mafstat.treePropertyIn_;
      tree_.reset();
      parametersOut_ = mafstat.parametersOut_;
      reestimateBrLen_ = mafstat.reestimateBrLen_;
      propGapsToKeep_ = mafstat.propGapsToKeep_;
      gapsAsUnresolved_ = mafstat.gapsAsUnresolved_;
      initParameters_ = mafstat.initParameters_;
      fixedParameters_ = mafstat.fixedParameters_;
      useClock_ = mafstat.useClock_;
      reparametrize_ = mafstat.reparametrize_;
      return *this;
    }
     
  public:
    std::string getShortName() const { return "MLModelFit"; }
    std::string getFullName() const { return "Maximum Likelihood Model Fitting"; }
    void compute(const MafBlock& block);
    std::vector<std::string> getSupportedTags() const { 
      std::vector<std::string> tags;
      tags.push_back("NbIterations");
      tags.insert(tags.end(), parametersOut_.begin(), parametersOut_.end());
      return tags;
    }

    static const std::string NO_PROPERTY;
  
  private:
    void init_();
};

} //end of namespace bpp.

#endif //_MAXIMUMLIKELIHOODDISTANCEESTIMATIONMAFSTATISTICS_H_

