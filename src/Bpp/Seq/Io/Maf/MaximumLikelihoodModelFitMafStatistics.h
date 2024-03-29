// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _MAXIMUMLIKELIHOODMODELFITMAFSTATISTICS_H_
#define _MAXIMUMLIKELIHOODMODELFITMAFSTATISTICS_H_

#include <Bpp/Seq/Io/Maf/MafStatistics.h>
#include <Bpp/Seq/Container/SiteContainer.h>

// From bpp-phyl:
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Phyl/Model/FrequencySet/NucleotideFrequencySet.h>
#include <Bpp/Phyl/Likelihood/AutonomousSubstitutionProcess.h>
#include <Bpp/Phyl/Tree/Tree.h>

// From bpp-core
// #include <Bpp/Numeric/Prob/DiscreteDistribution.h>

namespace bpp
{
/**
 * @brief Fit a substitution model.
 *
 * All nucleotide substitution models and rate distributions are supported.
 * Only time-homogeneous models are allowed though.
 */
class MaximumLikelihoodModelFitMafStatistics :
  public AbstractMafStatistics
{
private:
  std::shared_ptr<AutonomousSubstitutionProcessInterface> process_;
  std::string treePropertyIn_;
  std::shared_ptr<const TreeTemplate<Node>> tree_;
  std::vector<std::string> parametersOut_;
  bool reestimateBrLen_;
  double propGapsToKeep_; // Exclude sites with too many gaps
  bool gapsAsUnresolved_;  // For most models, should be yes as they do not allow for gap characters
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
   * @param process The substitution process to use.
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
      std::shared_ptr<AutonomousSubstitutionProcessInterface> process,
      const std::string& treePropertyIn,
      const std::vector<std::string>& parametersOut,
      const ParameterList& fixedParameters,
      bool reestimateBrLen = true,
      double propGapsToKeep = 0,
      bool gapsAsUnresolved = true,
      bool useClock = false,
      bool reparametrize = false) :
    AbstractMafStatistics(),
    process_(process),
    treePropertyIn_(treePropertyIn), tree_(), parametersOut_(parametersOut),
    reestimateBrLen_(reestimateBrLen), propGapsToKeep_(propGapsToKeep), gapsAsUnresolved_(gapsAsUnresolved),
    useClock_(useClock), reparametrize_(reparametrize),
    initParameters_(), fixedParameters_(fixedParameters)
  {
    init_();
  }

  /**
   * @brief Build a model fit maf statistic.
   *
   * This analysis use the same input tree for all blocks.
   *
   * @param process The substitution process.
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
      std::shared_ptr<AutonomousSubstitutionProcessInterface> process,
      std::shared_ptr<const TreeTemplate<Node>> tree,
      const std::vector<std::string>& parametersOut,
      const ParameterList& fixedParameters,
      bool reestimateBrLen = true,
      double propGapsToKeep = 0,
      bool gapsAsUnresolved = true,
      bool useClock = false,
      bool reparametrize = false) :
    AbstractMafStatistics(),
    process_(process),
    treePropertyIn_(NO_PROPERTY), tree_(tree), parametersOut_(parametersOut),
    reestimateBrLen_(reestimateBrLen), propGapsToKeep_(propGapsToKeep), gapsAsUnresolved_(gapsAsUnresolved),
    useClock_(useClock), reparametrize_(reparametrize),
    initParameters_(), fixedParameters_(fixedParameters)
  {
    init_();
  }

private:
  MaximumLikelihoodModelFitMafStatistics(const MaximumLikelihoodModelFitMafStatistics& mafstat) :
    AbstractMafStatistics(),
    process_(),
    treePropertyIn_(mafstat.treePropertyIn_), tree_(), parametersOut_(mafstat.parametersOut_),
    reestimateBrLen_(mafstat.reestimateBrLen_), propGapsToKeep_(mafstat.propGapsToKeep_), gapsAsUnresolved_(mafstat.gapsAsUnresolved_),
    useClock_(mafstat.useClock_), reparametrize_(mafstat.reparametrize_),
    initParameters_(mafstat.initParameters_), fixedParameters_(mafstat.fixedParameters_)
  {}

  MaximumLikelihoodModelFitMafStatistics& operator=(const MaximumLikelihoodModelFitMafStatistics& mafstat)
  {
    process_.reset();
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
  std::vector<std::string> getSupportedTags() const
  {
    std::vector<std::string> tags;
    tags.push_back("NbIterations");
    tags.insert(tags.end(), parametersOut_.begin(), parametersOut_.end());
    return tags;
  }

  static const std::string NO_PROPERTY;

private:
  void init_();
};
} // end of namespace bpp.

#endif // _MAXIMUMLIKELIHOODMODELFITMAFSTATISTICS_H_
