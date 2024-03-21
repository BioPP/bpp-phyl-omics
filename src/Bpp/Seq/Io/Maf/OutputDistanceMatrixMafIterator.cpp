// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "OutputDistanceMatrixMafIterator.h"

using namespace bpp;
using namespace std;

void OutputDistanceMatrixMafIterator::writeBlock_(std::ostream& out, const MafBlock& block) const
{
  //First get the tree for this block:
  if (!block.hasProperty(distProperty_))
    throw Exception("OutputDistanceMatrixMafIterator::writeBlock. No property available for " + distProperty_);
  try {
    if (extendedSeqNames_) {
      const DistanceMatrix& mat = dynamic_cast<const DistanceMatrix&>(block.getProperty(distProperty_));
      writer_.writeDistanceMatrix(mat, out);
    } else {
      DistanceMatrix mat(dynamic_cast<const DistanceMatrix&>(block.getProperty(distProperty_)));
      vector<string> names = mat.getNames();
      stripNames_(names);
      mat.setNames(names);
      writer_.writeDistanceMatrix(mat, out);
    }
  } catch (bad_cast& e) {
    throw Exception("OutputDistanceMatrixMafIterator::writeBlock. A property was found for '" + distProperty_ + "' but does not appear to contain a distance matrix.");
  }

  out << endl << endl;
}

void OutputDistanceMatrixMafIterator::stripNames_(vector<string>& names) const
{
  for (vector<string>::iterator it = names.begin(); it != names.end(); ++it) {
    size_t pos = it->find('.');
    if (pos != string::npos) {
      *it = it->substr(0, pos);
    }
  }
}

