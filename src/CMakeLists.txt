# CMake script for Bio++ Phylogenetics Library
# Author: Sylvain Gaillard and Julien Dutheil
# Created: 08/11/2012

# File list
SET(CPP_FILES
  Bpp/Seq/Io/Maf/CountDistanceEstimationMafIterator.cpp
  Bpp/Seq/Io/Maf/DistanceBasedPhylogenyReconstructionMafIterator.cpp
  Bpp/Seq/Io/Maf/TreeManipulationMafIterators.cpp
  Bpp/Seq/Io/Maf/OutputTreeMafIterator.cpp
  Bpp/Seq/Io/Maf/CountClustersMafStatistics.cpp
  )

SET(H_FILES
  Bpp/Seq/Io/Maf/AbstractDistanceEstimationMafIterator.h
  Bpp/Seq/Io/Maf/CountDistanceEstimationMafIterator.h
  Bpp/Seq/Io/Maf/AbstractPhylogenyReconstructionMafIterator.h
  Bpp/Seq/Io/Maf/DistanceBasedPhylogenyReconstructionMafIterator.h
  Bpp/Seq/Io/Maf/TreeManipulationMafIterators.h
  Bpp/Seq/Io/Maf/OutputTreeMafIterator.h
  Bpp/Seq/Io/Maf/CountClustersMafStatistics.h
  )

# Build the static lib
ADD_LIBRARY(bppphylomics-static STATIC ${CPP_FILES})
SET_TARGET_PROPERTIES(bppphylomics-static
  PROPERTIES OUTPUT_NAME bpp-phyl-omics
  CLEAN_DIRECT_OUTPUT 1
  )
TARGET_LINK_LIBRARIES(bppphylomics-static ${LIBS})

# Build the shared lib
ADD_LIBRARY(bppphylomics-shared SHARED ${CPP_FILES})
SET_TARGET_PROPERTIES(bppphylomics-shared
  PROPERTIES OUTPUT_NAME bpp-phyl-omics
  CLEAN_DIRECT_OUTPUT 1
  VERSION ${BPPSEQOMICS_VERSION}
  SOVERSION ${BPPSEQOMICS_VERSION_MAJOR}
  )
TARGET_LINK_LIBRARIES(bppphylomics-shared ${LIBS})

# Install libs
INSTALL(TARGETS bppphylomics-static bppphylomics-shared DESTINATION lib${LIB_SUFFIX})

# Install headers
INSTALL(DIRECTORY Bpp/ DESTINATION include/Bpp FILES_MATCHING PATTERN "*.h" PATTERN "*.all")

INSTALL(CODE "EXECUTE_PROCESS(COMMAND python ${CMAKE_SOURCE_DIR}/genIncludes.py ${CMAKE_PREFIX_PATH}/include/Bpp/Seq)")

