# CMake script for Bio++ Phylo(gen)-omics Library
# Authors:
#   Sylvain Gaillard
#   Julien Dutheil
#   Francois Gindraud (2017)
# Created: 08/11/2012

set (CPP_FILES
  Bpp/Seq/Io/Maf/CountClustersMafStatistics.cpp
  Bpp/Seq/Io/Maf/CountDistanceEstimationMafIterator.cpp
  Bpp/Seq/Io/Maf/DistanceBasedPhylogenyReconstructionMafIterator.cpp
  Bpp/Seq/Io/Maf/MaximumLikelihoodDistanceEstimationMafIterator.cpp
  Bpp/Seq/Io/Maf/MaximumLikelihoodModelFitMafStatistics.cpp
  Bpp/Seq/Io/Maf/DistanceBasedPhylogenyReconstructionMafIterator.cpp
  Bpp/Seq/Io/Maf/TreeManipulationMafIterators.cpp
  Bpp/Seq/Io/Maf/FilterTreeMafIterator.cpp
  Bpp/Seq/Io/Maf/OutputTreeMafIterator.cpp
  Bpp/Seq/Io/Maf/OutputDistanceMatrixMafIterator.cpp
  Bpp/Seq/Io/Maf/CountClustersMafStatistics.cpp
  )

# Build the static lib
add_library (${PROJECT_NAME}-static STATIC ${CPP_FILES})
target_include_directories (${PROJECT_NAME}-static PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )
set_target_properties (${PROJECT_NAME}-static PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
target_link_libraries (${PROJECT_NAME}-static ${BPP_LIBS_STATIC})
target_compile_options (${PROJECT_NAME}-static
  PUBLIC ${public-compile-options}
  PRIVATE ${private-compile-options}
  )

# Build the shared lib
add_library (${PROJECT_NAME}-shared SHARED ${CPP_FILES})
target_include_directories (${PROJECT_NAME}-shared PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )
set_target_properties (${PROJECT_NAME}-shared
  PROPERTIES OUTPUT_NAME ${PROJECT_NAME}
  MACOSX_RPATH 1
  VERSION ${${PROJECT_NAME}_VERSION}
  SOVERSION ${${PROJECT_NAME}_VERSION_MAJOR}
  )
target_link_libraries (${PROJECT_NAME}-shared ${BPP_LIBS_SHARED})
target_compile_options (${PROJECT_NAME}-shared
  PUBLIC ${public-compile-options}
  PRIVATE ${private-compile-options}
  )

# Install libs and headers
install (
  TARGETS ${PROJECT_NAME}-static ${PROJECT_NAME}-shared
  EXPORT ${PROJECT_NAME}-targets
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX}
  )
install (
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
  DESTINATION include
  FILES_MATCHING PATTERN "*.h"
  )
# Install cmake file representing targets
install (EXPORT ${PROJECT_NAME}-targets DESTINATION ${cmake-package-location})
