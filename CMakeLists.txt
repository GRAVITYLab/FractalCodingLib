cmake_minimum_required( VERSION 2.8 )

project( FractalCodingLib )

# The required packages
#find_package(CUDA)

# CUDA
#set(CUDA_DIR "${CUDA_TOOLKIT_ROOT_DIR}")
#set(CUDA_INC_DIR "${CUDA_DIR}/include")
#set(CUDA_LIB_DIR "${CUDA_DIR}/lib64")

# CUDA SDK
#set(CUDASDK_DIR "${CUDA_SDK_ROOT_DIR}")
#set(CUDASDK_INC_DIR "${CUDASDK_DIR}/common/inc")
#set(CUDASDK_LIB_DIR "${CUDASDK_DIR}/lib")

# Set ITL path
set(ITL_DIR "/home/abon/code/chaudhua/install" CACHE PATH "The path to ITL")

# Set OSUFlow path
set( OSUFLOW_INSTALL_DIR "/home/abon/install/osuflow/osuflow" CACHE PATH "Installed path of osuflow" )

# Set Matlab engine
#set( MATENGINE_INC_DIR "/home/abon/install/matlabr2011b/extern/include" CACHE PATH "Path to Matlab engine" )

# Set path to my own libraries
set( MY_INSTALL_DIR "/home/abon/code/chaudhua/install" CACHE PATH "Installed path of my libraries" )

# set OpenCV Path
#set( OPENCV_DIR "/homes/chaudhua/install" CACHE PATH "The installed path of OpenCV" )
#set( OPENCV_INC "${OPENCV_DIR}/include" )
#set( OPENCV_INC2 "${OPENCV_DIR}/include/opencv" )
#set( OPENCV_INC3 "${OPENCV_DIR}/include/opencv2" )
#set( OPENCV_LIB "${OPENCV_DIR}/lib")

# Set VTK path
#set(VTK_DIR "../.." CACHE PATH "The path to VTK")

#####
#include_directories ( "${CUDA_INC_DIR}"
#		      		  "${CUDASDK_INC_DIR}" )
#####
include_directories ( "${OSUFLOW_INSTALL_DIR}/include" )
#####
include_directories( "${ITL_DIR}/include" )
#####
#include_directories ( "${MATENGINE_INC_DIR}" )
#####
#include_directories( 
#		${OPENCV_INC}
#		${OPENCV_INC2}
#		${OPENCV_INC3}
#		)
#####
#include_directories ( 
#		"${VTK_DIR}/include/vtk-5.10"
#	 )
#####
#include_directories ( "${MY_INSTALL_DIR}/FieldPartitionLib/include" )
#####
include_directories ( "${PROJECT_SOURCE_DIR}/include" )
#cuda_include_directories ( "${PROJECT_SOURCE_DIR}/include" )

#cuda_add_library( FractalCodingLib STATIC 
add_library( FractalCodingLib STATIC 
                  ${PROJECT_SOURCE_DIR}/src/FEL_domain_core.cpp
				  ${PROJECT_SOURCE_DIR}/src/FEL_spandecomposition.cpp
                  ${PROJECT_SOURCE_DIR}/src/FEL_encoder_incremental.cpp
                  #${PROJECT_SOURCE_DIR}/src/FEL_encoder_localsearch.cpp
                  #${PROJECT_SOURCE_DIR}/src/FEL_encoder_dct.cpp
                  ${PROJECT_SOURCE_DIR}/src/FEL_decoder.cpp
				  #${PROJECT_SOURCE_DIR}/src/FEL_decoder_gpu.cu
		  		  #${PROJECT_SOURCE_DIR}/src/FEL_core_gpu.cu
		 		  ${PROJECT_SOURCE_DIR}/src/FEL_codebookentry.cpp
		  		  ${PROJECT_SOURCE_DIR}/src/FEL_encodedhistogram.cpp
		  		  ${PROJECT_SOURCE_DIR}/src/FEL_encodedsimplehistogram.cpp
                  ${PROJECT_SOURCE_DIR}/src/FEL_segmenttreenode.cpp
		  		  ${PROJECT_SOURCE_DIR}/src/FEL_segmenttree2.cpp
                  #${PROJECT_SOURCE_DIR}/src/FEL_encoder.cpp
                  ${PROJECT_SOURCE_DIR}/src/FEL_transformation.cpp
                  #${PROJECT_SOURCE_DIR}/src/FEL_matlab.cpp
 )
     
# Installation Instructions
# lib is interpreted as relative to CMAKE_INSTALL_PREFIX
FILE(GLOB headers "${PROJECT_SOURCE_DIR}/include/*.h")
install( FILES ${headers} DESTINATION include )
install( TARGETS FractalCodingLib DESTINATION lib )
