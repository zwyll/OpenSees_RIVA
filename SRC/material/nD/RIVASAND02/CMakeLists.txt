#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
    PRIVATE
      RIVASAND02.cpp
    PUBLIC
      RIVASAND02.h
      RIVASAND02Kernel.h
      RIVASAND02BaseKernel.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})
