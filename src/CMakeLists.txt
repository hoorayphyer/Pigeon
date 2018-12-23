include_directories(${MPI_INCLUDE_PATH} ${GSL_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${Silo_INCLUDE_DIR} ../include )

if (CUDA_FOUND)
  add_definitions(-DCUDA_ENABLED)
  include_directories(${CUDA_INCLUDE_DIRS})
  message("Cuda library found, enabling the GPU part of the code.")
else()
  message("Cuda library not found, skipping the GPU part of the code.")
endif()

if(LAB_ROOT_DIR)
  include("${LAB_ROOT_DIR}/src/CMakeLists.txt")
  include_directories(BEFORE "${LAB_ROOT_DIR}/include")
  set(shadows "${LAB_SHADOWS}")
  list(FIND shadows "UserControl.cpp" Lab_User_Control_Found)
  if(Lab_User_Control_Found LESS 0)
    set(Lab_User_Control_Found NO)
  else()
    set(Lab_User_Control_Found YES)
  endif()
  list(REMOVE_ITEM shadows "UserControl.cpp")
  set(Lab_Aperture_src ${shadows})
  message("${Yellow}Invoked LAB ${LAB_NAME} from ${LAB_ROOT_DIR}, with following shadowed${ColorReset}")
  foreach(ARG ${LAB_SHADOWS})
      message("${Yellow}    ${ARG}${ColorReset}")
  endforeach()
else()
  message("${Yellow}No lab is invoked.${ColorReset}")
endif()

set(Aperture_src "Fields.cpp" "DataExporter.cpp" "mpi_wrapper.cpp" "Domain.cpp" "FiniteDiff.cpp" "AperParams.cpp" "FieldUpdater.cpp" "GoogleTestAuxiliaries.cpp" "Particles.cpp" "ParticleSorter.cpp" "ParticlePusher.cpp" "CurrentDepositer.cpp" "ParticleDepositer.cpp" "Aperture.cpp" "InitialConditions.cpp" "ParticleUpdater.cpp" "Logger.cpp" "ParticleInjector.cpp" "ParticleCommunicator.cpp" "InfoCollector.cpp" "AperData_IO.cpp" "PairProducer.cpp" "Rebalancer.cpp" "Rng.cpp" "DynamicAdjustCoordinator.cpp" "ProcessRoster.cpp" "InfoOutputter.cpp" "EnsParticleGatherer.cpp" "AnnihManager.cpp" "FileSystem.cpp")

# modify Aperture_src using
if(LAB_ROOT_DIR)
  foreach(ARG ${Lab_Aperture_src})
    list(REMOVE_ITEM Aperture_src ${ARG})
    list(APPEND Aperture_src "${LAB_ROOT_DIR}/src/${ARG}")
  endforeach()
endif()

add_library(AperLib STATIC ${Aperture_src})
target_link_libraries(AperLib ${MPI_LIBRARIES} ${Boost_LIBRARIES} ${Silo_LIBRARIES} ${GSL_LIBRARIES})

if (CUDA_FOUND)
  target_link_libraries(AperLib ${CUDA_LIBRARIES})
endif()

if(MPI_COMPILE_FLAGS)
  set_target_properties(AperLib PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(AperLib PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

add_subdirectory(cuda)

if(LAB_ROOT_DIR AND Lab_User_Control_Found)
  set(UserControl_src "${LAB_ROOT_DIR}/src/UserControl.cpp")
else()
  set(UserControl_src "UserControl.cpp")
endif()
list(APPEND UserControl_src "UserControlTemplates.cpp")
add_library(UserControl STATIC ${UserControl_src})

add_executable(aperture "main.cpp")
target_link_libraries(aperture AperLib UserControl)

if(LAB_ROOT_DIR)
  set_target_properties(aperture PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${LAB_ROOT_DIR})

  add_custom_target( LABREVIEW ALL
    ${CMAKE_COMMAND} -E cmake_echo_color --yellow "The following are shadowed by LAB ${LAB_NAME}")

  #TODO: display shadowed headers also
  foreach(ARG ${Lab_Aperture_src})
    add_custom_command(TARGET LABREVIEW POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --yellow "    ${ARG}")
  endforeach()
  if(Lab_User_Control_Found)
    add_custom_command(TARGET LABREVIEW POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --yellow "    UserControl.cpp")
  endif()

  add_dependencies( LABREVIEW aperture )
else()
  add_custom_target( LABREVIEW ALL
    ${CMAKE_COMMAND} -E cmake_echo_color --yellow "NO lab is invoked")
endif()

if (CUDA_FOUND)
  target_link_libraries(aperture Pulsar_dev)
endif()
