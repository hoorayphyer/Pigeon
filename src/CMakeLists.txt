include_directories( ../include )

# NOTE since all implementation cpp is separated from the instantiation cpp, just include the directory of src/ in Lab to achieve the automatic shadow that was impossible to achieve before

# if(LAB_ROOT_DIR)
#   include("${LAB_ROOT_DIR}/src/CMakeLists.txt")
#   include_directories(BEFORE "${LAB_ROOT_DIR}/include")
#   set(shadows "${LAB_SHADOWS}")
#   list(FIND shadows "UserControl.cpp" Lab_User_Control_Found)
#   if(Lab_User_Control_Found LESS 0)
#     set(Lab_User_Control_Found NO)
#   else()
#     set(Lab_User_Control_Found YES)
#   endif()
#   list(REMOVE_ITEM shadows "UserControl.cpp")
#   set(Lab_Aperture_src ${shadows})
#   message("${Yellow}Invoked LAB ${LAB_NAME} from ${LAB_ROOT_DIR}, with following shadowed${ColorReset}")
#   foreach(ARG ${LAB_SHADOWS})
#       message("${Yellow}    ${ARG}${ColorReset}")
#   endforeach()
# else()
#   message("${Yellow}No lab is invoked.${ColorReset}")
# endif()

#TODO precede particle/.cpp with lab/.cpp
macro(module_gen category name)
  set(module_target "${category}_${name}_module")
  if(EXISTS "aperture/${category}_${name}.cpp")
    set(src_file "aperture/${category}_${name}.cpp")
  else()
    set(src_file "${category}/${name}.cpp")
  endif()
  add_library( ${module_target} STATIC ${src_file} )
  target_include_directories(${module_target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  list( APPEND ${category}_modules ${module_target} )
  unset(src_file)
  unset(module_target)
endmacro()

add_library(mpipp_lib STATIC "parallel/mpi++.cpp" "parallel/mpi_communication.cpp" "parallel/mpi_datatype.cpp")
target_include_directories(mpipp_lib PUBLIC ${MPI_INCLUDE_PATH})
target_link_libraries(mpipp_lib PRIVATE ${MPI_LIBRARIES})
set_target_properties(mpipp_lib PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
set_target_properties(mpipp_lib PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")

function(link_mpipp target)
  # target_include_directories(${target} PUBLIC ${MPI_INCLUDE_PATH})
  target_link_libraries(${target} PRIVATE mpipp_lib)
  set_target_properties(${target} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  set_target_properties(${target} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endfunction()

macro(mpi_module_gen category name)
  module_gen( ${category} ${name} )
  link_mpipp( "${category}_${name}_module" )
endmacro()



add_library(silo_lib STATIC "io/silo++.cpp" "io/silo_optlist.cpp" "io/silo_operations.cpp")
target_include_directories(silo_lib PRIVATE ${Silo_INCLUDE_DIR} )
target_link_libraries( silo_lib PRIVATE ${Silo_LIBRARIES})
link_mpipp(silo_lib)

module_gen(field field_shape_interplay)
mpi_module_gen(field communication)

mpi_module_gen(particle migration)

module_gen(particle pair_producer)
module_gen(particle pusher)
module_gen(particle array)

add_library(particle_updater STATIC "particle/updater.cpp")
target_link_libraries(particle_updater PRIVATE ${particle_modules})

mpi_module_gen(aperture ensemble)

# TODO now src are scattered
# # modify Aperture_src using
# if(LAB_ROOT_DIR)
#   foreach(ARG ${Lab_Aperture_src})
#     list(REMOVE_ITEM Aperture_src ${ARG})
#     list(APPEND Aperture_src "${LAB_ROOT_DIR}/src/${ARG}")
#   endforeach()
# endif()

# if(LAB_ROOT_DIR AND Lab_User_Control_Found)
#   set(UserControl_src "${LAB_ROOT_DIR}/src/UserControl.cpp")
# else()
#   set(UserControl_src "UserControl.cpp")
# endif()
# list(APPEND UserControl_src "UserControlTemplates.cpp")
# add_library(UserControl STATIC ${UserControl_src})

include(../old_field_solver/CMakeLists.txt)

# add_executable(aperture "main.cpp" "aperture.cpp")
# target_link_libraries(aperture PRIVATE ${field_lib} ${particle_updater})
# link_mpipp(aperture)

# target_link_libraries(aperture AperLib UserControl)

# if(LAB_ROOT_DIR)
#   set_target_properties(aperture PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY ${LAB_ROOT_DIR})

#   add_custom_target( LABREVIEW ALL
#     ${CMAKE_COMMAND} -E cmake_echo_color --yellow "The following are shadowed by LAB ${LAB_NAME}")

#   #TODO: display shadowed headers also
#   foreach(ARG ${Lab_Aperture_src})
#     add_custom_command(TARGET LABREVIEW POST_BUILD
#       COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --yellow "    ${ARG}")
#   endforeach()
#   if(Lab_User_Control_Found)
#     add_custom_command(TARGET LABREVIEW POST_BUILD
#       COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --yellow "    UserControl.cpp")
#   endif()

#   add_dependencies( LABREVIEW aperture )
# else()
#   add_custom_target( LABREVIEW ALL
#     ${CMAKE_COMMAND} -E cmake_echo_color --yellow "NO lab is invoked")
# endif()
