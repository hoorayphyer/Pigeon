include_directories( ${Silo_INCLUDE_DIR} ../include )

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
  add_library( ${category}_${name}_module STATIC "aperture/${category}_${name}.cpp")
  target_include_directories(${category}_${name}_module PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  list( APPEND ${category}_modules ${category}_${name}_module )
endmacro()

function(link_mpi target)
  target_include_directories(${target} PUBLIC ${MPI_INCLUDE_PATH})
  target_link_libraries(${target} PRIVATE ${MPI_LIBRARIES})
  set_target_properties(${target} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  set_target_properties(${target} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endfunction()

set( particle_modules "" )

add_library(parallel_lib STATIC "parallel/locale.cpp" "parallel/mpi++.cpp" "parallel/mpi_communication.cpp")
target_include_directories(parallel_lib PUBLIC ${MPI_INCLUDE_PATH})
link_mpi(parallel_lib)


# add_library(field_lib STATIC "field/field.cpp" )


# module_gen(particle migration)
# link_mpi(particle_migration_module)
# target_link_libraries(particle_migration_module PRIVATE ${parallel_lib})

# module_gen(particle depositer)
# module_gen(particle pair_producer)
module_gen(particle pusher)
module_gen(particle array)

# add_library(particle_update_lib STATIC "particle/updater.cpp")
# target_link_libraries(particle_update_lib PRIVATE ${particle_modules})


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

# add_executable(aperture "main.cpp" "aperture.cpp")
# target_link_libraries(aperture PRIVATE ${field_lib} ${particle_update_lib})

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
