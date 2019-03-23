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

add_library(aperture_particle_updater STATIC "particle_updater.cpp")
target_include_directories(aperture_particle_updater PUBLIC ${CMAKE_SOURCE_DIR}/include)
set( LINK_MPIPP TRUE )

#TODO precede particle/.cpp with lab/.cpp
function(module_gen category name)
  set(THIS_MODULE_NAME "${category}_${name}_module")
  add_library( ${THIS_MODULE_NAME} STATIC )
  target_include_directories(${THIS_MODULE_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/include)
  target_sources(${THIS_MODULE_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src/${category}/${name}.cpp)

  if( ${ARGC} GREATER 2 )
    if ( ${ARGV2} )
      target_link_libraries(${THIS_MODULE_NAME} PRIVATE mpipp_lib)
    endif()
  endif()

  if( "${catetory}" STREQUAL "particle" )
    target_link_libraries(aperture_particle_updater PRIVATE ${THIS_MODULE_NAME})
  endif()

endfunction()

module_gen(field field_shape_interplay)
module_gen(field current_deposition LINK_MPIPP)
module_gen(field communication LINK_MPIPP)

module_gen(particle migration LINK_MPIPP)
module_gen(particle pair_producer)
module_gen(particle pusher)
module_gen(particle array)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/old_field_solver)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ensemble)

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
# target_link_libraries(aperture PRIVATE ${old_field_solver} ${aperture_particle_updater})

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
