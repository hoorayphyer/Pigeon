# NOTE since all implementation cpp is separated from the instantiation cpp, just include the directory of src/ in Lab to achieve the automatic shadow that was impossible to achieve before

if(NOT TARGET picgen)
  include(${CMAKE_SOURCE_DIR}/pulsar/CMakeLists.txt) #TODO use the actual lab root
endif()

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
function(module_gen category name)
  set(THIS_MODULE_NAME "${category}_${name}")
  add_library( ${THIS_MODULE_NAME} SHARED ${CMAKE_SOURCE_DIR}/picgen/${category}_${name}.cpp )
  target_include_directories(${THIS_MODULE_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/include)
  target_include_directories(${THIS_MODULE_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src)

  target_link_libraries(${THIS_MODULE_NAME} PRIVATE picgen)
  foreach( lib IN LISTS ARGN)
    target_link_libraries(${THIS_MODULE_NAME} PRIVATE ${lib})
  endforeach()
endfunction()

module_gen(field communication mpipp)
module_gen(field mesh_shape_interplay)
module_gen(field current_deposition)

module_gen(particle migration mpipp)
module_gen(particle scattering)
module_gen(particle forces)
module_gen(particle array)

add_library(particle_updater SHARED "particle_updater.cpp")
target_include_directories(particle_updater PUBLIC ${CMAKE_SOURCE_DIR}/include )
target_link_libraries(particle_updater PRIVATE particle_forces particle_scattering particle_array field_current_deposition field_mesh_shape_interplay picgen)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/old_field_solver)

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
# add_library(UserControl SHARED ${UserControl_src})

add_executable(pic "main.cpp")
target_include_directories(pic PRIVATE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include )
target_link_libraries(pic PRIVATE old_field_solver particle_updater particle_migration particle_forces particle_scattering dye mpipp picgen)

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
