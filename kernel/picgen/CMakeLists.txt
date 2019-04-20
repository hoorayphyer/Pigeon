function(module_gen category name)
  set(THIS_MODULE_NAME "${category}_${name}")
  add_library( ${THIS_MODULE_NAME} SHARED ${category}_${name}.cpp )
  target_include_directories(${THIS_MODULE_NAME} PRIVATE ${KERNEL_DIR}/include)
  target_include_directories(${THIS_MODULE_NAME} PRIVATE ${KERNEL_DIR}/src)

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
target_include_directories(particle_updater PUBLIC ${KERNEL_DIR}/include )
target_link_libraries(particle_updater PRIVATE particle_forces particle_scattering particle_array field_current_deposition field_mesh_shape_interplay picgen)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/old_field_solver)

add_executable(pic "main.cpp")
target_include_directories(pic PRIVATE ${KERNEL_DIR} ${KERNEL_DIR}/include )
target_link_libraries(pic PRIVATE old_field_solver particle_updater particle_migration particle_forces particle_scattering dye mpipp picgen)

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
