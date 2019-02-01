include_directories( ${MPI_INCLUDE_PATH} ${Silo_INCLUDE_DIR} ../include )

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
add_library(field_lib STATIC "field/field.cpp")
add_library(ptc_lib STATIC "particle/depositer.cpp" "particle/migration.cpp" "particle/pair_producer.cpp" "particle/particle_vector.cpp" "particle/pusher.cpp")


# TODO now src are scattered
# # modify Aperture_src using
# if(LAB_ROOT_DIR)
#   foreach(ARG ${Lab_Aperture_src})
#     list(REMOVE_ITEM Aperture_src ${ARG})
#     list(APPEND Aperture_src "${LAB_ROOT_DIR}/src/${ARG}")
#   endforeach()
# endif()

add_library(AperLib STATIC "aperture.cpp" )
target_link_libraries(AperLib PRIVATE ${field_lib} ${ptc_lib})

# if(MPI_COMPILE_FLAGS)
#   set_target_properties(AperLib PROPERTIES
#     COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
# endif()

# if(MPI_LINK_FLAGS)
#   set_target_properties(AperLib PROPERTIES
#     LINK_FLAGS "${MPI_LINK_FLAGS}")
# endif()

# if(LAB_ROOT_DIR AND Lab_User_Control_Found)
#   set(UserControl_src "${LAB_ROOT_DIR}/src/UserControl.cpp")
# else()
#   set(UserControl_src "UserControl.cpp")
# endif()
# list(APPEND UserControl_src "UserControlTemplates.cpp")
# add_library(UserControl STATIC ${UserControl_src})

add_executable(aperture "main.cpp")
# target_link_libraries(aperture AperLib UserControl)
target_link_libraries(aperture AperLib)

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
