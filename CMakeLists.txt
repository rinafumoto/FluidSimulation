cmake_minimum_required(VERSION 3.12)
#-------------------------------------------------------------------------------------------
# I'm going to use vcpk in most cases for our install of 3rd party libs
# this is going to check the environment variable for CMAKE_TOOLCHAIN_FILE and this must point to where
# vcpkg.cmake is in the University this is set in your .bash_profile to
# export CMAKE_TOOLCHAIN_FILE=/public/devel/2020/vcpkg/scripts/buildsystems/vcpkg.cmake
#-------------------------------------------------------------------------------------------
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE AND DEFINED ENV{CMAKE_TOOLCHAIN_FILE})
   set(CMAKE_TOOLCHAIN_FILE $ENV{CMAKE_TOOLCHAIN_FILE})
endif()
# Name of the project
project(FluidSimulationBuild)
# This is the name of the Exe change this and it will change everywhere
set(TargetName FluidSimulation)
# This will include the file NGLConfig.cmake, you need to add the location to this either using
# -DCMAKE_PREFIX_PATH=~/NGL or as a system environment variable.
find_package(NGL CONFIG REQUIRED)
# Instruct CMake to run moc automatically when needed (Qt projects only)
set(CMAKE_AUTOMOC ON)
# Run the uic tools and search in the ui folder
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOUIC_SEARCH_PATHS ${PROJECT_SOURCE_DIR}/ui)
# find Qt libs first we check for Version 6
find_package(Qt6 COMPONENTS OpenGL Widgets OpenGLWidgets QUIET )
if ( Qt6_FOUND )
    message("Found Qt6 Using that")
else()
    message("Found Qt5 Using that")
    find_package(Qt5 COMPONENTS OpenGL Widgets REQUIRED)
endif()
# use C++ 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
# Set the name of the executable we want to build
add_executable(${TargetName})

target_sources(${TargetName} PRIVATE ${PROJECT_SOURCE_DIR}/src/main.cpp
${PROJECT_SOURCE_DIR}/src/MainWindow.cpp
${PROJECT_SOURCE_DIR}/src/NGLScene.cpp
${PROJECT_SOURCE_DIR}/src/NGLSceneMouseControls.cpp
${PROJECT_SOURCE_DIR}/src/Fluid.cpp
${PROJECT_SOURCE_DIR}/include/MainWindow.h
${PROJECT_SOURCE_DIR}/include/Fluid.h
${PROJECT_SOURCE_DIR}/include/NGLScene.h
${PROJECT_SOURCE_DIR}/include/WindowParams.h
)
if ( Qt6_FOUND )
  target_link_libraries(${TargetName} PRIVATE  NGL Qt::Widgets Qt::OpenGL Qt::OpenGLWidgets)
else()
  target_link_libraries(${TargetName} PRIVATE  NGL Qt::Widgets Qt::OpenGL )
endif()

add_custom_target(CopyShadersAndfonts ALL
   COMMAND ${CMAKE_COMMAND} -E copy_directory
   ${CMAKE_CURRENT_SOURCE_DIR}/shaders
   $<TARGET_FILE_DIR:${TargetName}>/shaders
)
