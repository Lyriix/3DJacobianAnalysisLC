#CMakeLists.txt to build project using cmake
#run $cmake CMakeLists.txt

cmake_minimum_required(VERSION 2.6)
project(pgm)

set(QT_USE_OPENGL TRUE)
set(ITK_DIR "/home/charly/softwares/ITK/ITK_BIN_DIR/bin")

FIND_PACKAGE(Qt4 REQUIRED)
FIND_PACKAGE(OpenGL REQUIRED)
FIND_PACKAGE(ITK REQUIRED)


INCLUDE(${QT_USE_FILE})
INCLUDE(${ITK_USE_FILE})

MESSAGE(${ITK_USE_FILE})
MESSAGE(${ITK_LIBRARIES})

include_directories(${CMAKE_CURRENT_BINARY_DIR})

LINK_DIRECTORIES(/usr/local/lib)

file(
GLOB_RECURSE
source_files
src/*.[cht]pp
shaders/*.frag
shaders/*.vert
)

SET(CMAKE_BUILD_TYPE Debug)
ADD_DEFINITIONS( -Wall -Wextra -Wfloat-equal -Wunreachable-code -std=c++11)


SET(UI src/local/interface/mainwindow.ui)
SET(MOC src/lib/interface/application_qt.hpp
        src/local/interface/myWindow.hpp
        src/local/interface/myWidgetGL.hpp)

QT4_WRAP_CPP(MOC_GENERATED ${MOC})
QT4_WRAP_UI(UI_GENERATED ${UI})



add_executable(
  pgm
  ${source_files}
  ${UI_GENERATED}
  ${MOC_GENERATED}
)


TARGET_LINK_LIBRARIES(pgm -lm -ldl -lGLEW ${OPENGL_LIBRARIES} ${QT_LIBRARIES} ${QT_GL_LIBRARIES} ${QT_QTOPENGL_LIBRARY} ${ITK_LIBRARIES})


