# Install script for directory: /home/longxing/Rifdock/rifdock/apps/rosetta

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_librosetta" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_librosetta")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_librosetta"
         RPATH "/usr/local/lib:/usr/local/lib64:/home/longxing/rosetta-master/source/cmake/build_cxx11_omp")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/longxing/Rifdock/rifdock/build/apps/rosetta/test_librosetta")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_librosetta" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_librosetta")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_librosetta"
         OLD_RPATH "/usr/local/lib:/usr/local/lib64:/home/longxing/rosetta-master/source/cmake/build_cxx11_omp:/home/longxing/Rifdock/rifdock/build/apps/rosetta/riflib:"
         NEW_RPATH "/usr/local/lib:/usr/local/lib64:/home/longxing/rosetta-master/source/cmake/build_cxx11_omp")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_librosetta")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/rifgen" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/rifgen")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/rifgen"
         RPATH "/usr/local/lib:/usr/local/lib64:/home/longxing/rosetta-master/source/cmake/build_cxx11_omp")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/longxing/Rifdock/rifdock/build/apps/rosetta/rifgen")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/rifgen" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/rifgen")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/rifgen"
         OLD_RPATH "/usr/local/lib:/usr/local/lib64:/home/longxing/rosetta-master/source/cmake/build_cxx11_omp:/home/longxing/Rifdock/rifdock/build/apps/rosetta/riflib:"
         NEW_RPATH "/usr/local/lib:/usr/local/lib64:/home/longxing/rosetta-master/source/cmake/build_cxx11_omp")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/rifgen")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/rif_dock_test" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/rif_dock_test")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/rif_dock_test"
         RPATH "/usr/local/lib:/usr/local/lib64:/home/longxing/rosetta-master/source/cmake/build_cxx11_omp")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/longxing/Rifdock/rifdock/build/apps/rosetta/rif_dock_test")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/rif_dock_test" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/rif_dock_test")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/rif_dock_test"
         OLD_RPATH "/usr/local/lib:/usr/local/lib64:/home/longxing/rosetta-master/source/cmake/build_cxx11_omp:/home/longxing/Rifdock/rifdock/build/apps/rosetta/riflib:"
         NEW_RPATH "/usr/local/lib:/usr/local/lib64:/home/longxing/rosetta-master/source/cmake/build_cxx11_omp")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/rif_dock_test")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/scheme_make_bounding_grids" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/scheme_make_bounding_grids")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/scheme_make_bounding_grids"
         RPATH "/usr/local/lib:/usr/local/lib64:/home/longxing/rosetta-master/source/cmake/build_cxx11_omp")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/longxing/Rifdock/rifdock/build/apps/rosetta/scheme_make_bounding_grids")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/scheme_make_bounding_grids" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/scheme_make_bounding_grids")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/scheme_make_bounding_grids"
         OLD_RPATH "/usr/local/lib:/usr/local/lib64:/home/longxing/rosetta-master/source/cmake/build_cxx11_omp:/home/longxing/Rifdock/rifdock/build/apps/rosetta/riflib:"
         NEW_RPATH "/usr/local/lib:/usr/local/lib64:/home/longxing/rosetta-master/source/cmake/build_cxx11_omp")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/scheme_make_bounding_grids")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/longxing/Rifdock/rifdock/build/apps/rosetta/python/cmake_install.cmake")
  include("/home/longxing/Rifdock/rifdock/build/apps/rosetta/riflib/cmake_install.cmake")

endif()

