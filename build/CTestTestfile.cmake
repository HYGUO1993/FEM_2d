# CMake generated Testfile for 
# Source directory: C:/Users/guoho/Documents/GitHub/FEM_2d
# Build directory: C:/Users/guoho/Documents/GitHub/FEM_2d/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
  add_test(basic_tests "C:/Users/guoho/Documents/GitHub/FEM_2d/build/bin/Debug/unit_tests.exe")
  set_tests_properties(basic_tests PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/guoho/Documents/GitHub/FEM_2d/CMakeLists.txt;19;add_test;C:/Users/guoho/Documents/GitHub/FEM_2d/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
  add_test(basic_tests "C:/Users/guoho/Documents/GitHub/FEM_2d/build/bin/Release/unit_tests.exe")
  set_tests_properties(basic_tests PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/guoho/Documents/GitHub/FEM_2d/CMakeLists.txt;19;add_test;C:/Users/guoho/Documents/GitHub/FEM_2d/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
  add_test(basic_tests "C:/Users/guoho/Documents/GitHub/FEM_2d/build/bin/MinSizeRel/unit_tests.exe")
  set_tests_properties(basic_tests PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/guoho/Documents/GitHub/FEM_2d/CMakeLists.txt;19;add_test;C:/Users/guoho/Documents/GitHub/FEM_2d/CMakeLists.txt;0;")
elseif(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
  add_test(basic_tests "C:/Users/guoho/Documents/GitHub/FEM_2d/build/bin/RelWithDebInfo/unit_tests.exe")
  set_tests_properties(basic_tests PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/guoho/Documents/GitHub/FEM_2d/CMakeLists.txt;19;add_test;C:/Users/guoho/Documents/GitHub/FEM_2d/CMakeLists.txt;0;")
else()
  add_test(basic_tests NOT_AVAILABLE)
endif()
