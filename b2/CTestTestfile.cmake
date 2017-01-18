# CMake generated Testfile for 
# Source directory: /home/victor/repos/nbots
# Build directory: /home/victor/repos/nbots/b2
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(unit_tests "utests")
set_tests_properties(unit_tests PROPERTIES  ENVIRONMENT "NB_FONTS_DIR=/home/victor/repos/nbots/resources/fonts")
