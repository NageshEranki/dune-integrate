add_executable("dune-integrate" dune-integrate.cc)
add_executable("mk2" test.cc)
add_executable("oned" oned.cc)
target_link_dune_default_libraries("dune-integrate")
