file(REMOVE_RECURSE
  "lib..a"
  "lib..pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/..dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
