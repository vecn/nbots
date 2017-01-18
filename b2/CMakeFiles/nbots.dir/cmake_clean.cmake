file(REMOVE_RECURSE
  "libnbots.pdb"
  "libnbots.so"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/nbots.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
