file(REMOVE_RECURSE
  "libnbots.pdb"
  "libnbots.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/nbots_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
