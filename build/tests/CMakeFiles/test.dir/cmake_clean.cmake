FILE(REMOVE_RECURSE
  "CMakeFiles/test.dir/test.c.o"
  "libtest.pdb"
  "libtest.a"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang C)
  INCLUDE(CMakeFiles/test.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
