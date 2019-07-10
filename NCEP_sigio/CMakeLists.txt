esma_set_this()

include_directories (${include_${this}})

set (srcs
   sigio_module.f90 sigio_r_module.f90
  )

esma_add_library(${this} SRCS ${srcs})
set_target_properties (${this} PROPERTIES Fortran_MODULE_DIRECTORY ${include_${this}})

set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${BIG_ENDIAN} ${common_Fortran_flags} ${GEOS_Fortran_Release_FPE_Flags} ${ALIGNCOM}")

