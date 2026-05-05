esma_set_this()

include_directories (${include_${this}})

set (srcs
   sigio_module.f90 sigio_r_module.f90
  )

esma_add_library(${this} SRCS ${srcs})
set_target_properties (${this} PROPERTIES Fortran_MODULE_DIRECTORY ${include_${this}})

# NAG requires mismatch flag for intentional type-punning with C byte I/O routines
if (CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  target_compile_options(${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-wmismatch=BAFRREADL,BAFRWRITEL,BYTESWAP>)
endif()

set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${BIG_ENDIAN} ${common_Fortran_flags} ${GEOS_Fortran_Release_FPE_Flags} ${ALIGNCOM}")

