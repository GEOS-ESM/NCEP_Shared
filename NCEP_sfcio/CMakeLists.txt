esma_set_this()

include_directories (${include_NCEP_bacio_r4i4})
include_directories (${include_${this}})

set (srcs
   sfcio_module.f90 testsfcio.f90
  )

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES NCEP_bacio_r4i4)
set_target_properties (${this} PROPERTIES Fortran_MODULE_DIRECTORY ${include_${this}})

set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${common_Fortran_flags} ${GEOS_Fortran_Release_FPE_Flags} ${ALIGNCOM}")

set_source_files_properties (sfcio_module.f90 PROPERTIES COMPILE_FLAGS ${FREAL8})
set_source_files_properties (testsfcio.f90 PROPERTIES COMPILE_FLAGS ${FREAL8})

