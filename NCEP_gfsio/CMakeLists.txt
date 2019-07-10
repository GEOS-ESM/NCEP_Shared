esma_set_this ()

esma_add_library (${this}
  SRCS gfsio_module.f90 gfsio_solve_axb.F90
  DEPENDENCIES GMAO_mpeu NCEP_bacio_r4i4 NCEP_w3_r4i4 ${MKL_LIBRARIES}
  )

set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${BIG_ENDIAN} ${common_Fortran_flags} ${GEOS_Fortran_Release_FPE_Flags} ${ALIGNCOM}")
