esma_set_this ()

esma_add_library (${this}
  SRCS gfsio_modul.f90 gfsio_solve_axb.F90
  DEPENDENCIES GMAO_mpeu NCEP_bacio_r4i4 NCEP_w3_r4i4 ${MKL_LIBRARIES}
  )
