esma_set_this ()

esma_add_library (${this}
  SRCS nemsio_module.f90  nemsio_openclose.f90  nemsio_read.f90  nemsio_write.f90
  DEPENDENCIES NCEP_bacio_r4i4 NCEP_w3_r4i4 GMAO_mpeu
  )
