esma_check_if_debug()

esma_add_subdirectories(
  NCEP_crtm
  NCEP_sfcio
  NCEP_sigio
  NCEP_irsse
  NCEP_gfsio
  NCEP_nemsio
  )

if (EXISTS NCEP_sp)
  add_subdirectory (NCEP_sp NCEP_sp_r4i4)
  add_subdirectory (NCEP_sp NCEP_sp_r8i4)
  add_subdirectory (NCEP_sp NCEP_sp_r8i8)
endif ()

if (EXISTS NCEP_bacio)
  add_subdirectory (NCEP_bacio NCEP_bacio_r4i4)
  add_subdirectory (NCEP_bacio NCEP_bacio_r8i4)
  add_subdirectory (NCEP_bacio NCEP_bacio_r8i8)
endif ()
  
if (EXISTS NCEP_w3)
  add_subdirectory (NCEP_w3 NCEP_w3_r4i4)
  add_subdirectory (NCEP_w3 NCEP_w3_r8i4)
  add_subdirectory (NCEP_w3 NCEP_w3_r8i8)
endif ()
  
if (EXISTS NCEP_bufr)
  add_subdirectory (NCEP_bufr NCEP_bufr_r4i4)
  add_subdirectory (NCEP_bufr NCEP_bufr_r8i4)
  add_subdirectory (NCEP_bufr NCEP_bufr_r8i8)
endif ()
