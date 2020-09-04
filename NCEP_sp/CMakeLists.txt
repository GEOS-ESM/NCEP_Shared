esma_set_this ()
string (REPLACE NCEP_sp_ "" precision ${this})

# TODO: there should be a platform switch for spffte.f (e.g., SGI) TLC

set (srcs
  bll2ps.f ncpus.f spanaly.f spdz2uv.f speps.f spffte.f
  spgradq.f splaplac.f splat.f splegend.f sppad.f spsynth.f
  sptez.f sptezd.f sptezm.f sptezmv.f sptezv.f sptgpm.f sptgpmd.f
  sptgpmv.f sptgps.f sptgpsd.f sptgpsv.f sptgpt.f sptgptd.f sptgptv.f
  sptran.f  sptrand.f sptranf.f sptranf0.f sptranf1.f
  sptranfv.f sptranv.f sptrun.f
  sptrund.f sptrung.f sptrungv.f sptrunl.f sptrunm.f sptrunmv.f
  sptruns.f sptrunsv.f sptrunv.f spuv2dz.f spvar.f spwget.f
  )

set (CMAKE_Fortran_FLAGS_RELEASE "-O2 ${GEOS_Fortran_Release_FPE_Flags} ${EXTENDED_SOURCE} ${ALIGNCOM}")

esma_add_library (${this} SRCS ${srcs} DEPENDENCIES NCEP_w3_${precision})

if (APPLE)
   set(MKL_Fortran TRUE)
endif ()
find_package(MKL)
if (MKL_FOUND)
   set(BLA_VENDOR Intel10_64lp_seq)
   target_link_libraries(${this} PRIVATE ${MKL_LIBRARIES})
else ()
   find_package(LAPACK REQUIRED)
   target_link_libraries(${this} PRIVATE ${LAPACK_LIBRARIES})
endif ()

if (precision MATCHES "r4i4")
  # use default real/int precisions
elseif (precision MATCHES "r8i4")
  string(REPLACE " " ";" tmp ${FREAL8})
  foreach (flag ${tmp})
    target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
  endforeach ()
elseif (precision MATCHES "r8i8")
  string(REPLACE " " ";" tmp ${FREAL8})
  foreach (flag ${tmp})
    target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
  endforeach ()
  target_compile_options (${this} PRIVATE ${FINT8})
endif ()
    
