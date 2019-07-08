esma_set_this ()
string (REPLACE NCEP_w3_ "" precision ${this})

  
# TODO: there should be a platform switch for spffte.f (e.g., SGI) TLC

set (srcs
   errexit.f getgbens.f iw3unpbf.f sbytes.f w3fi59.f w3log.f errmsg.f
   getgbep.f ixgb.f skgb.f w3fi63.f w3movdat.f fparsei.f getgbex.f
   lengds.f start.f w3fi68.f w3pradat.f fparser.f getgbexm.f mova2i.f
   summary.f w3fi71.f w3reddat.f gbytes_char.f getgb.f orders.f
   w3cnvtov.f w3fi72.f w3rtveds.f gbytes.f getgbh.f pdsens.f w3difdat.f
   w3fi73.f w3rtvunp.f getbit.f getgbm.f pdseup.f w3doxdat.f w3fi74.f
   w3tagb.f getgb1re.f getgbmh.f putgbe.f w3ersunb.f w3fi75.f w3tovmnd.f
   getgb1r.f getgbmp.f putgben.f w3fa03.f w3fi76.f w3unpkb7.f getgb1s.f
   getgbp.f putgbens.f w3fa04.f w3fi82.f w3utcdat.f getgbe.f getgi.f
   putgbex.f w3fb04.f w3fi83.f getgbeh.f getgir.f putgb.f w3fc05.f
   w3fs13.f getgbem.f idsdef.f putgbn.f w3fi01.f w3fs21.f getgbemh.f
   instrument.f r63w72.f w3fi04.f w3fs26.f getgbemp.f iw3jdn.f sbyte.f
   w3fi58.f w3locdat.f 
   w3fb06.f w3fb08.f w3fb11.f w3fb12.f w3valdat.f w3fb07.f w3fb09.f 
   mersenne_twister.f
   )

set (CMAKE_Fortran_FLAGS_RELEASE "-O2 ${common_Fortran_flags} ${GEOS_Fortran_Release_FPE_Flags} ${FIXED_SOURCE} ${ALIGNCOM}")

esma_add_library (${this} SRCS ${srcs})
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

ecbuild_add_executable(
   TARGET ndate_${precision}.x
   SOURCES ndate.f
   LIBS NCEP_w3_${precision}
   )


