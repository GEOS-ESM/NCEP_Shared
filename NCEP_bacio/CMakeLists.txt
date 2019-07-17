esma_set_this ()
string (REPLACE NCEP_bacio_ "" precision ${this})

  
set (srcs
   baciof.f bafrio.f chk_endianc.f bacio.c byteswap.c
   )


set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${GEOS_Fortran_Release_FPE_Flags} ${ALIGNCOM}")

esma_add_library(${this} SRCS ${srcs})

set (CMAKE_C_FLAGS_RELEASE "-O -I . -DLINUX -Dfunder -DFortranByte=char -DFortranInt=int -DFortranLlong='long long' -O3 -fPIC")

# MAT Note that all the funder FortranByte, etc. in the GNUmake file are not used!
set (LOCAL_CDEFS LINUX)
foreach (flag ${LOCAL_CDEFS})
   target_compile_definitions (${this} PRIVATE $<$<COMPILE_LANGUAGE:C>:${flag}>)
endforeach ()

if (precision MATCHES "r4i4")
  # use default real/int precisions

  string(REPLACE " " ";" rtmp "_RPREC=4")
  foreach (flag ${rtmp})
    target_compile_definitions (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
    target_compile_definitions (${this} PRIVATE $<$<COMPILE_LANGUAGE:C>:${flag}>)
  endforeach ()

  string(REPLACE " " ";" itmp "_IPREC=4")
  foreach (flag ${itmp})
    target_compile_definitions (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
    target_compile_definitions (${this} PRIVATE $<$<COMPILE_LANGUAGE:C>:${flag}>)
  endforeach ()
  target_compile_definitions (${this} PRIVATE $<$<COMPILE_LANGUAGE:C>:CLIB4>)

elseif (precision MATCHES "r8i4")
  string(REPLACE " " ";" tmp ${FREAL8})
  foreach (flag ${tmp})
    target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
  endforeach ()

  string(REPLACE " " ";" rtmp "_RPREC=8")
  foreach (flag ${rtmp})
    target_compile_definitions (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
    target_compile_definitions (${this} PRIVATE $<$<COMPILE_LANGUAGE:C>:${flag}>)
  endforeach ()

  string(REPLACE " " ";" itmp "_IPREC=4")
  foreach (flag ${itmp})
    target_compile_definitions (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
    target_compile_definitions (${this} PRIVATE $<$<COMPILE_LANGUAGE:C>:${flag}>)
  endforeach ()
  target_compile_definitions (${this} PRIVATE $<$<COMPILE_LANGUAGE:C>:CLIB4>)

elseif (precision MATCHES "r8i8")
  string(REPLACE " " ";" tmp ${FREAL8})
  foreach (flag ${tmp})
    target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
  endforeach ()

  string(REPLACE " " ";" rtmp "_RPREC=8")
  foreach (flag ${rtmp})
    target_compile_definitions (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
    target_compile_definitions (${this} PRIVATE $<$<COMPILE_LANGUAGE:C>:${flag}>)
  endforeach ()

  string(REPLACE " " ";" tmp ${FINT8})
  foreach (flag ${tmp})
    target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
  endforeach ()

  string(REPLACE " " ";" itmp "_IPREC=8")
  foreach (flag ${itmp})
    target_compile_definitions (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
    target_compile_definitions (${this} PRIVATE $<$<COMPILE_LANGUAGE:C>:${flag}>)
  endforeach ()
  target_compile_definitions (${this} PRIVATE $<$<COMPILE_LANGUAGE:C>:CLIB8>)

endif ()
    
