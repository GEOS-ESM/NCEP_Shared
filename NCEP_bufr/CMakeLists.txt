esma_set_this ()
string (REPLACE NCEP_bufr_ "" kinds ${this})

string (SUBSTRING ${kinds} 0 2 REAL_KIND)
string (SUBSTRING ${kinds} 2 2 INTEGER_KIND)

set (srcs
   adn30.f
   arallocc.c
   arallocf.F
   ardllocc.c
   ardllocf.F
   atrcpt.f
   bfrini.f
   blocks.f
   bort.f
   bort2.f
   bort_exit.c
   bufrlib.h
   bufrlib.prm
   bvers.f
   cadn30.f
   capit.f
   ccbfl.c
   cfe.h
   chekstab.f
   chrtrn.f
   chrtrna.f
   cktaba.f
   closbf.f
   closmg.f
   cmpia.c
   cmpmsg.f
   cmpstia1.c
   cmpstia2.c
   cmsgini.f
   cnved4.f
   cobfl.c
   cobfl.h
   codflg.f
   conwin.f
   copybf.f
   copymg.f
   copysb.f
   cpbfdx.f
   cpdxmm.f
   cpmstabs.c
   cpymem.f
   cpyupd.f
   crbmg.c
   cread.c
   cread.h
   cwbmg.c
   datebf.f
   datelen.f
   digit.f
   dlloctbf.c
   drfini.f
   drstpl.f
   dumpbf.f
   dxdump.f
   dxinit.f
   dxmini.f
   elemdx.f
   errwrt.f
   exitbufr.F
   fstag.f
   getabdb.f
   getbmiss.f
   getcfmng.f
   getlens.f
   getntbe.f
   gets1loc.f
   gettagpr.f
   gettagre.f
   gettbh.f
   getvalnb.f
   getwin.f
   hold4wlc.f
   i4dy.f
   ibfms.f
   icbfms.f
   ichkstr.f
   icmpdx.f
   icopysb.f
   icvidx.c
   idn30.f
   idxmsg.f
   ifbget.f
   ifxy.f
   igetdate.f
   igetfxy.f
   igetmxby.f
   igetntbi.f
   igetntbl.f
   igetprm.f
   igetrfel.f
   igetsc.f
   igettdi.f
   imrkopr.f
   inctab.f
   inittbf.c
   invcon.f
   invmrg.f
   invtag.f
   invwin.f
   iok2cpy.f
   iokoper.f
   ipkm.f
   ipks.f
   ireadmg.f
   ireadmm.f
   ireadmt.F
   ireadns.f
   ireadsb.f
   irev.F
   isetprm.F
   ishrdx.f
   isize.f
   istdesc.f
   iupb.f
   iupbs01.f
   iupbs3.f
   iupm.f
   iupvs01.f
   jstchr.f
   jstnum.f
   lcmgdf.f
   lmsg.f
   lstjpb.f
   makestab.f
   maxout.f
   mesgbc.f
   mesgbf.f
   minimg.f
   moda_bitbuf.F
   moda_bitmaps.F
   moda_bufrmg.F
   moda_bufrsr.F
   moda_comprs.F
   moda_comprx.F
   moda_h4wlc.F
   moda_idrdm.F
   moda_ifopbf.F
   moda_ival.F
   moda_ivttmp.F
   moda_lushr.F
   moda_mgwa.F
   moda_mgwb.F
   moda_msgcwd.F
   moda_msglim.F
   moda_msgmem.F
   moda_mstabs.F
   moda_nmikrp.F
   moda_nrv203.F
   moda_nulbfr.F
   moda_rdmtb.F
   moda_rlccmn.F
   moda_s01cm.F
   moda_sc3bfr.F
   moda_stbfr.F
   moda_stcode.F
   moda_tababd.F
   moda_tables.F
   moda_ufbcpl.F
   moda_unptyp.F
   moda_usrbit.F
   moda_usrint.F
   moda_usrtmp.F
   moda_xtab.F
   modv_MAXCD.F
   modv_MAXJL.F
   modv_MAXMEM.F
   modv_MAXMSG.F
   modv_MAXSS.F
   modv_MAXTBA.F
   modv_MAXTBB.F
   modv_MAXTBD.F
   modv_MXBTM.F
   modv_MXBTMSE.F
   modv_MXCDV.F
   modv_MXCSB.F
   modv_MXDXTS.F
   modv_MXH4WLC.F
   modv_MXLCC.F
   modv_MXMSGL.F
   modv_MXMTBB.F
   modv_MXMTBD.F
   modv_MXMTBF.F
   modv_MXNRV.F
   modv_MXRST.F
   modv_MXS01V.F
   modv_MXTAMC.F
   modv_MXTCO.F
   modv_NFILES.F
   mrginv.f
   msgfull.f
   msgini.f
   msgupd.f
   msgwrt.f
   mstabs.h
   mtfnam.f
   mtinfo.f
   mvb.f
   nemdefs.f
   nemock.f
   nemspecs.f
   nemtab.f
   nemtba.f
   nemtbax.f
   nemtbb.f
   nemtbd.f
   nenubd.f
   nevn.f
   newwin.f
   nmsub.f
   nmwrd.f
   numbck.f
   nummtb.c
   numtab.f
   numtbd.f
   nvnwin.f
   nwords.f
   nxtwin.f
   openbf.F
   openbt.f
   openmb.f
   openmg.f
   pad.f
   padmsg.f
   parstr.f
   parusr.f
   parutg.f
   pkb.f
   pkbs1.f
   pkc.f
   pkftbv.f
   pktdd.f
   pkvs01.F
   pkx.f
   posapx.f
   rbytes.c
   rcstpl.f
   rdbfdx.f
   rdcmps.f
   rdmemm.f
   rdmems.f
   rdmgsb.f
   rdmsgw.f
   rdmtbb.f
   rdmtbd.f
   rdmtbf.f
   rdtree.f
   rdusdx.f
   readdx.f
   readerme.f
   readlc.f
   readmg.f
   readmm.f
   readns.f
   reads3.f
   readsb.f
   restd.c
   rewnbf.f
   rjust.f
   rsvfvm.f
   rtrcpt.f
   rtrcptb.f
   seqsdx.f
   setblock.f
   setbmiss.f
   setvalnb.f
   sntbbe.f
   sntbde.f
   sntbfe.f
   sorttbf.c
   srchtbf.c
   status.f
   stbfdx.f
   stdmsg.f
   stndrd.f
   stntbi.f
   stntbia.f
   strbtm.f
   strcln.f
   strcpt.f
   string.f
   strnum.f
   strsuc.f
   strtbfe.c
   stseq.c
   tabent.f
   tabsub.f
   trybump.f
   ufbcnt.f
   ufbcpy.f
   ufbcup.f
   ufbdmp.f
   ufbevn.f
   ufbget.f
   ufbin3.f
   ufbint.f
   ufbinx.f
   ufbmem.f
   ufbmex.f
   ufbmms.f
   ufbmns.f
   ufbovr.f
   ufbpos.f
   ufbqcd.f
   ufbqcp.f
   ufbrep.f
   ufbrms.f
   ufbrp.f
   ufbrw.f
   ufbseq.f
   ufbsp.f
   ufbstp.f
   ufbtab.f
   ufbtam.f
   ufdump.f
   upb.f
   upbb.f
   upc.f
   upds3.f
   upftbv.f
   ups.f
   uptdd.f
   usrtpl.f
   valx.f
   wrcmps.f
   wrdesc.c
   wrdlen.F
   wrdxtb.f
   writcp.f
   writdx.f
   writlc.f
   writsa.f
   writsb.f
   wrtree.f
   wtstat.f
  )

add_custom_command (
  OUTPUT bufrlib.prm
  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/bufrlib.PRM bufrlib.prm
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
esma_add_library (${this} SRCS ${srcs})

set (CMAKE_C_FLAGS_RELEASE "-O")

# This is equivalent to FOPT=$(FOPT2) in GNU Make
if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   string (REPLACE "${FOPT3}" "${FOPT2}" CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS_RELEASE})
   string (REPLACE "${OPTREPORT0}" "" CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS_RELEASE})
   string (REPLACE "${FTZ}" "" CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS_RELEASE})
   string (REPLACE "${ALIGN_ALL}" "" CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS_RELEASE})
   string (REPLACE "${NO_ALIAS}" "" CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS_RELEASE})
endif ()

if (REAL_KIND STREQUAL "r8")
  string (REPLACE " " ";" tmp ${FREAL8})
  foreach (flag ${tmp})
    target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
  endforeach()
endif()

if (INTEGER_KIND STREQUAL "i8")
  string (REPLACE " " ";" tmp ${FINT8})
  foreach (flag ${tmp})
    target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
  endforeach()
  target_compile_definitions (${this} PRIVATE -DF77_INTSIZE_8)
endif()


target_compile_definitions (${this} PRIVATE $<$<COMPILE_LANGUAGE:C>:STATIC_ALLOCATION>)
target_compile_definitions (${this} PRIVATE $<$<COMPILE_LANGUAGE:C>:MAXNC=600>)
target_compile_definitions (${this} PRIVATE $<$<COMPILE_LANGUAGE:C>:MXNAF=3>)
target_compile_definitions (${this} PRIVATE $<$<COMPILE_LANGUAGE:C>:UNDERSCORE>)

