#!/bin/make
##########################################################################
#!BOP
#
# !ROUTINE:        Makefile for 'splib' 
# 
#
# !DESCRIPTION: 
#       create 'splib' utility routines
#
# !REVISION HISTORY: 
#
# 97.12.16  Meta Sienkiewicz   Initial Makefile for bufron29
# 26Apr2001 da Silva           Made fvDAS compatible
# 14Jul2004 Todling            Expanded export to create various prec. libs;
#                              (the only one guarantee to work is r8i4)
# 
#!EOP
##########################################################################

COREROOT = ../../..
COREBIN  = $(COREROOT)/bin
CORELIB  = $(COREROOT)/lib
COREINC  = $(COREROOT)/include
COREETC  = $(COREROOT)/etc

SHELL		= sh
MAKE		= make
MAKEFILE	= Makefile
MAKECONF        = Makefile.conf

CPP		= /lib/cpp
CPPFLAGS	= -P
FPP	        = $(CPP)

FC              = f90
FFLAGS          =  
F90		= $(FC)
F90FLAGS	= $(FFLAGS)
F90PP		= $(FPP)
F90PPFLAGS	= $(FPPFLAGS)

LN		= ln -s
RM		= rm -f
CP		= cp -p
MKDIR  = /bin/mkdir -p

ZIP		= compress
ZIPTAG		= .Z
ZCAT		= zcat
TAR		= tar
AR		= ar cq

LIBESSL 	=
LIBSESSL 	=

#..................................................................
.SUFFIXES:
.SUFFIXES: .H .h .F90 .f90 .F .f .o

.f.o:
	$(FC) -c $(FFLAGS) $*.f

.F.o:
	@$(RM) $*.f
	$(FPP) $(FPPFLAGS) $(DEFS) $*.F $*.f
	$(FC) -c $(FFLAGS) $*.f
	@$(RM) $*.f

.f90.o:
	@if [ -f $*.[Ff] ]; then \
	   echo "make: WARNING! coexistent $*.[Ff] and $*.f90" ; fi
	$(F90) -c $(F90FLAGS) $*.f90

.F90.o:
	@ if [ -f $*.[Ff] ]; then \
	   echo "make: WARNING! coexistent $*.[Ff] and $*.F90" ; fi
	@ $(RM) $*.f90
	$(F90PP) $(F90PPFLAGS) $(DEFS) $*.F90 $*.f90
	$(F90) -c $(F90FLAGS) $*.f90
	@ $(RM) $*.f90


# Architecture/site specific costumization
# ----------------------------------------
include $(MAKECONF)


#...................................................................

LIB	 = libsp.a
LIBS     = libsp.a    libsp_r8i8.a    libsp_r4i4.a

SRCS  = $(OBSOBJ:.o=.f) 

OBJS  =               bll2ps.o ncpus.o spanaly.o spdz2uv.o speps.o spffte.o \
        spgradq.o splaplac.o splat.o splegend.o sppad.o spsynth.o sptez.o \
        sptezd.o sptezm.o sptezmv.o sptezv.o sptgpm.o sptgpmd.o sptgpmv.o \
	sptgps.o sptgpsd.o sptgpsv.o sptgpt.o sptgptd.o sptgptv.o sptran.o \
	sptrand.o sptranf.o sptranf0.o sptranf1.o sptranfv.o sptranv.o sptrun.o \
	sptrund.o sptrung.o sptrungv.o sptrunl.o sptrunm.o sptrunmv.o sptruns.o \
	sptrunsv.o sptrunv.o spuv2dz.o spvar.o spwget.o

#...................................................................

lib:    $(LIB)

export: $(LIBSESSL)
	$(MKDIR) $(CORELIB) 
	make clean
	make LIB=libsp.a
	$(CP) libsp.a                                 $(CORELIB)
	make clean
	make _RPREC=8 _IPREC=8 LIB=libsp_r8i8.a
	$(CP) libsp_r8i8.a                            $(CORELIB)
	make clean
	make _RPREC=4 _IPREC=4 LIB=libsp_r4i4.a
	$(CP) libsp_r4i4.a                            $(CORELIB)
	make clean

$(LIB): $(OBJS)
	$(AR) $@ $(OBJS)

#
# Added to present missing routines to transf when
# building in the NCEP IBM ... essl.a under /usr/lib is incomplete
# 
libessl2_r8i8.a,libessl2_r4i4.a: essl.o
	$(AR) $@ essl.o
libessl2.a: essl.o
	$(AR) $@ essl.o

cleanobjs:
	$(RM) -r *~ *.[ol] *.mod

clean:
	$(RM) -r *~ *.[axol] *.mod *.[fF]*~ Makefile*~ rii_files/ \
	         *.aux *.dvi *.ps  *.toc *.log 

distclean:
	make clean
	$(RM) Makefile.conf









