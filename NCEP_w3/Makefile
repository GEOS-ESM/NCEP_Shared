#!/bin/make
##########################################################################
#!BOP
#
# !ROUTINE:        Makefile for 'w3lib' 
# 
#
# !DESCRIPTION: 
#       create 'w3lib' utility routines
#
# !REVISION HISTORY: 
#
# 97.12.16  Meta Sienkiewicz   Initial Makefile for bufron29
# 26Apr2001 da Silva           Made fvDAS compatible
# 14Jul2004 Todling            Expanded export to create various prec. libs
# 15Jul2004 Todling            Added missing routines found in paqc/w3lib
#                               (attempting to have a single w3 lib)
# 17Feb2005 Todling            Moved all bacio stuff to bacio lib
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
MAKECONF	= Makefile.conf

CC		= cc
CFLAGS		= # -O2 -DIBM4
CPP		= /lib/cpp
CPPFLAGS	= -P
FPP	        = $(CPP)

FC              = f90
FFLAGS          = # -O2
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

#..................................................................
# Architecture/site specific costumization
# ----------------------------------------
include $(MAKECONF)


.SUFFIXES: .H .h .F90 .f90 .F .f .c .o

.c.o:
	$(CC) $(CFLAGS) -c $<

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

#...................................................................

LIB  =	libw3.a
LIBS =	libw3.a libw3_r8i8.a libw3_r4i4.a

SRCS  = $(OBJ_MOD:.o=.f) $(OBJS:.o=.f) $(OBJSCC:.o=.c)

OBJ_MOD	= 

OBJS  = gbytes_char.o getgb.o     getgbmh.o   putgbex.o   w3fi73.o \
	getgb1r.o   getgbmp.o   putgbn.o    w3fi74.o \
	getgb1re.o  getgbp.o    r63w72.o    w3fi75.o \
	getgb1s.o   getgi.o     sbyte.o     w3fi76.o \
	getgbe.o    getgir.o    sbytes.o    w3fi82.o \
	getgbeh.o   idsdef.o    skgb.o      w3fi83.o \
	getgbem.o   iw3jdn.o    w3difdat.o  w3fs26.o \
	getgbemh.o  ixgb.o      w3doxdat.o  w3log.o \
	getgbemp.o  lengds.o    w3fi01.o    w3movdat.o \
	getgbens.o  pdsens.o    w3fi58.o    w3reddat.o \
	getgbep.o   pdseup.o    w3fi59.o    w3tagb.o \
	errmsg.o    getgbex.o   putgb.o     w3fi63.o    \
	gbytes.o    getgbexm.o  putgbe.o    w3fi68.o    \
	mova2i.o    getgbh.o    putgben.o   w3fi71.o \
	getbit.o    getgbm.o    putgbens.o  w3fi72.o \
	errexit.o   fparsei.o   fparser.o   instrument.o  \
	start.o     summary.o   w3utcdat.o  w3fs21.o \
        w3locdat.o                          orders.o \
        w3fc05.o    w3fs13.o    w3fi04.o    w3unpkb7.o \
        iw3unpbf.o  w3fa03.o    w3fa04.o    w3rtveds.o \
        w3tovmnd.o  w3cnvtov.o  w3rtvunp.o  w3fb04.o \
        w3ersunb.o  w3pradat.o

#       w3fa03.o   w3fa04.o   w3fa07.o   w3fb04.o \
#       w3fc05.o   w3fi02.o   w3fi03.o   w3fi04.o \
#       w3fi43.o   w3fi64.o   w3fi65.o   w3fi66.o \
#       w3fi77.o   w3fi88.o   w3fs03.o   w3fs06.o \
#       w3fs11.o   w3fs13.o   w3fs15.o   w3fs22.o \
#       w3rtveds.o w3rtvunp.o w3toveds.o w3tovmnd.o \
#       w3tovunp.o w3unpk77.o w3unpkb7.o xdopen.o xstore.o

OBJSCC = 

#...................................................................

lib:    $(LIB)

export:
	$(MKDIR) $(CORELIB) 
	make clean
	make LIB=libw3.a
	$(CP) libw3.a                                $(CORELIB)
	make clean
	make _RPREC=8 _IPREC=8 LIB=libw3_r8i8.a
	$(CP) libw3_r8i8.a                           $(CORELIB)
	make clean
	make _RPREC=4 _IPREC=4 LIB=libw3_r4i4.a
	$(CP) libw3_r4i4.a                           $(CORELIB)
	make clean

$(LIB): $(OBJ_MOD) $(OBJS) $(OBJSCC)
	$(AR) $@ $(OBJ_MOD) $(OBJS) $(OBJSCC)

cleanobjs:
	$(RM) -r *~ *.[ol] *.mod

clean:
	$(RM) -r *~ *.[axol] *.mod *.[fF]*~ Makefile*~ rii_files/ \
	         *.aux *.dvi *.ps  *.toc *.log 

distclean:
	make clean
	$(RM) Makefile.conf


#.........................................................................


