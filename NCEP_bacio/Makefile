#!/bin/make
##########################################################################
#!BOP
#
# !ROUTINE:        Makefile for 'baciolib' 
# 
#
# !DESCRIPTION: 
#       create 'baciolib' utility routines
#
# !REVISION HISTORY: 
#
# 15Jul2004 Todling            Added missing routines found in paqc/w3lib
#                               (attempting to have a single w3 lib)
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
CFLAGS		= -O2
CPP		= /lib/cpp
CPPFLAGS	= -P
FPP	        = $(CPP)

FC              = f90
FFLAGS          = -O2
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

LIB  =	libbacio.a
LIBS =	libbacio_r8i8.a libbacio_r4i4.a

SRCS  = $(OBJ_MOD:.o=.f) $(OBJS:.o=.f) $(OBJSCC:.o=.c)

OBJ_MOD	= bacio_module.o

OBJS  = baciof.o 

OBJSCC = bacio.v1.3.o

#...................................................................

lib:    $(LIB)

export:
	$(MKDIR) $(CORELIB) $(COREINC)/bacio4
	-$(RM) clib.h
	$(LN) clib4.h clib.h
	$(MAKE) _RPREC=4 _IPREC=4 LIB=libbacio_r4i4.a
	$(CP) clib.h   $(COREINC)/bacio4
	$(CP) *.mod    $(COREINC)/bacio4
	-$(RM) clib.h *.mod
	$(CP) libbacio_r4i4.a $(CORELIB)
	$(MAKE) clean
	$(MKDIR) $(CORELIB) $(COREINC)/bacio8
	$(LN) clib8.h clib.h
	$(MAKE) _RPREC=8 _IPREC=8 LIB=libbacio_r8i8.a
	$(CP) clib.h   $(COREINC)/bacio8
	$(CP) *.mod    $(COREINC)/bacio8
	-$(RM) clib.h *.mod
	$(CP) libbacio_r8i8.a $(CORELIB)
	$(MAKE) clean

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
include Makefile.depend

