# Default configuration, possibily redefined in Makefile.conf
# -----------------------------------------------------------
                                                                                                                                               
COREROOT = ../../..
COREBIN  = $(COREROOT)/bin
CORELIB  = $(COREROOT)/lib
COREINC  = $(COREROOT)/include
COREETC  = $(COREROOT)/etc
                                                                                                                                               
SHELL           = /bin/sh
                                                                                                                                               
ARCH            = `uname -s`
_Ddft           =
_Dsys           = -Dsys$(ARCH)

# -------------
# General Rules
# -------------
                                                                                                                                               
CP              = /bin/cp -p
RM              = /bin/rm -f
MKDIR           = /bin/mkdir -p
AR              = ar
ARFLAGS         = cq
PROTEX          = protex -f # -l
ProTexMake      = protex -S # -l
LATEX           = pdflatex
DVIPS           = dvips
                                                                                                                                               
MAKE            = make
MAKEFILE        = Makefile
MAKECONF        = Makefile.conf

LIB 		= libsfcio.a

include $(MAKECONF)

SRCSf90	= sfcio_module.f90 testsfcio.f90

OBJS	= ${SRCSf90:.f90=.o}

# -----------------------
# Default compiling rules
# -----------------------
                                                                                                                                               
.SUFFIXES : .F90 .f90 .f .c
                                                                                                                                               
.F90.o  :
	@echo
	@echo '---> Compiling $<'
	$(FC) $(FFLAGS) -c $<
                                                                                                                                               
.f90.o  :
	@echo
	@echo '---> Compiling $<'
	$(FC) $(FFLAGS) -c $<
                                                                                                                                               
.f.o  :
	@echo
	@echo '---> Compiling $<'
	$(FC) $(FFLAGS_f) -c $<
                                                                                                                                               
.c.o  :
	@echo
	@echo '---> Compiling $<'
	$(CC) $(CPLAGS) -c $<


lib: $(LIB)
                                                                                                                                               
$(LIB): $(OBJS)
	$(RM) $(LIB)
	$(AR) $(ARFLAGS) $@ $(OBJS)

testsfcio.x: $(OBJS) $(LIB) testsfcio.o
	$(FC) $(LDFLAGS) -o testsfcio.x testsfcio.o $(LIB)

export: $(LIB)
	@if [ ! -d $(CORELIB) ]; then \
	   $(MKDIR) $(CORELIB) ;\
	fi
	        $(CP) $(LIB) $(CORELIB)/
	@if [ ! -d $(COREINC)/sfcio ]; then \
	   $(MKDIR) $(COREINC)/sfcio ;\
	fi
	$(CP) *.mod  $(COREINC)/sfcio

# ----------
# Make clean
# ----------
                                                                                                                                        
clean:
	@echo
	@echo '==== CLEAN ==================================================='
	- $(RM) $(EXE_FILE) *.o *.mod *.MOD *.lst *.a *.x
	- $(RM) loadmap.txt log.make.global_anl
                                                                                                                                        

