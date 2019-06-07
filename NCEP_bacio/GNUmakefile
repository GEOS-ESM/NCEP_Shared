#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 20oct2004  da Silva  First ESMA version.
# 12Apr2005  da Silva  Revised.
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../../..
endif

THIS := $(shell basename `pwd`)

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

esma_help help:
	@echo "Standard ESMA targets:"
	@echo "% make esma_install    (builds and install under ESMADIR)"
	@echo "% make esma_clean      (removes deliverables: *.[aox], etc)"
	@echo "% make esma_distclean  (leaves in the same state as cvs co)"
	@echo "% make esma_doc        (generates PDF, installs under ESMADIR)"
	@echo "% make esma_help       (this message)"
	@echo "Environment:"
	@echo "      ESMADIR = $(ESMADIR)"
	@echo "      BASEDIR = $(BASEDIR)"
	@echo "         ARCH = $(ARCH)"
	@echo "         SITE = $(SITE) "
	@echo "         THIS = $(THIS) "


#                  --------------------------------
#                   Recurse Make in Sub-directories
#                  --------------------------------

export ESMADIR

SUBDIRS = r4i4 r8i4 r8i8

TARGETS = esma_install esma_clean esma_doc \
          install clean doc

$(TARGETS): 
	@ t=$@; argv="$(SUBDIRS)" ;\
	  for d in $$argv; do			 \
	    ( $(MKDIR) $$d			;\
              cd $$d			        ;\
	      echo ""; echo Making $$t for $$d $(THIS) in $(PWD) ;\
	      $(MAKE) -f ../GNUmakefile_rnin -e $$t NTHIS=$(THIS) DIRN=$$d ) \
	  done
	$(MAKE) local_$@

local_esma_install local_install: $(LIB)
	@echo No local install in here

local_esma_clean local_clean:
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd]

esma_distclean distclean:
	-$(RM) -r *~ *.[aoxd] *.[Mm][Oo][Dd] $(SUBDIRS)

#                  --------------------
#                  User Defined Targets
#                  --------------------

r4i4_install:
	$(MKDIR) r4i4
	echo Making install for r4i4 in $(PWD)
	$(MAKE) -C r4i4 -f ../GNUmakefile_rnin -e install NTHIS=$(THIS) DIRN=r4i4

r8i4_install:
	$(MKDIR) r8i4
	echo Making install for r8i4 in $(PWD)
	$(MAKE) -C r8i4 -f ../GNUmakefile_rnin -e install NTHIS=$(THIS) DIRN=r8i4

r8i8_install:
	$(MKDIR) r8i8
	echo Making install for r8i8 in $(PWD)
	$(MAKE) -C r8i8 -f ../GNUmakefile_rnin -e install NTHIS=$(THIS) DIRN=r8i8

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
