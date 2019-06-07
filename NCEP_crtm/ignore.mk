ifdef __ignored__
__IGNORED_DIR__	?= __ignored__
VPATH += $(__IGNORED_DIR__)
-include $(__IGNORED_DIR__)/__.mk

$(__IGNORED_DIR__)/__.mk: $(__IGNORED_DIR__)
$(__ignored__): $(__IGNORED_DIR__)

$(__IGNORED_DIR__):
	@ echo "$(MAKE) >>> Making ($@) ... <<<"
	mkdir -p $@; touch $@/__.mk
	for t in $(__ignored__); do touch $@/$$t; done

local_distclean distclean: distclean-ignored
distclean-ignored:
	@ echo "$(MAKE) >>> Making ($@) ... <<<"
	rm -fr $(__IGNORED_DIR__)
endif
