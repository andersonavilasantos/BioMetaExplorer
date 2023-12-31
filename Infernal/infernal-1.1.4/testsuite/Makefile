top_srcdir     = ..
srcdir         = .


SQCLEVEL    = 1
DEVSQCLEVEL = 2

ESLDIR  = easel

# beautification magic stolen from git 
#
QUIET_SUBDIR0 = +${MAKE} -C #space separator after -c
QUIET_SUBDIR1 = 
ifndef V
	QUIET_CC      = @echo '    ' CC $@;
	QUIET_GEN     = @echo '    ' GEN $@;
	QUIET_AR      = @echo '    ' AR $@;
	QUIET_SUBDIR0 = +@subdir=
	QUIET_SUBDIR1 = ; echo '    ' SUBDIR  $$subdir; \
		        ${MAKE} -C $$subdir
endif

.PHONY: check clean distclean


check: 
	@command -v python3 >/dev/null 2>&1 || { echo >&2 "python3 is required for 'make check', but is not in your PATH. Aborting."; exit 1; }
	@echo 
	@echo Running Infernal test suite...
	@echo
	@${srcdir}/../${ESLDIR}/devkit/sqc ${SQCLEVEL} ${srcdir}/testsuite.sqc .. ${srcdir}/..

devcheck: 
	@command -v python3 >/dev/null 2>&1 || { echo >&2 "python3 is required for 'make devcheck', but is not in your PATH. Aborting."; exit 1; }
	@echo 
	@echo Running Infernal developer test suite...
	@echo
	@${srcdir}/../${ESLDIR}/devkit/sqc ${DEVSQCLEVEL} ${srcdir}/dev_testsuite.sqc .. ${srcdir}/..

cmalign:
	@echo 
	@echo Running Infernal cmalign developer test suite...
	@echo
	@grep align: ${srcdir}/dev_testsuite.sqc > ${srcdir}/tmp.sqc
	@${srcdir}/../${ESLDIR}/devkit/sqc ${DEVSQCLEVEL} ${srcdir}/tmp.sqc     .. ${srcdir}/..
	@rm ${srcdir}/tmp.sqc

cmbuild:
	@echo 
	@echo Running Infernal cmbuild developer test suite...
	@echo
	@grep build: ${srcdir}/dev_testsuite.sqc > ${srcdir}/tmp.sqc
	@${srcdir}/../${ESLDIR}/devkit/sqc ${DEVSQCLEVEL} ${srcdir}/tmp.sqc     .. ${srcdir}/..
	@rm ${srcdir}/tmp.sqc

cmcalibrate:
	@echo 
	@echo Running Infernal cmcalibrate developer test suite...
	@echo
	@grep calibrate: ${srcdir}/dev_testsuite.sqc > ${srcdir}/tmp.sqc
	@${srcdir}/../${ESLDIR}/devkit/sqc ${DEVSQCLEVEL} ${srcdir}/tmp.sqc     .. ${srcdir}/..
	@rm ${srcdir}/tmp.sqc

cmconvert:
	@echo 
	@echo Running Infernal cmconvert developer test suite...
	@echo
	@grep convert: ${srcdir}/dev_testsuite.sqc > ${srcdir}/tmp.sqc
	@${srcdir}/../${ESLDIR}/devkit/sqc ${DEVSQCLEVEL} ${srcdir}/tmp.sqc     .. ${srcdir}/..
	@rm ${srcdir}/tmp.sqc

cmemit:
	@echo 
	@echo Running Infernal cmemit developer test suite...
	@echo
	@grep emit: ${srcdir}/dev_testsuite.sqc > ${srcdir}/tmp.sqc
	@${srcdir}/../${ESLDIR}/devkit/sqc ${DEVSQCLEVEL} ${srcdir}/tmp.sqc     .. ${srcdir}/..
	@rm ${srcdir}/tmp.sqc

cmfetch:
	@echo 
	@echo Running Infernal cmfetch developer test suite...
	@echo
	@grep fetch: ${srcdir}/dev_testsuite.sqc > ${srcdir}/tmp.sqc
	@${srcdir}/../${ESLDIR}/devkit/sqc ${DEVSQCLEVEL} ${srcdir}/tmp.sqc     .. ${srcdir}/..
	@rm ${srcdir}/tmp.sqc

cmpress:
	@echo 
	@echo Running Infernal cmpress developer test suite...
	@echo
	@grep press: ${srcdir}/dev_testsuite.sqc > ${srcdir}/tmp.sqc
	@${srcdir}/../${ESLDIR}/devkit/sqc ${DEVSQCLEVEL} ${srcdir}/tmp.sqc     .. ${srcdir}/..
	@rm ${srcdir}/tmp.sqc

cmscan:
	@echo 
	@echo Running Infernal cmscan developer test suite...
	@echo
	@grep scan: ${srcdir}/dev_testsuite.sqc > ${srcdir}/tmp.sqc
	@${srcdir}/../${ESLDIR}/devkit/sqc ${DEVSQCLEVEL} ${srcdir}/tmp.sqc     .. ${srcdir}/..
	@rm ${srcdir}/tmp.sqc

cmsearch:
	@echo 
	@echo Running Infernal cmsearch developer test suite...
	@echo
	@grep search: ${srcdir}/dev_testsuite.sqc > ${srcdir}/tmp.sqc
	@${srcdir}/../${ESLDIR}/devkit/sqc ${DEVSQCLEVEL} ${srcdir}/tmp.sqc     .. ${srcdir}/..
	@rm ${srcdir}/tmp.sqc

cmstat:
	@echo 
	@echo Running Infernal cmstat developer test suite...
	@echo
	@grep stat: ${srcdir}/dev_testsuite.sqc > ${srcdir}/tmp.sqc
	@${srcdir}/../${ESLDIR}/devkit/sqc ${DEVSQCLEVEL} ${srcdir}/tmp.sqc     .. ${srcdir}/..
	@rm ${srcdir}/tmp.sqc

clean:
	-rm -f *.gcno *.gcda *.gcov
	-rm -f *.o *~ 
ifndef V
	@echo '     ' CLEAN testsuite
endif

distclean: clean
	-rm -f Makefile
	-rm -f esltmp[a-z][a-z][0-9][0-9]*


