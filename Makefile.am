bin_PROGRAMS = eprem

eprem_SOURCES = \
src/baseTypes.c \
src/configuration.c \
src/cubeShellInit.c \
src/cubeShellStruct.c \
src/energeticParticlesBoundary.c \
src/energeticParticles.c \
src/energeticParticlesInit.c \
src/energeticParticlesTypes.c \
src/eprem.c \
src/error.c \
src/flow.c \
src/geometry.c \
src/global.c \
src/mpiInit.c \
src/nearestGridNode.c \
src/observerOutput.c \
src/search.c \
src/simCore.c \
src/unifiedOutput.c \
src/baseTypes.h \
src/configuration.h \
src/cubeShellInit.h \
src/cubeShellStruct.h \
src/energeticParticlesBoundary.h \
src/energeticParticles.h \
src/energeticParticlesInit.h \
src/energeticParticlesTypes.h \
src/error.h \
src/flow.h \
src/geometry.h \
src/global.h \
src/mpiInit.h \
src/nearestGridNode.h \
src/observerOutput.h \
src/safeNetcdf.h \
src/searchTypes.h \
src/simCore.h \
src/timers.h \
src/unifiedOutput.h

EXTRA_DIST = setup.sh

distclean-local:
	find . -type f -name '*~' -exec rm {} \;

maintainer-clean-local:
	rm -f autoscan.log
	rm -f configure.scan

install-exec-hook:
	cp -u config.log ${prefix} &> /dev/null
	mkdir -p bin
	mv -f eprem bin/eprem-latest
