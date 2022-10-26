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
  src/unifiedOutput.c

distclean-local:
	find . -type f -name '*~' -exec rm {} \;

install:
	mv src/eprem bin/eprem
