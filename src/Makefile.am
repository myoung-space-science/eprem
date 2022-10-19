noinst_LIBRARIES = libeprem.a

libeprem_a_OTHER = 

libeprem_a_SOURCES = \
  baseTypes.c \
  configuration.c \
  cubeShellInit.c \
  cubeShellStruct.c \
  energeticParticlesBoundary.c \
  energeticParticles.c \
  energeticParticlesInit.c \
  energeticParticlesTypes.c \
  error.c \
  flow.c \
  geometry.c \
  global.c \
  mpiInit.c \
  nearestGridNode.c \
  observerOutput.c \
  search.c \
  simCore.c \
  unifiedOutput.c \
  $(libeprem_a_OTHER)

bin_PROGRAMS = eprem

eprem_SOURCES = eprem.c 

eprem_DEPENDENCIES = libeprem.a
eprem_LDADD = $(eprem_DEPENDENCIES)
nodist_EXTRA_eprem_SOURCES = dummy.cxx

