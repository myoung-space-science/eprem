/*-----------------------------------------------
 -- EMMREM: configuration.h
 ------------------------------------------------*/

/* The Earth-Moon-Mars Radiation Environment Module (EMMREM) software is */
/* free software; you can redistribute and/or modify the EMMREM sotware */
/* or any part of the EMMREM software under the terms of the GNU General */
/* Public License (GPL) as published by the Free Software Foundation; */
/* either version 2 of the License, or (at your option) any later */
/* version. Software that uses any portion of the EMMREM software must */
/* also be released under the GNU GPL license (version 2 of the GNU GPL */
/* license or a later version). A copy of this GNU General Public License */
/* may be obtained by writing to the Free Software Foundation, Inc., 59 */
/* Temple Place, Suite 330, Boston MA 02111-1307 USA or by viewing the */
/* license online at http://www.gnu.org/copyleft/gpl.html. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libconfig.h>
#include "global.h"
#include "configuration.h"
#include "mpiInit.h"
#include "error.h"

Config_t config;
config_t cfg;

const double third = 1.0/3.0;

void
initGlobalParameters( char* configFilename )
{

  getParams(configFilename);
  checkParams();
  setRuntimeConstants();

}


void
getParams( char* configFilename)
{

  // initialize the config parser structure and read external config file
  config_init(&cfg);

  if(! config_read_file(&cfg, configFilename)) {

    if (mpi_rank == 0)
      printf("\n\nError on line %i of %s: %s\n\n", config_error_line(&cfg), configFilename, config_error_text(&cfg));

    panic("Unable to read configuration file.");

  }

  config.useDegrees = readInt("useDegrees", 0, 0, 1);

  config.numNodesPerStream = readInt("numNodesPerStream",N_PROCS,N_PROCS,LARGEINT);
  config.numRowsPerFace = readInt("numRowsPerFace", 2, 1, LARGEINT);
  config.numColumnsPerFace = readInt("numColumnsPerFace", 2, 1, LARGEINT);
  config.numEnergySteps = readInt("numEnergySteps", 20, 2, LARGEINT);
  config.numMuSteps = readInt("numMuSteps", 20, 2, LARGEINT);

  config.rScale = readDouble("rScale", RSAU, SMALLFLOAT, LARGEFLOAT);
  config.flowMag = readDouble("flowMag", 400.0e5, SMALLFLOAT, LARGEFLOAT);
  config.mhdDensityAu = readDouble("mhdDensityAu", 8.30, SMALLFLOAT, LARGEFLOAT);
  config.mhdBAu = readDouble("mhdBAu", 1.60e-5, SMALLFLOAT, LARGEFLOAT);
  config.simStartTime = readDouble("simStartTime", 0.0, 0.0, LARGEFLOAT);
  config.tDel = readDouble("tDel", 0.01041666666667, SMALLFLOAT, LARGEFLOAT);
  config.simStopTime = readDouble("simStopTime", config.simStartTime + config.tDel, config.simStartTime, LARGEFLOAT);
  config.numEpSteps = readInt("numEpSteps", 30, 1, LARGEINT);
  config.aziSunStart = readDouble("aziSunStart", 0.0, 0.0, LARGEFLOAT);
  config.omegaSun = readDouble("omegaSun", 0.001429813, 0.0, LARGEFLOAT);
  config.lamo = readDouble("lamo", 1.0, SMALLFLOAT, LARGEFLOAT);
  config.dsh_min = readDouble("dsh_min", 5.0e-5, SMALLFLOAT, LARGEFLOAT);
  config.dsh_hel_min = readDouble("dsh_hel_min", 2.5e-4, SMALLFLOAT, LARGEFLOAT);
  config.kperxkpar = readDouble("kperxkpar", 0.01, SMALLFLOAT, LARGEFLOAT);
  config.mfpRadialPower = readDouble("mfpRadialPower", 2.0, -1.0 * LARGEFLOAT, LARGEFLOAT);
  config.rigidityPower = readDouble("rigidityPower", third, SMALLFLOAT, LARGEFLOAT);
  config.focusingLimit = readDouble("focusingLimit", 1.0, 0.0, 1.0);

  config.eMin = readDouble("eMin", 1.0, SMALLFLOAT, LARGEFLOAT);
  config.eMax = readDouble("eMax", 1000.0, config.eMin, LARGEFLOAT);
  config.useStochastic = readInt("useStochastic", 0, 0, 1);
  config.useEPBoundary = readInt("useEPBoundary", 1, 0, 1);
  config.checkSeedPopulation = readInt("checkSeedPopulation", 1, 0, 1);

  config.seedFunctionTest = readInt("seedFunctionTest", 0, 0, 1);

  config.fluxLimiter = readInt("fluxLimiter", 1, 0, 1);

  config.gammaEhigh = readDouble("gammaEhigh", 0.0, -1.0 * LARGEFLOAT, LARGEFLOAT);
  config.gammaElow = readDouble("gammaElow", 0.0, -1.0 * LARGEFLOAT, LARGEFLOAT);

  config.FailModeDump = readInt("FailModeDump", 1, 0, 1);

  config.outputFloat = readInt("outputFloat", 0, 0, 1);

  config.streamLegacyPrefix = readInt("streamLegacyPrefix", 0, 0, 1);
  config.pointLegacyPrefix = readInt("pointLegacyPrefix", 0, 0, 1);

  config.unifiedOutput = readInt("unifiedOutput", 1, 0, 1);
  config.unifiedOutputTime = readDouble("unifiedOutputTime", 0.0, 0.0, LARGEFLOAT);

  config.streamFluxOutput = readInt("streamFluxOutput", 0, 0, 1);
  config.streamFluxOutputTime = readDouble("streamFluxOutputTime", 0.0, 0.0, LARGEFLOAT);

  config.subTimeCouple = readInt("subTimeCouple", 0, 0, 1);

  config.epremDomain = readInt("epremDomain", 0, 0, 1);
  config.epremDomainOutputTime = readDouble("epremDomainOutputTime", 0.0, 0.0, LARGEFLOAT);

  config.unstructuredDomain = readInt("unstructuredDomain", 0, 0, 1);
  config.unstructuredDomainOutputTime = readDouble("unstructuredDomainOutputTime", 0.0, 0.0, LARGEFLOAT);

  config.useAdiabaticChange = readInt("useAdiabaticChange", 1, 0, 1);
  config.useAdiabaticFocus = readInt("useAdiabaticFocus", 1, 0, 1);
  config.useShellDiffusion = readInt("useShellDiffusion", 0, 0, 1);
  config.useParallelDiffusion = readInt("useParallelDiffusion", 1, 0, 1);
  config.useDrift = readInt("useDrift", 0, 0, 1);

  config.numSpecies = readInt("numSpecies", 1, 1, 100);
  Scalar_t defaultMass[1] = {1.0};
  if (config_lookup(&cfg, "mass") == NULL) {
    config.mass = readDoubleArray("mass", 1, 0, defaultMass, 1.0, LARGEFLOAT);
  } else {
    config.mass = readDoubleArray("mass", 1, config.numSpecies, defaultMass, 1.0, LARGEFLOAT);
  }
  Scalar_t defaultCharge[1] = {1.0};
  if (config_lookup(&cfg, "charge") == NULL) {
    config.charge = readDoubleArray("charge", 1, 0, defaultCharge, 1.0, LARGEFLOAT);
  } else {
    config.charge = readDoubleArray("charge", 1, config.numSpecies, defaultCharge, 1.0, LARGEFLOAT);
  }

  config.pointObserverOutput = readInt("pointObserverOutput", 0, 0, 1);
  config.pointObserverOutputTime = readDouble("pointObserverOutputTime", 0.0, 0.0, LARGEFLOAT);

  if (config.pointObserverOutput == 1) {
    config.numObservers = readInt("numObservers", 1, 0, 1000);
  } else {
    config.numObservers = readInt("numObservers", 0, 0, 1000);
  }
  if (config.numObservers > 0) {
    Scalar_t defaultObsR[1] = {config.rScale};
    Scalar_t defaultObsTheta[1] = {0.0};
    Scalar_t defaultObsPhi[1] = {0.0};
    Scalar_t *thetaArr, *phiArr;
    config.obsR = readDoubleArray("obsR", 1, config.numObservers, defaultObsR, config.rScale, LARGEFLOAT);
    if (config.useDegrees == 1) {
      thetaArr = readDoubleArray("obsTheta", 1, config.numObservers, defaultObsTheta, 0.0, RAD2DEG*PI);
      phiArr = readDoubleArray("obsPhi", 1, config.numObservers, defaultObsPhi, 0.0, RAD2DEG*TWO_PI);
    } else {
      thetaArr = readDoubleArray("obsTheta", 1, config.numObservers, defaultObsTheta, 0.0, PI);
      phiArr = readDoubleArray("obsPhi", 1, config.numObservers, defaultObsPhi, 0.0, TWO_PI);
    }
    if (config.useDegrees == 1) {
      config.obsTheta = (Scalar_t *)malloc(sizeof(double) * config.numObservers);
      config.obsPhi = (Scalar_t *)malloc(sizeof(double) * config.numObservers);
      for (int i=0; i<config.numObservers; i++) {
        config.obsTheta[i] = DEG2RAD*thetaArr[i];
        config.obsPhi[i]   = DEG2RAD*phiArr[i];
      }
    } else {
      config.obsTheta = thetaArr;
      config.obsPhi   = phiArr;
    }
  }

  config.idw_p = readDouble("idw_p", 3.0, SMALLFLOAT, LARGEFLOAT);

  config.mhdCouple = readInt("mhdCouple", 0, 0, 1);
  config.mhdNumFiles = readInt("mhdNumFiles", 0, 0, 32767);
  config.useMhdSteadyStateDt = readInt("useMhdSteadyStateDt", 1, 0, 1);
  config.mhdSteadyState = readInt("mhdSteadyState", 1, 0, 1);
  config.mhdDirectory = (char*)readString("mhdDirectory"," ");
  config.mhdDigits = readInt("mhdDigits", 3, 0, 32767);

  config.mhdCoupledTime = readInt("mhdCoupledTime", 1, 0, 1);
  config.mhdStartTime = readDouble("mhdStartTime", 0.0, 0.0, LARGEFLOAT);
  config.epEquilibriumCalcDuration = readDouble("epEquilibriumCalcDuration", 0.0, 0.0, LARGEFLOAT);
  config.preEruptionDuration = readDouble("preEruptionDuration", 0.0, 0.0, LARGEFLOAT);

  config.mhdRadialMin = readDouble("mhdRadialMin", 0.0, 0.0, LARGEFLOAT);
  config.mhdRadialMax = readDouble("mhdRadialMax", 0.0, 0.0, LARGEFLOAT);
  config.mhdVmin = readDouble("mhdVmin", 50.0e5, 0.0, LARGEFLOAT);

  config.mhdInitFromOuterBoundary = readInt("mhdInitFromOuterBoundary", 2, 0, 2);
  config.mhdInitMonteCarlo = readInt("mhdInitMonteCarlo", 0, 0, 1);
  config.mhdInitRadius = readDouble("mhdInitRadius", 0.0, 0.0, LARGEFLOAT);
  config.mhdInitTimeStep = readDouble("mhdInitTimeStep", 0.000011574074074, 0.0, LARGEFLOAT);

  config.useManualStreamSpawnLoc = readInt("useManualStreamSpawnLoc",0,0,1);
  Scalar_t defaultPos[1] = {0.0};
  if (config.useManualStreamSpawnLoc > 0){
    config.streamSpawnLocAzi = readDoubleArray("streamSpawnLocAzi", 1, 6*config.numRowsPerFace*config.numColumnsPerFace, defaultPos, 0.0, TWO_PI);
    config.streamSpawnLocZen = readDoubleArray("streamSpawnLocZen", 1, 6*config.numRowsPerFace*config.numColumnsPerFace, defaultPos, 0.0, PI);
  }

  config.parallelFlow = readDouble("parallelFlow", 0.0, 0.0, LARGEFLOAT);
  config.fieldAligned = readInt("fieldAligned", 0, 0, 1);

  config.epCalcStartTime = readDouble("epCalcStartTime", config.simStartTime, 0.0, LARGEFLOAT);

  config.mhdRotateSolution = readInt("mhdRotateSolution", 1, 0, 1);

  config.mhdBConvert = readDouble("mhdBConvert", 1.0, 0.0, LARGEFLOAT);
  config.mhdVConvert = readDouble("mhdVConvert", 1.0, 0.0, LARGEFLOAT);
  config.mhdRhoConvert = readDouble("mhdRhoConvert", 1.0, 0.0, LARGEFLOAT);
  config.mhdTimeConvert = readDouble("mhdTimeConvert", 1.0, 0.0, LARGEFLOAT);

  config.useBoundaryFunction = readInt("useBoundaryFunction", 1, 0, 1);
  config.boundaryFunctionInitDomain = readInt("boundaryFunctionInitDomain", 1, 0, 1);

  config.boundaryFunctAmplitude = readDouble("boundaryFunctAmplitude", 1.0, SMALLFLOAT, LARGEFLOAT);
  config.boundaryFunctXi = readDouble("boundaryFunctXi", 1.0, 0.0, LARGEFLOAT);
  config.boundaryFunctBeta = readDouble("boundaryFunctBeta", 2.0, 0.0, LARGEFLOAT);
  config.boundaryFunctR0 = readDouble("boundaryFunctR0", 1.0, config.rScale, LARGEFLOAT);
  config.boundaryFunctGamma = readDouble("boundaryFunctGamma", 2.0, 0.0, LARGEFLOAT);
  config.boundaryFunctEr = readDouble("boundaryFunctEr", 1.0, 0.0, LARGEFLOAT);
  config.boundaryFunctEcutoff = readDouble("boundaryFunctEcutoff", 1.0, 0.0, LARGEFLOAT);

  config.shockSolver = readInt("shockSolver", 0, 0, 1);
  config.shockDetectPercent = readDouble("shockDetectPercent", 1.0, 0.0, LARGEFLOAT);
  config.minInjectionEnergy = readDouble("minInjectionEnergy", 0.01, SMALLFLOAT, LARGEFLOAT);
  config.shockInjectionFactor = readDouble("shockInjectionFactor", 1.0, 0.0, LARGEFLOAT);

  config.idealShock = readInt("idealShock", 0, 0, 1);
  config.idealShockSharpness = readDouble("idealShockSharpness", 1.0, SMALLFLOAT, LARGEFLOAT);
  config.idealShockScaleLength = readDouble("idealShockScaleLength", 0.0046491, SMALLFLOAT, LARGEFLOAT);
  config.idealShockScale = readDouble("idealShockScale", config.idealShockSharpness / config.idealShockScaleLength, 0.0, LARGEFLOAT);
  config.idealShockJump = readDouble("idealShockJump", 4.0, SMALLFLOAT, LARGEFLOAT);
  config.idealShockFalloff = readDouble("idealShockFalloff", 0.0, 0.0, LARGEFLOAT);
  config.idealShockSpeed = readDouble("idealShockSpeed", 1500e5, SMALLFLOAT, LARGEFLOAT);
  config.idealShockInitTime = readDouble("idealShockInitTime", config.simStartTime, config.simStartTime, LARGEFLOAT);
  if (config.useDegrees == 1) {
    config.idealShockTheta = DEG2RAD * readDouble("idealShockTheta", 90.0, 0.0, 180.0);
    config.idealShockPhi = DEG2RAD * readDouble("idealShockPhi", 0.0, 0.0, 360.0);
    config.idealShockWidth = DEG2RAD * readDouble("idealShockWidth", 0.0, 0.0, 180.0);
    config.idealShockThetaWidth = DEG2RAD * readDouble("idealShockThetaWidth", RAD2DEG * config.idealShockWidth, 0.0, 180.0);
    config.idealShockPhiWidth = DEG2RAD * readDouble("idealShockPhiWidth", RAD2DEG * config.idealShockWidth, 0.0, 180.0);
  } else {
    config.idealShockTheta = readDouble("idealShockTheta", HALF_PI, 0.0, PI);
    config.idealShockPhi = readDouble("idealShockPhi", 0.0, 0.0, TWO_PI);
    config.idealShockWidth = readDouble("idealShockWidth", 0.0, 0.0, PI);
    config.idealShockThetaWidth = readDouble("idealShockThetaWidth", config.idealShockWidth, 0.0, PI);
    config.idealShockPhiWidth = readDouble("idealShockPhiWidth", config.idealShockWidth, 0.0, PI);
  }

  config.dumpFreq = readInt("dumpFreq",1, 0, 1000000);
  config.outputRestart = readInt("outputRestart",0, 0, 1000000);
  config.dumpOnAbort = readInt("dumpOnAbort",0, 0, 1);
  config.saveRestartFile = readInt("saveRestartFile",0, 0, 1);

  config.warningsFile = (char*)readString("warningsFile", "warningsXXX.txt");

  config.adiabaticChangeAlg = readInt("adiabaticChangeAlg", 1, 1, 3);
  config.adiabaticFocusAlg = readInt("adiabaticFocusAlg", 1, 1, 3);

}


Index_t readInt(char *key, Index_t defaultVal, Index_t minVal, Index_t maxVal) {

  Index_t val;

  if (! config_lookup_int(&cfg, key, &val) )
    val = defaultVal;

  checkIntBounds(key, val, minVal, maxVal);

  if (mpi_rank == 0)
    printf("%s: %i\n", key, (Index_t)val);

  return (Index_t)val;

}


Scalar_t readDouble(char *key, Scalar_t defaultVal, Scalar_t minVal, Scalar_t maxVal) {

  Scalar_t val;

  if (! config_lookup_float(&cfg, key, &val) )
    val = defaultVal;

  checkDoubleBounds(key, val, minVal, maxVal);

  if (mpi_rank == 0)
    printf("%s: %.4e\n", key, (Scalar_t)val);

  return (Scalar_t)val;

}


const char *readString(char *key, char *defaultVal) {

  const char *val;

  if (! config_lookup_string(&cfg, key, &val) )
    val = defaultVal;

  if (mpi_rank == 0)
    printf("%s: %s\n", key, val);

  return val;

}


Scalar_t *readDoubleArray(char *key, int defaultSize, int size, Scalar_t *defaultVal, Scalar_t minVal, Scalar_t maxVal) {

  Index_t i;
  Scalar_t *val;
  const config_setting_t *Arr;

  if (mpi_rank == 0) {
    printf("%s: ", key);
  }

  if (size > 0) {

    val = (Scalar_t *)malloc(sizeof(double) * size);
    Arr = config_lookup(&cfg, key);

    for (i = 0; i < size; i++) {
      val[i] = config_setting_get_float_elem(Arr, i);
      checkDoubleBounds(key, val[i], minVal, maxVal);
    }

    if (mpi_rank == 0) {
      printf("[");
      if (size == 1) {
        printf("%.4e]\n", val[0]);
      } else {
        for (i = 0; i < size-1; i++) {
          printf("%.4e, ", val[i]);
        }
        printf("%.4e]\n", val[i]);
      }
    }

    return val;

  } else {

    if (mpi_rank == 0){
      printf("[");
      if (defaultSize == 1) {
        printf("%.4e]\n", defaultVal[0]);
      } else {
        for (i = 0; i < defaultSize-1; i++) {
          printf("%.4e, ", defaultVal[i]);
        }
        printf("%.4e]\n", defaultVal[i]);
      }
    }

    return defaultVal;

  }

}


void
checkParams( void )
{
  // Enforce bounds on shock angles in radians.
  checkDoubleBounds("idealShockTheta", config.idealShockTheta, 0.0, PI);
  checkDoubleBounds("idealShockPhi", config.idealShockTheta, 0.0, 2.0 * PI);
  checkDoubleBounds("idealShockWidth", config.idealShockTheta, 0.0, PI);
  // Enforce bounds on observer angles in radians.
  for (int i=0; i<config.numObservers; i++) {
    checkDoubleBounds("obsTheta", config.obsTheta[i], 0.0, PI);
    checkDoubleBounds("obsPhi", config.obsTheta[i], 0.0, 2.0 * PI);
  }
}


void
checkIntBounds(char *key, Index_t val, Index_t minVal, Index_t maxVal)
{
  if ( (val < minVal) || (val > maxVal) ) {

    if (mpi_rank == 0) {
      printf("%s=%d is out of the acceptable range: [%d, %d]\n", key, (Index_t)val, minVal, maxVal);
      panic("the configuration reader detected an invalid value.\n");
    }

  }
}


void
checkDoubleBounds(char *key, Scalar_t val, Scalar_t minVal, Scalar_t maxVal)
{
  if ( (val < minVal) || (val > maxVal) ) {

    if (mpi_rank == 0) {
      printf("%s=%.4e is out of the acceptable range: [%.4e, %.4e]\n", key, (Scalar_t)val, minVal, maxVal);
      panic("the configuration reader detected an invalid value.\n");
    }

  }
}


void
setRuntimeConstants( void )
{

  FACE_ROWS = config.numRowsPerFace;
  FACE_COLS = config.numColumnsPerFace;
  NUM_SPECIES = config.numSpecies;
  NUM_ESTEPS = config.numEnergySteps;
  NUM_MUSTEPS = config.numMuSteps;
  NUM_OBS = config.numObservers;
  AdiabaticChangeAlg = config.adiabaticChangeAlg;
  AdiabaticFocusAlg = config.adiabaticFocusAlg;

  TOTAL_NUM_SHELLS = config.numNodesPerStream;

  config.simStartTimeDay = config.simStartTime / DAY;
  config.simStopTimeDay  = config.simStopTime / DAY;
  config.tDel            /= DAY;
  
  config.mhdUs           = ( config.flowMag / C );
  config.mhdNsAu         = ( config.mhdDensityAu / MHD_DENSITY_NORM );
  config.mhdBsAu         = ( config.mhdBAu / MHD_B_NORM );

}
