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

#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <mfhdf.h>
#include "readMHD.h"
#include "mpiInit.h"
#include "global.h"
#include "configuration.h"
#include "error.h"
#include "simCore.h"
#include "flow.h"
#include "observerOutput.h"
#include "timers.h"
#include "mhdIO.h"

Scalar_t *mhdTime;
Scalar_t *mhdHelTime;

Scalar_t phiOffset;
Scalar_t phiHelOffset;

Index_t mhdFileIndex0;
Index_t mhdFileIndex1;
Index_t mhdHelFileIndex0;
Index_t mhdHelFileIndex1;

Index_t mhdFileIndex_loaded0=-9999;
Index_t mhdFileIndex_loaded1=-9999;
Index_t mhdHelFileIndex_loaded0=-9999;
Index_t mhdHelFileIndex_loaded1=-9999;

Index_t mhdMallocFlag;
Index_t mhdEqFileFlag;

int32_t mhdDimMin[1] = {0};

int32_t mhdBppDimMax[1];
int32_t mhdBptDimMax[1];
int32_t mhdBprDimMax[1];
int32_t mhdHelBppDimMax[1];
int32_t mhdHelBptDimMax[1];
int32_t mhdHelBprDimMax[1];

int32_t mhdBtpDimMax[1];
int32_t mhdBttDimMax[1];
int32_t mhdBtrDimMax[1];
int32_t mhdHelBtpDimMax[1];
int32_t mhdHelBttDimMax[1];
int32_t mhdHelBtrDimMax[1];

int32_t mhdBrpDimMax[1];
int32_t mhdBrtDimMax[1];
int32_t mhdBrrDimMax[1];
int32_t mhdHelBrpDimMax[1];
int32_t mhdHelBrtDimMax[1];
int32_t mhdHelBrrDimMax[1];

int32_t mhdVppDimMax[1];
int32_t mhdVptDimMax[1];
int32_t mhdVprDimMax[1];
int32_t mhdHelVppDimMax[1];
int32_t mhdHelVptDimMax[1];
int32_t mhdHelVprDimMax[1];

int32_t mhdVtpDimMax[1];
int32_t mhdVttDimMax[1];
int32_t mhdVtrDimMax[1];
int32_t mhdHelVtpDimMax[1];
int32_t mhdHelVttDimMax[1];
int32_t mhdHelVtrDimMax[1];

int32_t mhdVrpDimMax[1];
int32_t mhdVrtDimMax[1];
int32_t mhdVrrDimMax[1];
int32_t mhdHelVrpDimMax[1];
int32_t mhdHelVrtDimMax[1];
int32_t mhdHelVrrDimMax[1];

int32_t mhdDpDimMax[1];
int32_t mhdDtDimMax[1];
int32_t mhdDrDimMax[1];
int32_t mhdHelDpDimMax[1];
int32_t mhdHelDtDimMax[1];
int32_t mhdHelDrDimMax[1];

float * mhdBppDim;
float * mhdBptDim;
float * mhdBprDim;
float * mhdHelBppDim;
float * mhdHelBptDim;
float * mhdHelBprDim;

float * mhdBtpDim;
float * mhdBttDim;
float * mhdBtrDim;
float * mhdHelBtpDim;
float * mhdHelBttDim;
float * mhdHelBtrDim;

float * mhdBrpDim;
float * mhdBrtDim;
float * mhdBrrDim;
float * mhdHelBrpDim;
float * mhdHelBrtDim;
float * mhdHelBrrDim;

float * mhdVppDim;
float * mhdVptDim;
float * mhdVprDim;
float * mhdHelVppDim;
float * mhdHelVptDim;
float * mhdHelVprDim;

float * mhdVtpDim;
float * mhdVttDim;
float * mhdVtrDim;
float * mhdHelVtpDim;
float * mhdHelVttDim;
float * mhdHelVtrDim;

float * mhdVrpDim;
float * mhdVrtDim;
float * mhdVrrDim;
float * mhdHelVrpDim;
float * mhdHelVrtDim;
float * mhdHelVrrDim;

float * mhdDpDim;
float * mhdDtDim;
float * mhdDrDim;
float * mhdHelDpDim;
float * mhdHelDtDim;
float * mhdHelDrDim;


float * mhdBp_0;
float * mhdBt_0;
float * mhdBr_0;
float * mhdVp_0;
float * mhdVt_0;
float * mhdVr_0;
float * mhdD_0;
float * mhdHelBp_0;
float * mhdHelBt_0;
float * mhdHelBr_0;
float * mhdHelVp_0;
float * mhdHelVt_0;
float * mhdHelVr_0;
float * mhdHelD_0;

float * mhdBp_1;
float * mhdBt_1;
float * mhdBr_1;
float * mhdVp_1;
float * mhdVt_1;
float * mhdVr_1;
float * mhdD_1;
float * mhdHelBp_1;
float * mhdHelBt_1;
float * mhdHelBr_1;
float * mhdHelVp_1;
float * mhdHelVt_1;
float * mhdHelVr_1;
float * mhdHelD_1;

MPI_Win mhdBp_0_win;
MPI_Win mhdBt_0_win;
MPI_Win mhdBr_0_win;
MPI_Win mhdVp_0_win;
MPI_Win mhdVt_0_win;
MPI_Win mhdVr_0_win;
MPI_Win mhdD_0_win;
MPI_Win mhdHelBp_0_win;
MPI_Win mhdHelBt_0_win;
MPI_Win mhdHelBr_0_win;
MPI_Win mhdHelVp_0_win;
MPI_Win mhdHelVt_0_win;
MPI_Win mhdHelVr_0_win;
MPI_Win mhdHelD_0_win;

MPI_Win mhdBp_1_win;
MPI_Win mhdBt_1_win;
MPI_Win mhdBr_1_win;
MPI_Win mhdVp_1_win;
MPI_Win mhdVt_1_win;
MPI_Win mhdVr_1_win;
MPI_Win mhdD_1_win;
MPI_Win mhdHelBp_1_win;
MPI_Win mhdHelBt_1_win;
MPI_Win mhdHelBr_1_win;
MPI_Win mhdHelVp_1_win;
MPI_Win mhdHelVt_1_win;
MPI_Win mhdHelVr_1_win;
MPI_Win mhdHelD_1_win;

int coupleStarted=0;
char file_extension[5];

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     ERR(intn status)                                     /*--*/
/*--                                                              --*/
/*--  This function checks for an error when loading an SD        --*/
/*--  element.                                                    --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/
  if ((status != 0) && (mpi_rank == 0)) {
    printf("\n\nERROR WITH AN SD! status = %d\n",status);
  }
}/*-------- END ERR(intn status)  ----------------------------------*/
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     mhdFetchCouplingInfo(void)                           /*--*/
/*--                                                              --*/
/*--  This function loads coupling info from the MHD directories  --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  FILE *rfile;

  int max = 10000;
  char line[10000];

  char *name = NULL;
  char *value = NULL;
  char delims[] = ": ";

  char mhdHelInfoFilename[MAX_STRING_SIZE] = "mhd_helio_run_info.txt";
  char mhdHelInfoFilenameWithPath[MAX_STRING_SIZE];

  // -- coronal coupling --//

  // -- heliospheric coupling --//
  if (config.mhdHelCouple > 0) {

    // build the path to the info file
    sprintf(mhdHelInfoFilenameWithPath, "%s%s", config.mhdHelDirectory, mhdHelInfoFilename);

    // attempt to open the info file
    rfile = fopen(mhdHelInfoFilenameWithPath, "r");
    if (rfile==NULL) {
      printf("ERROR - Could not open file \"%s\"\nReverting to defaults\n", mhdHelInfoFilenameWithPath);
    } else {
      if (mpi_rank == 0) {
        printf("Reading parameters from \"%s\"\n", mhdHelInfoFilenameWithPath);
      }
    }

    // read each line and check for known parameters
    // this assumes that each line contains `<name>: <value>` or `<name>:<value>`
    while (fgets(line, max, rfile) != NULL) {
      name = strtok(line, delims);
      if (name != NULL) {
        value = strtok(NULL, delims);
        if ((value != NULL) && (strcmp(name, "coronal_helio_phi_offset") == 0)) {
          phiHelOffset = (Scalar_t)atof(value);
        }
      }
    }
    if (mpi_rank == 0) {
      printf("phiHelOffset: %f\n", phiHelOffset);
    }

    // close the file
    fclose(rfile);

    name = NULL;
    value = NULL;

  }

}/*-------- END mhdFetchCouplingInfo()  ----------------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     mhdFetchFileList(void)                               /*--*/
/*--                                                              --*/
/*--  This function loads the list timesteps and filenames to use --*/
/*--  with MHD files.                                             --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  FILE * rfile;

  int max = 10000;
  char line[10000];

  int nFileLines;

  char *result = NULL;
  char delims[] = " \t\n\r";

  char mhdTimeFilename[MAX_STRING_SIZE] = "mhdTime.txt";
  char mhdTimeFilenameWithPath[MAX_STRING_SIZE];

  Scalar_t initialTime, timeVar;

  MPI_Aint N;
  MPI_Aint size;
  int disp_unit;

  // build the path to the time file
  sprintf(mhdTimeFilenameWithPath, "%s%s", config.mhdDirectory, mhdTimeFilename);

  // attempt to open the time file
  rfile = fopen(mhdTimeFilenameWithPath, "r");
  if (rfile==NULL) {
    printf("ERROR - Could not open file \"%s\"\n", mhdTimeFilenameWithPath);
    panic("Can't find MHD time list.");
  }

  // read the number of lines in the time file list
  nFileLines = 0;
  while (fgets(line, max, rfile) != NULL) {
    nFileLines++;
  }
  mhdTime = (Scalar_t *)malloc(sizeof(Scalar_t) * nFileLines);

  // reset the file pointer to the beginning of the file
  rewind(rfile);

  // read in the time file list
  for (int t=0; t<nFileLines; t++) {
    if (fgets(line, max, rfile) != NULL) {
      result = strtok(line, delims);
      if (result != NULL) {
        timeVar = (Scalar_t)atof(result) * config.mhdTimeConvert / DAY;
        if (t == 0) {initialTime = timeVar;}
        mhdTime[t] = config.mhdStartTime / DAY + (timeVar - initialTime);
      }
      result = strtok(NULL, delims);
    }
  }

  // close the file
  fclose(rfile);

  // store the number of files
  config.mhdNumFiles = nFileLines;

  result = NULL;

  // if coupling to the heliospheric domain
  if (config.mhdHelCouple > 0) {

    // build the path to the time file
    sprintf(mhdTimeFilenameWithPath, "%s%s", config.mhdHelDirectory, mhdTimeFilename);

    // attempt to open the time file
    rfile = fopen(mhdTimeFilenameWithPath, "r");
    if (rfile==NULL) {
      printf("ERROR - Could not open file \"%s\"\n", mhdTimeFilenameWithPath);
      panic("Can't find MHD Helio time list.");
    }

    // read the number of lines in the time file list
    nFileLines = 0;
    while (fgets(line, max, rfile) != NULL) {
      nFileLines++;
    }
    mhdHelTime = (Scalar_t *)malloc(sizeof(Scalar_t) * nFileLines);

    // reset the file pointer to the beginning of the file
    rewind(rfile);

    // read in the time file list
    for (int t=0; t<nFileLines; t++) {
      if (fgets(line, max, rfile) != NULL) {
        result = strtok(line, delims);
        if (result != NULL) {
          timeVar = (Scalar_t)atof(result) * config.mhdTimeConvert / DAY;
          if (t == 0) {initialTime = timeVar;}
          mhdHelTime[t] = config.mhdStartTime / DAY + (timeVar - initialTime);
        }
        result = strtok(NULL, delims);
      }
    }

    // close the file
    fclose(rfile);

    // store the number of files
    config.mhdHelNumFiles = nFileLines;

    result = NULL;

  }

  // Malloc arrays. (Only malloc once).
  if (mhdMallocFlag == 0) {

    // Set file type to hdf4.  RMC: Eventually this needs to either be
    // an input flag, or auto-detected (the current autodetection in
    // mhdIO.c does not work on some systems).  This will be moved from here
    // to a more logical place eventually.

    hdf5_input = 0;
    strncpy(file_extension,".hdf",strlen(".hdf")+1);

    // The size of the index arrays doesn't change in time in this version.
    mhdReadFieldIndex();
    // heliospheric coupling
    if (config.mhdHelCouple > 0)
      mhdHelReadFieldIndex();

    if(mpi_rank_shared==0){
      N=(int)mhdBprDimMax[0]*(int)mhdBptDimMax[0]*(int)mhdBppDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdBp_0, &mhdBp_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdBp_1, &mhdBp_1_win);
      N=(int)mhdBtrDimMax[0]*(int)mhdBttDimMax[0]*(int)mhdBtpDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdBt_0, &mhdBt_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdBt_1, &mhdBt_1_win);
      N=(int)mhdBrrDimMax[0]*(int)mhdBrtDimMax[0]*(int)mhdBrpDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdBr_0, &mhdBr_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdBr_1, &mhdBr_1_win);
      N=(int)mhdVprDimMax[0]*(int)mhdVptDimMax[0]*(int)mhdVppDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdVp_0, &mhdVp_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdVp_1, &mhdVp_1_win);
      N=(int)mhdVtrDimMax[0]*(int)mhdVttDimMax[0]*(int)mhdVtpDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdVt_0, &mhdVt_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdVt_1, &mhdVt_1_win);
      N=(int)mhdVrrDimMax[0]*(int)mhdVrtDimMax[0]*(int)mhdVrpDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdVr_0, &mhdVr_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdVr_1, &mhdVr_1_win);
      N=(int)mhdDrDimMax[0]*(int)mhdDtDimMax[0]*(int)mhdDpDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdD_0, &mhdD_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdD_1, &mhdD_1_win);

      // heliospheric coupling
      if (config.mhdHelCouple > 0) {

        N=(int)mhdHelBprDimMax[0]*(int)mhdHelBptDimMax[0]*(int)mhdHelBppDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelBp_0, &mhdHelBp_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelBp_1, &mhdHelBp_1_win);
        N=(int)mhdHelBtrDimMax[0]*(int)mhdHelBttDimMax[0]*(int)mhdHelBtpDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelBt_0, &mhdHelBt_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelBt_1, &mhdHelBt_1_win);
        N=(int)mhdHelBrrDimMax[0]*(int)mhdHelBrtDimMax[0]*(int)mhdHelBrpDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelBr_0, &mhdHelBr_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelBr_1, &mhdHelBr_1_win);
        N=(int)mhdHelVprDimMax[0]*(int)mhdHelVptDimMax[0]*(int)mhdHelVppDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelVp_0, &mhdHelVp_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelVp_1, &mhdHelVp_1_win);
        N=(int)mhdHelVtrDimMax[0]*(int)mhdHelVttDimMax[0]*(int)mhdHelVtpDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelVt_0, &mhdHelVt_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelVt_1, &mhdHelVt_1_win);
        N=(int)mhdHelVrrDimMax[0]*(int)mhdHelVrtDimMax[0]*(int)mhdHelVrpDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelVr_0, &mhdHelVr_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelVr_1, &mhdHelVr_1_win);
        N=(int)mhdHelDrDimMax[0]*(int)mhdHelDtDimMax[0]*(int)mhdHelDpDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelD_0, &mhdHelD_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelD_1, &mhdHelD_1_win);

      }

    }else{
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdBp_0, &mhdBp_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdBp_1, &mhdBp_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdBt_0, &mhdBt_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdBt_1, &mhdBt_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdBr_0, &mhdBr_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdBr_1, &mhdBr_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdVp_0, &mhdVp_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdVp_1, &mhdVp_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdVt_0, &mhdVt_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdVt_1, &mhdVt_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdVr_0, &mhdVr_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdVr_1, &mhdVr_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdD_0, &mhdD_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdD_1, &mhdD_1_win);
      MPI_Win_shared_query(mhdBp_0_win, 0, &size, &disp_unit, &mhdBp_0);
      MPI_Win_shared_query(mhdBt_0_win, 0, &size, &disp_unit, &mhdBt_0);
      MPI_Win_shared_query(mhdBr_0_win, 0, &size, &disp_unit, &mhdBr_0);
      MPI_Win_shared_query(mhdVp_0_win, 0, &size, &disp_unit, &mhdVp_0);
      MPI_Win_shared_query(mhdVt_0_win, 0, &size, &disp_unit, &mhdVt_0);
      MPI_Win_shared_query(mhdVr_0_win, 0, &size, &disp_unit, &mhdVr_0);
      MPI_Win_shared_query(mhdD_0_win, 0, &size, &disp_unit, &mhdD_0);
      MPI_Win_shared_query(mhdBp_1_win, 0, &size, &disp_unit, &mhdBp_1);
      MPI_Win_shared_query(mhdBt_1_win, 0, &size, &disp_unit, &mhdBt_1);
      MPI_Win_shared_query(mhdBr_1_win, 0, &size, &disp_unit, &mhdBr_1);
      MPI_Win_shared_query(mhdVp_1_win, 0, &size, &disp_unit, &mhdVp_1);
      MPI_Win_shared_query(mhdVt_1_win, 0, &size, &disp_unit, &mhdVt_1);
      MPI_Win_shared_query(mhdVr_1_win, 0, &size, &disp_unit, &mhdVr_1);
      MPI_Win_shared_query(mhdD_1_win, 0, &size, &disp_unit, &mhdD_1);

      // heliospheric coupling
      if (config.mhdHelCouple > 0) {

        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelBp_0, &mhdHelBp_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelBp_1, &mhdHelBp_1_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelBt_0, &mhdHelBt_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelBt_1, &mhdHelBt_1_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelBr_0, &mhdHelBr_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelBr_1, &mhdHelBr_1_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelVp_0, &mhdHelVp_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelVp_1, &mhdHelVp_1_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelVt_0, &mhdHelVt_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelVt_1, &mhdHelVt_1_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelVr_0, &mhdHelVr_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelVr_1, &mhdHelVr_1_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelD_0, &mhdHelD_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &mhdHelD_1, &mhdHelD_1_win);
        MPI_Win_shared_query(mhdHelBp_0_win, 0, &size, &disp_unit, &mhdHelBp_0);
        MPI_Win_shared_query(mhdHelBt_0_win, 0, &size, &disp_unit, &mhdHelBt_0);
        MPI_Win_shared_query(mhdHelBr_0_win, 0, &size, &disp_unit, &mhdHelBr_0);
        MPI_Win_shared_query(mhdHelVp_0_win, 0, &size, &disp_unit, &mhdHelVp_0);
        MPI_Win_shared_query(mhdHelVt_0_win, 0, &size, &disp_unit, &mhdHelVt_0);
        MPI_Win_shared_query(mhdHelVr_0_win, 0, &size, &disp_unit, &mhdHelVr_0);
        MPI_Win_shared_query(mhdHelD_0_win, 0, &size, &disp_unit, &mhdHelD_0);
        MPI_Win_shared_query(mhdHelBp_1_win, 0, &size, &disp_unit, &mhdHelBp_1);
        MPI_Win_shared_query(mhdHelBt_1_win, 0, &size, &disp_unit, &mhdHelBt_1);
        MPI_Win_shared_query(mhdHelBr_1_win, 0, &size, &disp_unit, &mhdHelBr_1);
        MPI_Win_shared_query(mhdHelVp_1_win, 0, &size, &disp_unit, &mhdHelVp_1);
        MPI_Win_shared_query(mhdHelVt_1_win, 0, &size, &disp_unit, &mhdHelVt_1);
        MPI_Win_shared_query(mhdHelVr_1_win, 0, &size, &disp_unit, &mhdHelVr_1);
        MPI_Win_shared_query(mhdHelD_1_win, 0, &size, &disp_unit, &mhdHelD_1);

      }

    }
    MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdBp_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdBt_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdBr_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdVp_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdVt_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdVr_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdD_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdBp_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdBt_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdBr_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdVp_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdVt_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdVr_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdD_1_win);

    // heliospheric coupling
    if (config.mhdHelCouple > 0) {

      MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdHelBp_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdHelBt_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdHelBr_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdHelVp_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdHelVt_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdHelVr_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdHelD_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdHelBp_1_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdHelBt_1_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdHelBr_1_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdHelVp_1_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdHelVt_1_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdHelVr_1_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, mhdHelD_1_win);

    }

    mhdMallocFlag = 1;

  }


}/*-------- END mhdFetchFileList()  --------------------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     cleanupMPIWindows(void)                              /*--*/
/*--                                                              --*/
/*--  This function cleans up MPI windows.                        --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

   if (config.mhdCouple > 0) {
      MPI_Win_unlock_all(mhdBp_0_win);
      MPI_Win_unlock_all(mhdBt_0_win);
      MPI_Win_unlock_all(mhdBr_0_win);
      MPI_Win_unlock_all(mhdVp_0_win);
      MPI_Win_unlock_all(mhdVt_0_win);
      MPI_Win_unlock_all(mhdVr_0_win);
      MPI_Win_unlock_all(mhdD_0_win);
      MPI_Win_unlock_all(mhdBp_1_win);
      MPI_Win_unlock_all(mhdBt_1_win);
      MPI_Win_unlock_all(mhdBr_1_win);
      MPI_Win_unlock_all(mhdVp_1_win);
      MPI_Win_unlock_all(mhdVt_1_win);
      MPI_Win_unlock_all(mhdVr_1_win);
      MPI_Win_unlock_all(mhdD_1_win);

      MPI_Win_free(&mhdBp_0_win);
      MPI_Win_free(&mhdBt_0_win);
      MPI_Win_free(&mhdBr_0_win);
      MPI_Win_free(&mhdVp_0_win);
      MPI_Win_free(&mhdVt_0_win);
      MPI_Win_free(&mhdVr_0_win);
      MPI_Win_free(&mhdD_0_win);
      MPI_Win_free(&mhdBp_1_win);
      MPI_Win_free(&mhdBt_1_win);
      MPI_Win_free(&mhdBr_1_win);
      MPI_Win_free(&mhdVp_1_win);
      MPI_Win_free(&mhdVt_1_win);
      MPI_Win_free(&mhdVr_1_win);
      MPI_Win_free(&mhdD_1_win);
  }

  // heliospheric coupling
  if (config.mhdHelCouple > 0) {

    MPI_Win_unlock_all(mhdHelBp_0_win);
    MPI_Win_unlock_all(mhdHelBt_0_win);
    MPI_Win_unlock_all(mhdHelBr_0_win);
    MPI_Win_unlock_all(mhdHelVp_0_win);
    MPI_Win_unlock_all(mhdHelVt_0_win);
    MPI_Win_unlock_all(mhdHelVr_0_win);
    MPI_Win_unlock_all(mhdHelD_0_win);
    MPI_Win_unlock_all(mhdHelBp_1_win);
    MPI_Win_unlock_all(mhdHelBt_1_win);
    MPI_Win_unlock_all(mhdHelBr_1_win);
    MPI_Win_unlock_all(mhdHelVp_1_win);
    MPI_Win_unlock_all(mhdHelVt_1_win);
    MPI_Win_unlock_all(mhdHelVr_1_win);
    MPI_Win_unlock_all(mhdHelD_1_win);

    MPI_Win_free(&mhdHelBp_0_win);
    MPI_Win_free(&mhdHelBt_0_win);
    MPI_Win_free(&mhdHelBr_0_win);
    MPI_Win_free(&mhdHelVp_0_win);
    MPI_Win_free(&mhdHelVt_0_win);
    MPI_Win_free(&mhdHelVr_0_win);
    MPI_Win_free(&mhdHelD_0_win);
    MPI_Win_free(&mhdHelBp_1_win);
    MPI_Win_free(&mhdHelBt_1_win);
    MPI_Win_free(&mhdHelBr_1_win);
    MPI_Win_free(&mhdHelVp_1_win);
    MPI_Win_free(&mhdHelVt_1_win);
    MPI_Win_free(&mhdHelVr_1_win);
    MPI_Win_free(&mhdHelD_1_win);

  }

}/*-------- END cleanupMPIWindows()  -------------------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     mhdGetInterpData ( Scalar_t dt )                     /*--*/
/*--                                                              --*/
/*-- This function checks if the simulation time is right to read --*/
/*-- or copy MHD data, and does so if needed.                     --*/
/*-- It also sets the interpolation factors 's_cor' and 's_hel'   --*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
{

  int i,N;
  int need_sync_0,need_sync_1;
  Scalar_t time_interp;

//
// Set the time we are interpolating to (the current time plus the time step).
//

  time_interp = t_global + dt;

//
// Check if MHD is being used yet.
//
  if ((t_global < mhdTime[0]) && (config.mhdSteadyState == 0)){
    mhdGridStatus = MHD_DEFAULT;
    return;
  }
  else{
    mhdGridStatus = MHD_COUPLED;
  }

//
// Find the two MHD files that bound the current time_interp, and set
// interpolation factors.
//

  if (time_interp <= mhdTime[0]){

    mhdFileIndex0 = 0;
    mhdFileIndex1 = 0;

    s_cor = 0.0;

  } else if (time_interp >= mhdTime[config.mhdNumFiles-1]){

    mhdFileIndex0 = config.mhdNumFiles-1;
    mhdFileIndex1 = config.mhdNumFiles-1;

    s_cor = 0.0;

  } else {

    for (i=1; i<config.mhdNumFiles; i++){
      if (mhdTime[i] > time_interp ){
         mhdFileIndex0 = i-1;
         mhdFileIndex1 = i;
         break;
      }
    }

    s_cor = ( time_interp - mhdTime[mhdFileIndex0] ) /
            ( mhdTime[mhdFileIndex1] - mhdTime[mhdFileIndex0]);

  }

//
// Read or copy MHD data as needed.
//

  if (mhdFileIndex0 != mhdFileIndex_loaded0){ // if state0 needs to be updated
    MPI_Barrier(comm_shared);
    if (mhdFileIndex0 == mhdFileIndex_loaded1){ // if state0 used is already stored previously in state1
      if(mpi_rank_shared==0) {
        N=(int)mhdBprDimMax[0]*(int)mhdBptDimMax[0]*(int)mhdBppDimMax[0];
        memcpy(&mhdBp_0[0],&mhdBp_1[0],N*sizeof(float));
        N=(int)mhdBtrDimMax[0]*(int)mhdBttDimMax[0]*(int)mhdBtpDimMax[0];
        memcpy(&mhdBt_0[0],&mhdBt_1[0],N*sizeof(float));
        N=(int)mhdBrrDimMax[0]*(int)mhdBrtDimMax[0]*(int)mhdBrpDimMax[0];
        memcpy(&mhdBr_0[0],&mhdBr_1[0],N*sizeof(float));
        N=(int)mhdVprDimMax[0]*(int)mhdVptDimMax[0]*(int)mhdVppDimMax[0];
        memcpy(&mhdVp_0[0],&mhdVp_1[0],N*sizeof(float));
        N=(int)mhdVtrDimMax[0]*(int)mhdVttDimMax[0]*(int)mhdVtpDimMax[0];
        memcpy(&mhdVt_0[0],&mhdVt_1[0],N*sizeof(float));
        N=(int)mhdVrrDimMax[0]*(int)mhdVrtDimMax[0]*(int)mhdVrpDimMax[0];
        memcpy(&mhdVr_0[0],&mhdVr_1[0],N*sizeof(float));
        N=(int)mhdDrDimMax[0]*(int)mhdDtDimMax[0]*(int)mhdDpDimMax[0];
        memcpy(&mhdD_0[0],&mhdD_1[0],N*sizeof(float));
      }
    } else {
      if(mpi_rank_shared==0) {
        mhdReadData(mhdFileIndex0,
                &mhdBp_0, &mhdBt_0, &mhdBr_0,
                &mhdVp_0, &mhdVt_0, &mhdVr_0,
                &mhdD_0);
      }
    }
    mhdFileIndex_loaded0 = mhdFileIndex0;
    need_sync_0 = 1;
  } else {
    need_sync_0 = 0;
  }

  if (mhdFileIndex1 != mhdFileIndex_loaded1){
    MPI_Barrier(comm_shared);
    if(mpi_rank_shared==0) {
      mhdReadData(mhdFileIndex1,
                &mhdBp_1, &mhdBt_1, &mhdBr_1,
                &mhdVp_1, &mhdVt_1, &mhdVr_1,
                &mhdD_1);
    }
    mhdFileIndex_loaded1 = mhdFileIndex1;
    need_sync_1 = 1;
  } else {
    need_sync_1 = 0;
  }

  if (need_sync_0 == 1){
    MPI_Barrier(comm_shared);
    MPI_Win_sync(mhdBp_0_win);
    MPI_Win_sync(mhdBt_0_win);
    MPI_Win_sync(mhdBr_0_win);
    MPI_Win_sync(mhdVp_0_win);
    MPI_Win_sync(mhdVt_0_win);
    MPI_Win_sync(mhdVr_0_win);
    MPI_Win_sync(mhdD_0_win);
  }

  if (need_sync_1 == 1){
    MPI_Barrier(comm_shared);
    MPI_Win_sync(mhdBp_1_win);
    MPI_Win_sync(mhdBt_1_win);
    MPI_Win_sync(mhdBr_1_win);
    MPI_Win_sync(mhdVp_1_win);
    MPI_Win_sync(mhdVt_1_win);
    MPI_Win_sync(mhdVr_1_win);
    MPI_Win_sync(mhdD_1_win);
  }

  if (need_sync_0 + need_sync_1 > 0) MPI_Barrier(comm_shared);
//
// Now do everything with the heliosphere.
//
  if (config.mhdHelCouple > 0) {

//
// Find the two MHD HELIO files that bound the current time_interp, and set
// interpolation factors.
//
    if (time_interp < mhdHelTime[0]){

      mhdHelFileIndex0 = 0;
      mhdHelFileIndex1 = 0;

      s_hel = 0.0;

    } else if (time_interp > mhdHelTime[config.mhdHelNumFiles-1]){

      mhdHelFileIndex0 = config.mhdHelNumFiles-1;
      mhdHelFileIndex1 = config.mhdHelNumFiles-1;

      s_hel = 0.0;

    } else {

      for (i=1; i<config.mhdHelNumFiles; i++){
        if (mhdHelTime[i] > time_interp){
          mhdHelFileIndex0 = i-1;
          mhdHelFileIndex1 = i;
          break;
        }
      }

      s_hel = ( time_interp - mhdHelTime[mhdHelFileIndex0] ) /
              ( mhdHelTime[mhdHelFileIndex1] - mhdHelTime[mhdHelFileIndex0]);

    }

//
// Read or copy MHD HELIO data as needed.
//
    if (mhdHelFileIndex0 != mhdHelFileIndex_loaded0){// If new data needs to be read into state 0
      MPI_Barrier(comm_shared);
      if (mhdHelFileIndex0 == mhdHelFileIndex_loaded1)// If data already exists in previous state 1
      {
        if(mpi_rank_shared==0) {
          N=(int)mhdHelBprDimMax[0]*(int)mhdHelBptDimMax[0]*(int)mhdHelBppDimMax[0];
          memcpy(&mhdHelBp_0[0],&mhdHelBp_1[0],N*sizeof(float));
          N=(int)mhdHelBtrDimMax[0]*(int)mhdHelBttDimMax[0]*(int)mhdHelBtpDimMax[0];
          memcpy(&mhdHelBt_0[0],&mhdHelBt_1[0],N*sizeof(float));
          N=(int)mhdHelBrrDimMax[0]*(int)mhdHelBrtDimMax[0]*(int)mhdHelBrpDimMax[0];
          memcpy(&mhdHelBr_0[0],&mhdHelBr_1[0],N*sizeof(float));
          N=(int)mhdHelVprDimMax[0]*(int)mhdHelVptDimMax[0]*(int)mhdHelVppDimMax[0];
          memcpy(&mhdHelVp_0[0],&mhdHelVp_1[0],N*sizeof(float));
          N=(int)mhdHelVtrDimMax[0]*(int)mhdHelVttDimMax[0]*(int)mhdHelVtpDimMax[0];
          memcpy(&mhdHelVt_0[0],&mhdHelVt_1[0],N*sizeof(float));
          N=(int)mhdHelVrrDimMax[0]*(int)mhdHelVrtDimMax[0]*(int)mhdHelVrpDimMax[0];
          memcpy(&mhdHelVr_0[0],&mhdHelVr_1[0],N*sizeof(float));
          N=(int)mhdHelDrDimMax[0]*(int)mhdHelDtDimMax[0]*(int)mhdHelDpDimMax[0];
          memcpy(&mhdHelD_0[0],&mhdHelD_1[0],N*sizeof(float));
        }
      }
      else {// Read in new data
        if(mpi_rank_shared==0) {
          mhdHelReadData(mhdHelFileIndex0,
                  &mhdHelBp_0, &mhdHelBt_0, &mhdHelBr_0,
                  &mhdHelVp_0, &mhdHelVt_0, &mhdHelVr_0,
                  &mhdHelD_0);
        }
      }
      mhdHelFileIndex_loaded0 = mhdHelFileIndex0;
      need_sync_0 = 1;
    } else {// Nothing changed
      need_sync_0 = 0;
    }

    if (mhdHelFileIndex1 != mhdHelFileIndex_loaded1){ //  Check if state1 needs to be read in
      MPI_Barrier(comm_shared);
      if(mpi_rank_shared==0) {
        mhdHelReadData(mhdHelFileIndex1,
                  &mhdHelBp_1, &mhdHelBt_1, &mhdHelBr_1,
                  &mhdHelVp_1, &mhdHelVt_1, &mhdHelVr_1,
                  &mhdHelD_1);
      }
      mhdHelFileIndex_loaded1 = mhdHelFileIndex1;
      need_sync_1 = 1;
    } else {
      need_sync_1 = 0;
    }

    if (need_sync_0 == 1){
      MPI_Barrier(comm_shared);
      MPI_Win_sync(mhdHelBp_0_win);
      MPI_Win_sync(mhdHelBt_0_win);
      MPI_Win_sync(mhdHelBr_0_win);
      MPI_Win_sync(mhdHelVp_0_win);
      MPI_Win_sync(mhdHelVt_0_win);
      MPI_Win_sync(mhdHelVr_0_win);
      MPI_Win_sync(mhdHelD_0_win);
    }

    if (need_sync_1 == 1){
      MPI_Barrier(comm_shared);
      MPI_Win_sync(mhdHelBp_1_win);
      MPI_Win_sync(mhdHelBt_1_win);
      MPI_Win_sync(mhdHelBr_1_win);
      MPI_Win_sync(mhdHelVp_1_win);
      MPI_Win_sync(mhdHelVt_1_win);
      MPI_Win_sync(mhdHelVr_1_win);
      MPI_Win_sync(mhdHelD_1_win);
    }

    if (need_sync_0 + need_sync_1 > 0) MPI_Barrier(comm_shared);

  }
//
// Tell everyone what we just did.
//
  if (mpi_rank==0) {
    printf("  --> MHD(C) Couple: Time: %14.8e  Time+Dt: %14.8e  s_cor: %10.8f\n",
          t_global,
          time_interp,s_cor);
    printf("  -->                idx0: %03d  idx1: %03d  mhdTime0: %14.8e  mhdTime1: %14.8e\n",
          mhdFileIndex0+1,
          mhdFileIndex1+1,
          mhdTime[mhdFileIndex0],
          mhdTime[mhdFileIndex1]);
    if (config.mhdHelCouple > 0){
      printf("  --> MHD(H) Couple: Time: %14.8e  Time+Dt: %14.8e  s_hel: %10.8f\n",
          t_global,
          time_interp,s_hel);
      printf("  -->                idx0: %03d  idx1: %03d  mhdTime0: %14.8e  mhdTime1: %14.8e\n",
          mhdHelFileIndex0+1,
          mhdHelFileIndex1+1,
          mhdHelTime[mhdHelFileIndex0],
          mhdHelTime[mhdHelFileIndex1]);
    }
  }

}
/*----------------- END mhdGetInterpData(dt)  ------------------------*/
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--*/     void                                                   /*--*/
/*--*/     mhdReadData(Index_t fileIndex,                         /*--*/
/*--*/            float *mhdBp[], float *mhdBt[], float *mhdBr[], /*--*/
/*--*/            float *mhdVp[], float *mhdVt[], float *mhdVr[], /*--*/
/*--*/            float *mhdD[])                                  /*--*/
/*--*/                                                            /*--*/
/*--                                                                --*/
/*--This function reads the MHD data from a HDF file.              --*/
/*--------------------------------------------------------------------*/
{/*-------------------------------------------------------------------*/
  char fileNames[7][MAX_STRING_SIZE];

  double timer_tmp = 0;

  timer_tmp = MPI_Wtime();

  if (config.mhdDigits == 3) {

    sprintf(fileNames[0], "%sbp%03d%s", config.mhdDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[1], "%sbt%03d%s", config.mhdDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[2], "%sbr%03d%s", config.mhdDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[3], "%svp%03d%s", config.mhdDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[4], "%svt%03d%s", config.mhdDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[5], "%svr%03d%s", config.mhdDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[6], "%srho%03d%s", config.mhdDirectory, fileIndex + 1, file_extension);

  } else {

    sprintf(fileNames[0], "%sbp%06d%s", config.mhdDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[1], "%sbt%06d%s", config.mhdDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[2], "%sbr%06d%s", config.mhdDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[3], "%svp%06d%s", config.mhdDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[4], "%svt%06d%s", config.mhdDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[5], "%svr%06d%s", config.mhdDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[6], "%srho%06d%s", config.mhdDirectory, fileIndex + 1, file_extension);

  }

  // reading in Bp
  mhdReadDatafromFile(fileNames[0], mhdBp);

  // reading in Bt
  mhdReadDatafromFile(fileNames[1], mhdBt);

  // reading in Br
  mhdReadDatafromFile(fileNames[2], mhdBr);

  // reading in Vp
  mhdReadDatafromFile(fileNames[3], mhdVp);

  // reading in Vt
  mhdReadDatafromFile(fileNames[4], mhdVt);

  // reading in Vr
  mhdReadDatafromFile(fileNames[5], mhdVr);

  // reading in D
  mhdReadDatafromFile(fileNames[6], mhdD);

  if (mpi_rank == 0) printf("  --> IO MHD: Coronal sequence %03d read.\n",fileIndex+1);

  timer_mhd_io = timer_mhd_io + (MPI_Wtime() - timer_tmp);

}/*-------- END mhdReadData()  -----------------------*/
/*------------------------------------------------------------------*/


/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--*/     void                                                   /*--*/
/*--*/     mhdHelReadData(Index_t fileIndex,                      /*--*/
/*--*/            float *mhdBp[], float *mhdBt[], float *mhdBr[], /*--*/
/*--*/            float *mhdVp[], float *mhdVt[], float *mhdVr[], /*--*/
/*--*/            float *mhdD[])                                  /*--*/
/*--*/                                                            /*--*/
/*--                                                                --*/
/*--This function checks if the simulation time is right to load    --*/
/*--another timestep file from the MHD HEL data, and does if so.    --*/
/*--------------------------------------------------------------------*/
{/*-------------------------------------------------------------------*/

  char fileNames[7][MAX_STRING_SIZE];

  double timer_tmp = 0;

  timer_tmp = MPI_Wtime();

  if (config.mhdHelDigits == 3) {

    sprintf(fileNames[0], "%sbp%03d%s", config.mhdHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[1], "%sbt%03d%s", config.mhdHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[2], "%sbr%03d%s", config.mhdHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[3], "%svp%03d%s", config.mhdHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[4], "%svt%03d%s", config.mhdHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[5], "%svr%03d%s", config.mhdHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[6], "%srho%03d%s", config.mhdHelDirectory, fileIndex + 1, file_extension);
  }
  else
  {
    sprintf(fileNames[0], "%sbp%06d%s", config.mhdHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[1], "%sbt%06d%s", config.mhdHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[2], "%sbr%06d%s", config.mhdHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[3], "%svp%06d%s", config.mhdHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[4], "%svt%06d%s", config.mhdHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[5], "%svr%06d%s", config.mhdHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[6], "%srho%06d%s", config.mhdHelDirectory, fileIndex + 1, file_extension);
  }

  // reading in Bp
  mhdReadDatafromFile(fileNames[0], mhdBp);

  // reading in Bt
  mhdReadDatafromFile(fileNames[1], mhdBt);

  // reading in Br
  mhdReadDatafromFile(fileNames[2], mhdBr);

  // reading in Vp
  mhdReadDatafromFile(fileNames[3], mhdVp);

  // reading in Vt
  mhdReadDatafromFile(fileNames[4], mhdVt);

  // reading in Vr
  mhdReadDatafromFile(fileNames[5], mhdVr);

  // reading in D
  mhdReadDatafromFile(fileNames[6], mhdD);

  if (mpi_rank == 0) printf("  --> IO MHD: Helio sequence %03d read.\n",fileIndex+1);

  timer_mhd_io = timer_mhd_io + (MPI_Wtime() - timer_tmp);

}/*-------- END mhdHelReadData()  ----------------------------------*/
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     mhdReadFieldIndex(void)                              /*--*/
/*--                                                              --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  char fileNames[7][MAX_STRING_SIZE];

  Scalar_t rMin, rMax, rTemp;

  double timer_tmp=0;

  timer_tmp = MPI_Wtime();

  if (config.mhdDigits == 3) {

    sprintf(fileNames[0], "%sbp%03d%s", config.mhdDirectory, 1, file_extension);
    sprintf(fileNames[1], "%sbt%03d%s", config.mhdDirectory, 1, file_extension);
    sprintf(fileNames[2], "%sbr%03d%s", config.mhdDirectory, 1, file_extension);
    sprintf(fileNames[3], "%svp%03d%s", config.mhdDirectory, 1, file_extension);
    sprintf(fileNames[4], "%svt%03d%s", config.mhdDirectory, 1, file_extension);
    sprintf(fileNames[5], "%svr%03d%s", config.mhdDirectory, 1, file_extension);
    sprintf(fileNames[6], "%srho%03d%s", config.mhdDirectory, 1, file_extension);

  } else {

    sprintf(fileNames[0], "%sbp%06d%s", config.mhdDirectory, 1, file_extension);
    sprintf(fileNames[1], "%sbt%06d%s", config.mhdDirectory, 1, file_extension);
    sprintf(fileNames[2], "%sbr%06d%s", config.mhdDirectory, 1, file_extension);
    sprintf(fileNames[3], "%svp%06d%s", config.mhdDirectory, 1, file_extension);
    sprintf(fileNames[4], "%svt%06d%s", config.mhdDirectory, 1, file_extension);
    sprintf(fileNames[5], "%svr%06d%s", config.mhdDirectory, 1, file_extension);
    sprintf(fileNames[6], "%srho%06d%s", config.mhdDirectory, 1, file_extension);

  }

  // Read in dimensions --> Allocate the mesh storage --> Read the mesh

    // mhdBpp
    mhdReadMeshDimensions(fileNames[0], "dim3", 0, &mhdBppDimMax[0]);
    mhdBppDim = (float *)malloc(sizeof(float) * (int)(mhdBppDimMax[0]));
    mhdReadMesh(fileNames[0], "dim3", 0, &mhdBppDim);

    // mhdBptDim
    mhdReadMeshDimensions(fileNames[0], "dim2", 1, &mhdBptDimMax[0]);
    mhdBptDim = (float *)malloc(sizeof(float) * (int)(mhdBptDimMax[0]));
    mhdReadMesh(fileNames[0], "dim2", 1, &mhdBptDim);

    // mhdBprDim
    mhdReadMeshDimensions(fileNames[0], "dim1", 2, &mhdBprDimMax[0]);
    mhdBprDim = (float *)malloc(sizeof(float) * (int)(mhdBprDimMax[0]));
    mhdReadMesh(fileNames[0], "dim1", 2, &mhdBprDim);

    // grab the min and max for the r
    rMin = mhdBprDim[0];
    rMax = mhdBprDim[mhdBprDimMax[0] - 1];

    // mhdBtpDim
    mhdReadMeshDimensions(fileNames[1], "dim3", 0, &mhdBtpDimMax[0]);
    mhdBtpDim = (float *)malloc(sizeof(float) * (int)(mhdBtpDimMax[0]));
    mhdReadMesh(fileNames[1], "dim3", 0, &mhdBtpDim);

    // mhdBttDim
    mhdReadMeshDimensions(fileNames[1], "dim2", 1, &mhdBttDimMax[0]);
    mhdBttDim = (float *)malloc(sizeof(float) * (int)(mhdBttDimMax[0]));
    mhdReadMesh(fileNames[1], "dim2", 1, &mhdBttDim);

    // mhdBtrDim
    mhdReadMeshDimensions(fileNames[1], "dim1", 2, &mhdBtrDimMax[0]);
    mhdBtrDim = (float *)malloc(sizeof(float) * (int)(mhdBtrDimMax[0]));
    mhdReadMesh(fileNames[1], "dim1", 2, &mhdBtrDim);

    // grab the min and max for the r
    rTemp = mhdBtrDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = mhdBtrDim[mhdBtrDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // mhdBrpDim
    mhdReadMeshDimensions(fileNames[2], "dim3", 0, &mhdBrpDimMax[0]);
    mhdBrpDim = (float *)malloc(sizeof(float) * (int)(mhdBrpDimMax[0]));
    mhdReadMesh(fileNames[2], "dim3", 0, &mhdBrpDim);

    // mhdBrtDim
    mhdReadMeshDimensions(fileNames[2], "dim2", 1, &mhdBrtDimMax[0]);
    mhdBrtDim = (float *)malloc(sizeof(float) * (int)(mhdBrtDimMax[0]));
    mhdReadMesh(fileNames[2], "dim2", 1, &mhdBrtDim);

    // mhdBrrDim
    mhdReadMeshDimensions(fileNames[2], "dim1", 2, &mhdBrrDimMax[0]);
    mhdBrrDim = (float *)malloc(sizeof(float) * (int)(mhdBrrDimMax[0]));
    mhdReadMesh(fileNames[2], "dim1", 2, &mhdBrrDim);

    // grab the min and max for the r
    rTemp = mhdBrrDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = mhdBrrDim[mhdBrrDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // mhdVppDim
    mhdReadMeshDimensions(fileNames[3], "dim3", 0, &mhdVppDimMax[0]);
    mhdVppDim = (float *)malloc(sizeof(float) * (int)(mhdVppDimMax[0]));
    mhdReadMesh(fileNames[3], "dim3", 0, &mhdVppDim);

    // mhdVptDim
    mhdReadMeshDimensions(fileNames[3], "dim2", 1, &mhdVptDimMax[0]);
    mhdVptDim = (float *)malloc(sizeof(float) * (int)(mhdVptDimMax[0]));
    mhdReadMesh(fileNames[3], "dim2", 1, &mhdVptDim);

    // mhdVprDim
    mhdReadMeshDimensions(fileNames[3], "dim1", 2, &mhdVprDimMax[0]);
    mhdVprDim = (float *)malloc(sizeof(float) * (int)(mhdVprDimMax[0]));
    mhdReadMesh(fileNames[3], "dim1", 2, &mhdVprDim);

    // grab the min and max for the r
    rTemp = mhdVprDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = mhdVprDim[mhdVprDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // mhdVtpDim
    mhdReadMeshDimensions(fileNames[4], "dim3", 0, &mhdVtpDimMax[0]);
    mhdVtpDim = (float *)malloc(sizeof(float) * (int)(mhdVtpDimMax[0]));
    mhdReadMesh(fileNames[4], "dim3", 0, &mhdVtpDim);

    // mhdVttDim
    mhdReadMeshDimensions(fileNames[4], "dim2", 1, &mhdVttDimMax[0]);
    mhdVttDim = (float *)malloc(sizeof(float) * (int)(mhdVttDimMax[0]));
    mhdReadMesh(fileNames[4], "dim2", 1, &mhdVttDim);

    // mhdVtrDim
    mhdReadMeshDimensions(fileNames[4], "dim1", 2, &mhdVtrDimMax[0]);
    mhdVtrDim = (float *)malloc(sizeof(float) * (int)(mhdVtrDimMax[0]));
    mhdReadMesh(fileNames[4], "dim1", 2, &mhdVtrDim);

    // grab the min and max for the r
    rTemp = mhdVtrDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = mhdVtrDim[mhdVtrDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // mhdVrpDim
    mhdReadMeshDimensions(fileNames[5], "dim3", 0, &mhdVrpDimMax[0]);
    mhdVrpDim = (float *)malloc(sizeof(float) * (int)(mhdVrpDimMax[0]));
    mhdReadMesh(fileNames[5], "dim3", 0, &mhdVrpDim);

    // mhdVrtDim
    mhdReadMeshDimensions(fileNames[5], "dim2", 1, &mhdVrtDimMax[0]);
    mhdVrtDim = (float *)malloc(sizeof(float) * (int)(mhdVrtDimMax[0]));
    mhdReadMesh(fileNames[5], "dim2", 1, &mhdVrtDim);

    // mhdVrrDim
    mhdReadMeshDimensions(fileNames[5], "dim1", 2, &mhdVrrDimMax[0]);
    mhdVrrDim = (float *)malloc(sizeof(float) * (int)(mhdVrrDimMax[0]));
    mhdReadMesh(fileNames[5], "dim1", 2, &mhdVrrDim);

    // grab the min and max for the r
    rTemp = mhdVrrDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = mhdVrrDim[mhdVrrDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // mhdDpDim
    mhdReadMeshDimensions(fileNames[6], "dim3", 0, &mhdDpDimMax[0]);
    mhdDpDim = (float *)malloc(sizeof(float) * (int)(mhdDpDimMax[0]));
    mhdReadMesh(fileNames[6], "dim3", 0, &mhdDpDim);

    // mhdDtDim
    mhdReadMeshDimensions(fileNames[6], "dim2", 1, &mhdDtDimMax[0]);
    mhdDtDim = (float *)malloc(sizeof(float) * (int)(mhdDtDimMax[0]));
    mhdReadMesh(fileNames[6], "dim2", 1, &mhdDtDim);

    // mhdDrDim
    mhdReadMeshDimensions(fileNames[6], "dim1", 2, &mhdDrDimMax[0]);
    mhdDrDim = (float *)malloc(sizeof(float) * (int)(mhdDrDimMax[0]));
    mhdReadMesh(fileNames[6], "dim1", 2, &mhdDrDim);

    // grab the min and max for the r
    rTemp = mhdDrDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = mhdDrDim[mhdDrDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

  // set rScale, and mhdRadialMin/Max
  config.rScale = rMin * RSAU;
  config.mhdRadialMin = rMin;
  config.mhdRadialMax = rMax;

  timer_mhd_io = timer_mhd_io + (MPI_Wtime() - timer_tmp);

}/*-------- END mhdReadFieldIndex()  -----------------------*/
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     mhdHelReadFieldIndex(void)                           /*--*/
/*--                                                              --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  char fileNames[7][MAX_STRING_SIZE];
  Scalar_t rMin, rMax, rTemp;

  double timer_tmp = 0;

  timer_tmp = MPI_Wtime();

  if (config.mhdHelDigits == 3)
  {
    sprintf(fileNames[0], "%sbp%03d%s", config.mhdHelDirectory, 1, file_extension);
    sprintf(fileNames[1], "%sbt%03d%s", config.mhdHelDirectory, 1, file_extension);
    sprintf(fileNames[2], "%sbr%03d%s", config.mhdHelDirectory, 1, file_extension);
    sprintf(fileNames[3], "%svp%03d%s", config.mhdHelDirectory, 1, file_extension);
    sprintf(fileNames[4], "%svt%03d%s", config.mhdHelDirectory, 1, file_extension);
    sprintf(fileNames[5], "%svr%03d%s", config.mhdHelDirectory, 1, file_extension);
    sprintf(fileNames[6], "%srho%03d%s", config.mhdHelDirectory,1, file_extension);
  }
  else
  {
    sprintf(fileNames[0], "%sbp%06d%s", config.mhdHelDirectory,  1, file_extension);
    sprintf(fileNames[1], "%sbt%06d%s", config.mhdHelDirectory,  1, file_extension);
    sprintf(fileNames[2], "%sbr%06d%s", config.mhdHelDirectory,  1, file_extension);
    sprintf(fileNames[3], "%svp%06d%s", config.mhdHelDirectory,  1, file_extension);
    sprintf(fileNames[4], "%svt%06d%s", config.mhdHelDirectory,  1, file_extension);
    sprintf(fileNames[5], "%svr%06d%s", config.mhdHelDirectory,  1, file_extension);
    sprintf(fileNames[6], "%srho%06d%s", config.mhdHelDirectory, 1, file_extension);
  }

    // mhdHelBpp
    mhdReadMeshDimensions(fileNames[0], "dim3", 0, &mhdHelBppDimMax[0]);
    mhdHelBppDim = (float *)malloc(sizeof(float) * (int)(mhdHelBppDimMax[0]));
    mhdReadMesh(fileNames[0], "dim3", 0, &mhdHelBppDim);

    // mhdHelBptDim
    mhdReadMeshDimensions(fileNames[0], "dim2", 1, &mhdHelBptDimMax[0]);
    mhdHelBptDim = (float *)malloc(sizeof(float) * (int)(mhdHelBptDimMax[0]));
    mhdReadMesh(fileNames[0], "dim2", 1, &mhdHelBptDim);

    // mhdHelBprDim
    mhdReadMeshDimensions(fileNames[0], "dim1", 2, &mhdHelBprDimMax[0]);
    mhdHelBprDim = (float *)malloc(sizeof(float) * (int)(mhdHelBprDimMax[0]));
    mhdReadMesh(fileNames[0], "dim1", 2, &mhdHelBprDim);

    // grab the min and max for the r
    rMin = mhdHelBprDim[0];
    rMax = mhdHelBprDim[mhdHelBprDimMax[0] - 1];

    // mhdHelBtpDim
    mhdReadMeshDimensions(fileNames[1], "dim3", 0, &mhdHelBtpDimMax[0]);
    mhdHelBtpDim = (float *)malloc(sizeof(float) * (int)(mhdHelBtpDimMax[0]));
    mhdReadMesh(fileNames[1], "dim3", 0, &mhdHelBtpDim);

    // mhdHelBttDim
    mhdReadMeshDimensions(fileNames[1], "dim2", 1, &mhdHelBttDimMax[0]);
    mhdHelBttDim = (float *)malloc(sizeof(float) * (int)(mhdHelBttDimMax[0]));
    mhdReadMesh(fileNames[1], "dim2", 1, &mhdHelBttDim);

    // mhdHelBtrDim
    mhdReadMeshDimensions(fileNames[1], "dim1", 2, &mhdHelBtrDimMax[0]);
    mhdHelBtrDim = (float *)malloc(sizeof(float) * (int)(mhdHelBtrDimMax[0]));
    mhdReadMesh(fileNames[1], "dim1", 2, &mhdHelBtrDim);

    // grab the min and max for the r
    rTemp = mhdHelBtrDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = mhdHelBtrDim[mhdHelBtrDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // mhdHelBrpDim
    mhdReadMeshDimensions(fileNames[2], "dim3", 0, &mhdHelBrpDimMax[0]);
    mhdHelBrpDim = (float *)malloc(sizeof(float) * (int)(mhdHelBrpDimMax[0]));
    mhdReadMesh(fileNames[2], "dim3", 0, &mhdHelBrpDim);

    // mhdHelBrtDim
    mhdReadMeshDimensions(fileNames[2], "dim2", 1, &mhdHelBrtDimMax[0]);
    mhdHelBrtDim = (float *)malloc(sizeof(float) * (int)(mhdHelBrtDimMax[0]));
    mhdReadMesh(fileNames[2], "dim2", 1, &mhdHelBrtDim);

    // mhdHelBrrDim
    mhdReadMeshDimensions(fileNames[2], "dim1", 2, &mhdHelBrrDimMax[0]);
    mhdHelBrrDim = (float *)malloc(sizeof(float) * (int)(mhdHelBrrDimMax[0]));
    mhdReadMesh(fileNames[2], "dim1", 2, &mhdHelBrrDim);

    // grab the min and max for the r
    rTemp = mhdHelBrrDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = mhdHelBrrDim[mhdHelBrrDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // mhdHelVppDim
    mhdReadMeshDimensions(fileNames[3], "dim3", 0, &mhdHelVppDimMax[0]);
    mhdHelVppDim = (float *)malloc(sizeof(float) * (int)(mhdHelVppDimMax[0]));
    mhdReadMesh(fileNames[3], "dim3", 0, &mhdHelVppDim);

    // mhdHelVptDim
    mhdReadMeshDimensions(fileNames[3], "dim2", 1, &mhdHelVptDimMax[0]);
    mhdHelVptDim = (float *)malloc(sizeof(float) * (int)(mhdHelVptDimMax[0]));
    mhdReadMesh(fileNames[3], "dim2", 1, &mhdHelVptDim);

    // mhdHelVprDim
    mhdReadMeshDimensions(fileNames[3], "dim1", 2, &mhdHelVprDimMax[0]);
    mhdHelVprDim = (float *)malloc(sizeof(float) * (int)(mhdHelVprDimMax[0]));
    mhdReadMesh(fileNames[3], "dim1", 2, &mhdHelVprDim);

    // grab the min and max for the r
    rTemp = mhdHelVprDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = mhdHelVprDim[mhdHelVprDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // mhdHelVtpDim
    mhdReadMeshDimensions(fileNames[4], "dim3", 0, &mhdHelVtpDimMax[0]);
    mhdHelVtpDim = (float *)malloc(sizeof(float) * (int)(mhdHelVtpDimMax[0]));
    mhdReadMesh(fileNames[4], "dim3", 0, &mhdHelVtpDim);

    // mhdHelVttDim
    mhdReadMeshDimensions(fileNames[4], "dim2", 1, &mhdHelVttDimMax[0]);
    mhdHelVttDim = (float *)malloc(sizeof(float) * (int)(mhdHelVttDimMax[0]));
    mhdReadMesh(fileNames[4], "dim2", 1, &mhdHelVttDim);

    // mhdHelVtrDim
    mhdReadMeshDimensions(fileNames[4], "dim1", 2, &mhdHelVtrDimMax[0]);
    mhdHelVtrDim = (float *)malloc(sizeof(float) * (int)(mhdHelVtrDimMax[0]));
    mhdReadMesh(fileNames[4], "dim1", 2, &mhdHelVtrDim);

    // grab the min and max for the r
    rTemp = mhdHelVtrDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = mhdHelVtrDim[mhdHelVtrDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // mhdHelVrpDim
    mhdReadMeshDimensions(fileNames[5], "dim3", 0, &mhdHelVrpDimMax[0]);
    mhdHelVrpDim = (float *)malloc(sizeof(float) * (int)(mhdHelVrpDimMax[0]));
    mhdReadMesh(fileNames[5], "dim3", 0, &mhdHelVrpDim);

    // mhdHelVrtDim
    mhdReadMeshDimensions(fileNames[5], "dim2", 1, &mhdHelVrtDimMax[0]);
    mhdHelVrtDim = (float *)malloc(sizeof(float) * (int)(mhdHelVrtDimMax[0]));
    mhdReadMesh(fileNames[5], "dim2", 1, &mhdHelVrtDim);

    // mhdHelVrrDim
    mhdReadMeshDimensions(fileNames[5], "dim1", 2, &mhdHelVrrDimMax[0]);
    mhdHelVrrDim = (float *)malloc(sizeof(float) * (int)(mhdHelVrrDimMax[0]));
    mhdReadMesh(fileNames[5], "dim1", 2, &mhdHelVrrDim);

    // grab the min and max for the r
    rTemp = mhdHelVrrDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = mhdHelVrrDim[mhdHelVrrDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // mhdHelDpDim
    mhdReadMeshDimensions(fileNames[6], "dim3", 0, &mhdHelDpDimMax[0]);
    mhdHelDpDim = (float *)malloc(sizeof(float) * (int)(mhdHelDpDimMax[0]));
    mhdReadMesh(fileNames[6], "dim3", 0, &mhdHelDpDim);

    // mhdHelDtDim
    mhdReadMeshDimensions(fileNames[6], "dim2", 1, &mhdHelDtDimMax[0]);
    mhdHelDtDim = (float *)malloc(sizeof(float) * (int)(mhdHelDtDimMax[0]));
    mhdReadMesh(fileNames[6], "dim2", 1, &mhdHelDtDim);

    // mhdHelDrDim
    mhdReadMeshDimensions(fileNames[6], "dim1", 2, &mhdHelDrDimMax[0]);
    mhdHelDrDim = (float *)malloc(sizeof(float) * (int)(mhdHelDrDimMax[0]));
    mhdReadMesh(fileNames[6], "dim1", 2, &mhdHelDrDim);

    // grab the min and max for the r
    rTemp = mhdHelDrDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = mhdHelDrDim[mhdHelDrDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

  // set rScale, and mhdHelRadialMin/Max
  config.mhdHelRadialMin = rMin;
  config.mhdHelRadialMax = rMax;

  timer_mhd_io = timer_mhd_io + (MPI_Wtime() - timer_tmp);

}/*-------- END mhdHelReadFieldIndex()  ----------------------------*/
/*------------------------------------------------------------------*/
