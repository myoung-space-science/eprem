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

#include <mfhdf.h>
#include "readMAS.h"
#include "mpiInit.h"
#include "global.h"
#include "configuration.h"
#include "error.h"
#include "simCore.h"
#include "flow.h"
#include "observerOutput.h"
#include "timers.h"

Scalar_t *masTime;
Scalar_t *masHelTime;

Scalar_t phiOffset;
Scalar_t phiHelOffset;

Index_t masFileIndex0;
Index_t masFileIndex1;
Index_t masHelFileIndex0;
Index_t masHelFileIndex1;

Index_t masFileIndex_loaded0=-9999;
Index_t masFileIndex_loaded1=-9999;
Index_t masHelFileIndex_loaded0=-9999;
Index_t masHelFileIndex_loaded1=-9999;

Index_t masMallocFlag;
Index_t masEqFileFlag;

int32 masDimMin[1] = {0};

int32 masBppDimMax[1];
int32 masBptDimMax[1];
int32 masBprDimMax[1];
int32 masHelBppDimMax[1];
int32 masHelBptDimMax[1];
int32 masHelBprDimMax[1];

int32 masBtpDimMax[1];
int32 masBttDimMax[1];
int32 masBtrDimMax[1];
int32 masHelBtpDimMax[1];
int32 masHelBttDimMax[1];
int32 masHelBtrDimMax[1];

int32 masBrpDimMax[1];
int32 masBrtDimMax[1];
int32 masBrrDimMax[1];
int32 masHelBrpDimMax[1];
int32 masHelBrtDimMax[1];
int32 masHelBrrDimMax[1];

int32 masVppDimMax[1];
int32 masVptDimMax[1];
int32 masVprDimMax[1];
int32 masHelVppDimMax[1];
int32 masHelVptDimMax[1];
int32 masHelVprDimMax[1];

int32 masVtpDimMax[1];
int32 masVttDimMax[1];
int32 masVtrDimMax[1];
int32 masHelVtpDimMax[1];
int32 masHelVttDimMax[1];
int32 masHelVtrDimMax[1];

int32 masVrpDimMax[1];
int32 masVrtDimMax[1];
int32 masVrrDimMax[1];
int32 masHelVrpDimMax[1];
int32 masHelVrtDimMax[1];
int32 masHelVrrDimMax[1];

int32 masDpDimMax[1];
int32 masDtDimMax[1];
int32 masDrDimMax[1];
int32 masHelDpDimMax[1];
int32 masHelDtDimMax[1];
int32 masHelDrDimMax[1];

float * masBppDim;
float * masBptDim;
float * masBprDim;
float * masHelBppDim;
float * masHelBptDim;
float * masHelBprDim;

float * masBtpDim;
float * masBttDim;
float * masBtrDim;
float * masHelBtpDim;
float * masHelBttDim;
float * masHelBtrDim;

float * masBrpDim;
float * masBrtDim;
float * masBrrDim;
float * masHelBrpDim;
float * masHelBrtDim;
float * masHelBrrDim;

float * masVppDim;
float * masVptDim;
float * masVprDim;
float * masHelVppDim;
float * masHelVptDim;
float * masHelVprDim;

float * masVtpDim;
float * masVttDim;
float * masVtrDim;
float * masHelVtpDim;
float * masHelVttDim;
float * masHelVtrDim;

float * masVrpDim;
float * masVrtDim;
float * masVrrDim;
float * masHelVrpDim;
float * masHelVrtDim;
float * masHelVrrDim;

float * masDpDim;
float * masDtDim;
float * masDrDim;
float * masHelDpDim;
float * masHelDtDim;
float * masHelDrDim;


float * masBp_0;
float * masBt_0;
float * masBr_0;
float * masVp_0;
float * masVt_0;
float * masVr_0;
float * masD_0;
float * masHelBp_0;
float * masHelBt_0;
float * masHelBr_0;
float * masHelVp_0;
float * masHelVt_0;
float * masHelVr_0;
float * masHelD_0;

float * masBp_1;
float * masBt_1;
float * masBr_1;
float * masVp_1;
float * masVt_1;
float * masVr_1;
float * masD_1;
float * masHelBp_1;
float * masHelBt_1;
float * masHelBr_1;
float * masHelVp_1;
float * masHelVt_1;
float * masHelVr_1;
float * masHelD_1;

MPI_Win masBp_0_win;
MPI_Win masBt_0_win;
MPI_Win masBr_0_win;
MPI_Win masVp_0_win;
MPI_Win masVt_0_win;
MPI_Win masVr_0_win;
MPI_Win masD_0_win;
MPI_Win masHelBp_0_win;
MPI_Win masHelBt_0_win;
MPI_Win masHelBr_0_win;
MPI_Win masHelVp_0_win;
MPI_Win masHelVt_0_win;
MPI_Win masHelVr_0_win;
MPI_Win masHelD_0_win;

MPI_Win masBp_1_win;
MPI_Win masBt_1_win;
MPI_Win masBr_1_win;
MPI_Win masVp_1_win;
MPI_Win masVt_1_win;
MPI_Win masVr_1_win;
MPI_Win masD_1_win;
MPI_Win masHelBp_1_win;
MPI_Win masHelBt_1_win;
MPI_Win masHelBr_1_win;
MPI_Win masHelVp_1_win;
MPI_Win masHelVt_1_win;
MPI_Win masHelVr_1_win;
MPI_Win masHelD_1_win;

int coupleStarted=0;


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
/*--*/     masFetchCouplingInfo(void)                           /*--*/
/*--                                                              --*/
/*--  This function loads coupling info from the MAS directories  --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  FILE *rfile;

  int max = 10000;
  char line[10000];

  int nFileLines;

  char *name = NULL;
  char *value = NULL;
  char delims[] = ": ";

  char masHelInfoFilename[MAX_STRING_SIZE] = "mas_helio_run_info.txt";
  char masHelInfoFilenameWithPath[MAX_STRING_SIZE];

  // -- coronal coupling --//

  // -- heliospheric coupling --//
  if (config.masHelCouple > 0) {

    // build the path to the info file
    sprintf(masHelInfoFilenameWithPath, "%s%s", config.masHelDirectory, masHelInfoFilename);

    // attempt to open the info file
    rfile = fopen(masHelInfoFilenameWithPath, "r");
    if (rfile==NULL) {
      printf("ERROR - Could not open file \"%s\"\nReverting to defaults\n", masHelInfoFilenameWithPath);
    } else {
      if (mpi_rank == 0) {
        printf("Reading parameters from \"%s\"\n", masHelInfoFilenameWithPath);
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

}/*-------- END masFetchCouplingInfo()  ----------------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     masFetchFileList(void)                               /*--*/
/*--                                                              --*/
/*--  This function loads the list timesteps and filenames to use --*/
/*--  with MAS files.                                             --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  FILE * rfile;

  int max = 10000;
  char line[10000];

  int nFileLines;

  char *result = NULL;
  char delims[] = " \t\n\r";

  char masTimeFilename[MAX_STRING_SIZE] = "masTime.txt";
  char masTimeFilenameWithPath[MAX_STRING_SIZE];

  Scalar_t initialTime, timeVar;

  MPI_Aint N;
  MPI_Aint size;
  int disp_unit;

  // build the path to the time file
  sprintf(masTimeFilenameWithPath, "%s%s", config.masDirectory, masTimeFilename);

  // attempt to open the time file
  rfile = fopen(masTimeFilenameWithPath, "r");
  if (rfile==NULL) {
    printf("ERROR - Could not open file \"%s\"\n", masTimeFilenameWithPath);
    panic("Can't find MAS time list.");
  }

  // read the number of lines in the time file list
  nFileLines = 0;
  while (fgets(line, max, rfile) != NULL) {
    nFileLines++;
  }
  masTime = (Scalar_t *)malloc(sizeof(Scalar_t) * nFileLines);

  // reset the file pointer to the beginning of the file
  rewind(rfile);

  // read in the time file list
  for (int t=0; t<nFileLines; t++) {
    if (fgets(line, max, rfile) != NULL) {
      result = strtok(line, delims);
      if (result != NULL) {
        timeVar = (Scalar_t)atof(result) * MAS_TIME_CONVERT / DAY;
        if (t == 0) {initialTime = timeVar;}
        masTime[t] = config.masStartTime / DAY + (timeVar - initialTime);
      }
      result = strtok(NULL, delims);
    }
  }

  // close the file
  fclose(rfile);

  // store the number of files
  config.masNumFiles = nFileLines;

  result = NULL;

  // if coupling to the heliospheric domain
  if (config.masHelCouple > 0) {
        
    // build the path to the time file
    sprintf(masTimeFilenameWithPath, "%s%s", config.masHelDirectory, masTimeFilename);

    // attempt to open the time file
    rfile = fopen(masTimeFilenameWithPath, "r");
    if (rfile==NULL) {
      printf("ERROR - Could not open file \"%s\"\n", masTimeFilenameWithPath);
      panic("Can't find MAS Helio time list.");
    }

    // read the number of lines in the time file list
    nFileLines = 0;
    while (fgets(line, max, rfile) != NULL) {
      nFileLines++;
    }
    masHelTime = (Scalar_t *)malloc(sizeof(Scalar_t) * nFileLines);

    // reset the file pointer to the beginning of the file
    rewind(rfile);

    // read in the time file list
    for (int t=0; t<nFileLines; t++) {
      if (fgets(line, max, rfile) != NULL) {
        result = strtok(line, delims);
        if (result != NULL) {
          timeVar = (Scalar_t)atof(result) * MAS_TIME_CONVERT / DAY;
          if (t == 0) {initialTime = timeVar;}
          masHelTime[t] = config.masStartTime / DAY + (timeVar - initialTime);
        }
        result = strtok(NULL, delims);
      }
    }

    // close the file
    fclose(rfile);

    // store the number of files
    config.masHelNumFiles = nFileLines;

    result = NULL;

  }
  
  // Malloc arrays. (Only malloc once).
  if (masMallocFlag == 0) {

    // luckily, the size of the index arrays doesn't change in time
    masReadFieldIndex();
    // heliospheric coupling
    if (config.masHelCouple > 0)
      masHelReadFieldIndex();

    if(mpi_rank_shared==0){
      N=(int)masBprDimMax[0]*(int)masBptDimMax[0]*(int)masBppDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masBp_0, &masBp_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masBp_1, &masBp_1_win);
      N=(int)masBtrDimMax[0]*(int)masBttDimMax[0]*(int)masBtpDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masBt_0, &masBt_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masBt_1, &masBt_1_win);
      N=(int)masBrrDimMax[0]*(int)masBrtDimMax[0]*(int)masBrpDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masBr_0, &masBr_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masBr_1, &masBr_1_win);
      N=(int)masVprDimMax[0]*(int)masVptDimMax[0]*(int)masVppDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masVp_0, &masVp_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masVp_1, &masVp_1_win);
      N=(int)masVtrDimMax[0]*(int)masVttDimMax[0]*(int)masVtpDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masVt_0, &masVt_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masVt_1, &masVt_1_win);
      N=(int)masVrrDimMax[0]*(int)masVrtDimMax[0]*(int)masVrpDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masVr_0, &masVr_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masVr_1, &masVr_1_win);
      N=(int)masDrDimMax[0]*(int)masDtDimMax[0]*(int)masDpDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masD_0, &masD_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masD_1, &masD_1_win);
      
      // heliospheric coupling
      if (config.masHelCouple > 0) {
        
        N=(int)masHelBprDimMax[0]*(int)masHelBptDimMax[0]*(int)masHelBppDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBp_0, &masHelBp_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBp_1, &masHelBp_1_win);
        N=(int)masHelBtrDimMax[0]*(int)masHelBttDimMax[0]*(int)masHelBtpDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBt_0, &masHelBt_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBt_1, &masHelBt_1_win);
        N=(int)masHelBrrDimMax[0]*(int)masHelBrtDimMax[0]*(int)masHelBrpDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBr_0, &masHelBr_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBr_1, &masHelBr_1_win);
        N=(int)masHelVprDimMax[0]*(int)masHelVptDimMax[0]*(int)masHelVppDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVp_0, &masHelVp_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVp_1, &masHelVp_1_win);
        N=(int)masHelVtrDimMax[0]*(int)masHelVttDimMax[0]*(int)masHelVtpDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVt_0, &masHelVt_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVt_1, &masHelVt_1_win);
        N=(int)masHelVrrDimMax[0]*(int)masHelVrtDimMax[0]*(int)masHelVrpDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVr_0, &masHelVr_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVr_1, &masHelVr_1_win);
        N=(int)masHelDrDimMax[0]*(int)masHelDtDimMax[0]*(int)masHelDpDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelD_0, &masHelD_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelD_1, &masHelD_1_win);
        
      }

    }else{
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masBp_0, &masBp_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masBp_1, &masBp_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masBt_0, &masBt_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masBt_1, &masBt_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masBr_0, &masBr_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masBr_1, &masBr_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masVp_0, &masVp_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masVp_1, &masVp_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masVt_0, &masVt_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masVt_1, &masVt_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masVr_0, &masVr_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masVr_1, &masVr_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masD_0, &masD_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masD_1, &masD_1_win);
      MPI_Win_shared_query(masBp_0_win, 0, &size, &disp_unit, &masBp_0);
      MPI_Win_shared_query(masBt_0_win, 0, &size, &disp_unit, &masBt_0);
      MPI_Win_shared_query(masBr_0_win, 0, &size, &disp_unit, &masBr_0);
      MPI_Win_shared_query(masVp_0_win, 0, &size, &disp_unit, &masVp_0);
      MPI_Win_shared_query(masVt_0_win, 0, &size, &disp_unit, &masVt_0);
      MPI_Win_shared_query(masVr_0_win, 0, &size, &disp_unit, &masVr_0);
      MPI_Win_shared_query(masD_0_win, 0, &size, &disp_unit, &masD_0);
      MPI_Win_shared_query(masBp_1_win, 0, &size, &disp_unit, &masBp_1);
      MPI_Win_shared_query(masBt_1_win, 0, &size, &disp_unit, &masBt_1);
      MPI_Win_shared_query(masBr_1_win, 0, &size, &disp_unit, &masBr_1);
      MPI_Win_shared_query(masVp_1_win, 0, &size, &disp_unit, &masVp_1);
      MPI_Win_shared_query(masVt_1_win, 0, &size, &disp_unit, &masVt_1);
      MPI_Win_shared_query(masVr_1_win, 0, &size, &disp_unit, &masVr_1);
      MPI_Win_shared_query(masD_1_win, 0, &size, &disp_unit, &masD_1);

      // heliospheric coupling
      if (config.masHelCouple > 0) {
        
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBp_0, &masHelBp_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBp_1, &masHelBp_1_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBt_0, &masHelBt_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBt_1, &masHelBt_1_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBr_0, &masHelBr_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBr_1, &masHelBr_1_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVp_0, &masHelVp_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVp_1, &masHelVp_1_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVt_0, &masHelVt_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVt_1, &masHelVt_1_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVr_0, &masHelVr_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVr_1, &masHelVr_1_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelD_0, &masHelD_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelD_1, &masHelD_1_win);
        MPI_Win_shared_query(masHelBp_0_win, 0, &size, &disp_unit, &masHelBp_0);
        MPI_Win_shared_query(masHelBt_0_win, 0, &size, &disp_unit, &masHelBt_0);
        MPI_Win_shared_query(masHelBr_0_win, 0, &size, &disp_unit, &masHelBr_0);
        MPI_Win_shared_query(masHelVp_0_win, 0, &size, &disp_unit, &masHelVp_0);
        MPI_Win_shared_query(masHelVt_0_win, 0, &size, &disp_unit, &masHelVt_0);
        MPI_Win_shared_query(masHelVr_0_win, 0, &size, &disp_unit, &masHelVr_0);
        MPI_Win_shared_query(masHelD_0_win, 0, &size, &disp_unit, &masHelD_0);
        MPI_Win_shared_query(masHelBp_1_win, 0, &size, &disp_unit, &masHelBp_1);
        MPI_Win_shared_query(masHelBt_1_win, 0, &size, &disp_unit, &masHelBt_1);
        MPI_Win_shared_query(masHelBr_1_win, 0, &size, &disp_unit, &masHelBr_1);
        MPI_Win_shared_query(masHelVp_1_win, 0, &size, &disp_unit, &masHelVp_1);
        MPI_Win_shared_query(masHelVt_1_win, 0, &size, &disp_unit, &masHelVt_1);
        MPI_Win_shared_query(masHelVr_1_win, 0, &size, &disp_unit, &masHelVr_1);
        MPI_Win_shared_query(masHelD_1_win, 0, &size, &disp_unit, &masHelD_1);
        
      }

    }
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masBp_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masBt_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masBr_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masVp_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masVt_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masVr_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masD_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masBp_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masBt_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masBr_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masVp_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masVt_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masVr_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masD_1_win);

    // heliospheric coupling
    if (config.masHelCouple > 0) {
      
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelBp_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelBt_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelBr_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelVp_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelVt_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelVr_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelD_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelBp_1_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelBt_1_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelBr_1_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelVp_1_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelVt_1_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelVr_1_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelD_1_win);

    }
      
    masMallocFlag = 1;

  }


}/*-------- END masFetchFileList()  --------------------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     cleanupMPIWindows(void)                              /*--*/
/*--                                                              --*/
/*--  This function cleans up MPI windows.                        --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

   if (config.masCouple > 0) {
      MPI_Win_unlock_all(masBp_0_win);
      MPI_Win_unlock_all(masBt_0_win);
      MPI_Win_unlock_all(masBr_0_win);
      MPI_Win_unlock_all(masVp_0_win);
      MPI_Win_unlock_all(masVt_0_win);
      MPI_Win_unlock_all(masVr_0_win);
      MPI_Win_unlock_all(masD_0_win);
      MPI_Win_unlock_all(masBp_1_win);
      MPI_Win_unlock_all(masBt_1_win);
      MPI_Win_unlock_all(masBr_1_win);
      MPI_Win_unlock_all(masVp_1_win);
      MPI_Win_unlock_all(masVt_1_win);
      MPI_Win_unlock_all(masVr_1_win);
      MPI_Win_unlock_all(masD_1_win);

      MPI_Win_free(&masBp_0_win);
      MPI_Win_free(&masBt_0_win);
      MPI_Win_free(&masBr_0_win);
      MPI_Win_free(&masVp_0_win);
      MPI_Win_free(&masVt_0_win);
      MPI_Win_free(&masVr_0_win);
      MPI_Win_free(&masD_0_win);
      MPI_Win_free(&masBp_1_win);
      MPI_Win_free(&masBt_1_win);
      MPI_Win_free(&masBr_1_win);
      MPI_Win_free(&masVp_1_win);
      MPI_Win_free(&masVt_1_win);
      MPI_Win_free(&masVr_1_win);
      MPI_Win_free(&masD_1_win);
  }
  
  // heliospheric coupling
  if (config.masHelCouple > 0) {
    
    MPI_Win_unlock_all(masHelBp_0_win);
    MPI_Win_unlock_all(masHelBt_0_win);
    MPI_Win_unlock_all(masHelBr_0_win);
    MPI_Win_unlock_all(masHelVp_0_win);
    MPI_Win_unlock_all(masHelVt_0_win);
    MPI_Win_unlock_all(masHelVr_0_win);
    MPI_Win_unlock_all(masHelD_0_win);
    MPI_Win_unlock_all(masHelBp_1_win);
    MPI_Win_unlock_all(masHelBt_1_win);
    MPI_Win_unlock_all(masHelBr_1_win);
    MPI_Win_unlock_all(masHelVp_1_win);
    MPI_Win_unlock_all(masHelVt_1_win);
    MPI_Win_unlock_all(masHelVr_1_win);
    MPI_Win_unlock_all(masHelD_1_win);

    MPI_Win_free(&masHelBp_0_win);
    MPI_Win_free(&masHelBt_0_win);
    MPI_Win_free(&masHelBr_0_win);
    MPI_Win_free(&masHelVp_0_win);
    MPI_Win_free(&masHelVt_0_win);
    MPI_Win_free(&masHelVr_0_win);
    MPI_Win_free(&masHelD_0_win);
    MPI_Win_free(&masHelBp_1_win);
    MPI_Win_free(&masHelBt_1_win);
    MPI_Win_free(&masHelBr_1_win);
    MPI_Win_free(&masHelVp_1_win);
    MPI_Win_free(&masHelVt_1_win);
    MPI_Win_free(&masHelVr_1_win);
    MPI_Win_free(&masHelD_1_win);

    
  }

}/*-------- END cleanupMPIWindows()  --------------------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     masGetInterpData ( Scalar_t dt )                     /*--*/
/*--                                                              --*/
/*-- This function checks if the simulation time is right to read --*/
/*-- or copy MAS data, and does so if needed.                     --*/
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
// Check if MAS is being used yet.  
//
  if ((t_global < masTime[0]) && (config.masSteadyState == 0)){
    mhdGridStatus = MHD_DEFAULT;
    return;
  } 
  else{
    mhdGridStatus = MHD_MAS;
  }
  
//
// Find the two MAS files that bound the current time_interp, and set 
// interpolation factors.
//
  if (time_interp <= masTime[0]){
    
    masFileIndex0 = 0;
    masFileIndex1 = 0;

    s_cor = 0.0;
    
  } else if (time_interp >= masTime[config.masNumFiles-1]){
    
    masFileIndex0 = config.masNumFiles-1;
    masFileIndex1 = config.masNumFiles-1;
    
    s_cor = 0.0;
    
  } else {

    for (i=1; i<config.masNumFiles; i++){
      if (masTime[i] > time_interp ){
         masFileIndex0 = i-1;
         masFileIndex1 = i;
         break;
      } 
    }

    s_cor = ( time_interp - masTime[masFileIndex0] ) /
            ( masTime[masFileIndex1] - masTime[masFileIndex0]);
            
  }
  
  
//  
// Read or copy MAS data as needed.
//

  if (masFileIndex0 != masFileIndex_loaded0){
    if (masFileIndex0 == masFileIndex_loaded1){
      if(mpi_rank_shared==0) {   
        N=(int)masBprDimMax[0]*(int)masBptDimMax[0]*(int)masBppDimMax[0];
        for(i=0;i<N;i++){masBp_0[i] = masBp_1[i];}
        N=(int)masBtrDimMax[0]*(int)masBttDimMax[0]*(int)masBtpDimMax[0];
        for(i=0;i<N;i++){masBt_0[i] = masBt_1[i];}
        N=(int)masBrrDimMax[0]*(int)masBrtDimMax[0]*(int)masBrpDimMax[0];
        for(i=0;i<N;i++){masBr_0[i] = masBr_1[i];}
        N=(int)masVprDimMax[0]*(int)masVptDimMax[0]*(int)masVppDimMax[0];
        for(i=0;i<N;i++){masVp_0[i] = masVp_1[i];}
        N=(int)masVtrDimMax[0]*(int)masVttDimMax[0]*(int)masVtpDimMax[0];
        for(i=0;i<N;i++){masVt_0[i] = masVt_1[i];}
        N=(int)masVrrDimMax[0]*(int)masVrtDimMax[0]*(int)masVrpDimMax[0];
        for(i=0;i<N;i++){masVr_0[i] = masVr_1[i];}
        N=(int)masDrDimMax[0]*(int)masDtDimMax[0]*(int)masDpDimMax[0];
        for(i=0;i<N;i++){ masD_0[i] = masD_1[i];}
      }
    } else {
      if(mpi_rank_shared==0) { 
        masReadData(masFileIndex0,
                &masBp_0, &masBt_0, &masBr_0,
                &masVp_0, &masVt_0, &masVr_0,
                &masD_0);
      }
    }
    masFileIndex_loaded0 = masFileIndex0;
    need_sync_0 = 1;
  } else {
    need_sync_0 = 0;
  }
  
  if (masFileIndex1 != masFileIndex_loaded1){

    if(mpi_rank_shared==0) { 
      masReadData(masFileIndex1,
                &masBp_1, &masBt_1, &masBr_1,
                &masVp_1, &masVt_1, &masVr_1,
                &masD_1);
    }
    masFileIndex_loaded1 = masFileIndex1;  
    need_sync_1 = 1;
  } else {
    need_sync_1 = 0;
  }
    
  if (need_sync_0 == 1){  
    MPI_Win_sync(masBp_0_win);
    MPI_Win_sync(masBt_0_win);
    MPI_Win_sync(masBr_0_win);
    MPI_Win_sync(masVp_0_win);
    MPI_Win_sync(masVt_0_win);
    MPI_Win_sync(masVr_0_win);
    MPI_Win_sync(masD_0_win);
  }
  
  if (need_sync_1 == 1){  
    MPI_Win_sync(masBp_1_win);
    MPI_Win_sync(masBt_1_win);
    MPI_Win_sync(masBr_1_win);
    MPI_Win_sync(masVp_1_win);
    MPI_Win_sync(masVt_1_win);
    MPI_Win_sync(masVr_1_win);
    MPI_Win_sync(masD_1_win);
  }
  
  if (need_sync_0+need_sync_1 > 0){
    MPI_Barrier(comm_shared); 
  }
  
//
// Now do everything with the heliosphere.
//
  if (config.masHelCouple > 0) {
    
//
// Find the two MAS HELIO files that bound the current time_interp, and set 
// interpolation factors.
//
    if (time_interp < masHelTime[0]){
  	  
      masHelFileIndex0 = 0;
      masHelFileIndex1 = 0;

    } else if (time_interp > masHelTime[config.masHelNumFiles-1]){
  	  
      masHelFileIndex0 = config.masHelNumFiles-1;
      masHelFileIndex1 = config.masHelNumFiles-1;

    } else {
      
      for (i=1; i<config.masHelNumFiles; i++){
        if (masHelTime[i] > time_interp){
          masHelFileIndex0 = i-1;
          masHelFileIndex1 = i;
          break;
        } 
      }
    }
//  
// Read or copy MAS HELIO data as needed.
//
    if (masHelFileIndex0 != masHelFileIndex_loaded0){
      if (masHelFileIndex0 == masHelFileIndex_loaded1){
        if(mpi_rank_shared==0) {   
          N=(int)masHelBprDimMax[0]*(int)masHelBptDimMax[0]*(int)masHelBppDimMax[0];
          for(i=0;i<N;i++){masHelBp_0[i] = masHelBp_1[i];}
          N=(int)masHelBtrDimMax[0]*(int)masHelBttDimMax[0]*(int)masHelBtpDimMax[0];
          for(i=0;i<N;i++){masHelBt_0[i] = masHelBt_1[i];}
          N=(int)masHelBrrDimMax[0]*(int)masHelBrtDimMax[0]*(int)masHelBrpDimMax[0];
          for(i=0;i<N;i++){masHelBr_0[i] = masHelBr_1[i];}
          N=(int)masHelVprDimMax[0]*(int)masHelVptDimMax[0]*(int)masHelVppDimMax[0];
          for(i=0;i<N;i++){masHelVp_0[i] = masHelVp_1[i];}
          N=(int)masHelVtrDimMax[0]*(int)masHelVttDimMax[0]*(int)masHelVtpDimMax[0];
          for(i=0;i<N;i++){masHelVt_0[i] = masHelVt_1[i];}
          N=(int)masHelVrrDimMax[0]*(int)masHelVrtDimMax[0]*(int)masHelVrpDimMax[0];
          for(i=0;i<N;i++){masHelVr_0[i] = masHelVr_1[i];}
          N=(int)masHelDrDimMax[0]*(int)masHelDtDimMax[0]*(int)masHelDpDimMax[0];
          for(i=0;i<N;i++){ masHelD_0[i] = masHelD_1[i];}
        }
      } else {
        if(mpi_rank_shared==0) { 
          masHelReadData(masHelFileIndex0,
                  &masHelBp_0, &masHelBt_0, &masHelBr_0,
                  &masHelVp_0, &masHelVt_0, &masHelVr_0,
                  &masHelD_0);
        }
      }
      masHelFileIndex_loaded0 = masHelFileIndex0;
      need_sync_0 = 1;
    } else {
      need_sync_0 = 0;
    }
    
    if (masHelFileIndex1 != masHelFileIndex_loaded1){
  
      if(mpi_rank_shared==0) { 
        masHelReadData(masHelFileIndex1,
                  &masHelBp_1, &masHelBt_1, &masHelBr_1,
                  &masHelVp_1, &masHelVt_1, &masHelVr_1,
                  &masHelD_1);
      }
      masHelFileIndex_loaded1 = masHelFileIndex1;  
      need_sync_1 = 1;
    } else {
      need_sync_1 = 0;
    }
      
    if (need_sync_0 == 1){  
      MPI_Win_sync(masHelBp_0_win);
      MPI_Win_sync(masHelBt_0_win);
      MPI_Win_sync(masHelBr_0_win);
      MPI_Win_sync(masHelVp_0_win);
      MPI_Win_sync(masHelVt_0_win);
      MPI_Win_sync(masHelVr_0_win);
      MPI_Win_sync(masHelD_0_win);
    }
    
    if (need_sync_1 == 1){  
      MPI_Win_sync(masHelBp_1_win);
      MPI_Win_sync(masHelBt_1_win);
      MPI_Win_sync(masHelBr_1_win);
      MPI_Win_sync(masHelVp_1_win);
      MPI_Win_sync(masHelVt_1_win);
      MPI_Win_sync(masHelVr_1_win);
      MPI_Win_sync(masHelD_1_win);
    }
    
    if (need_sync_0+need_sync_1 > 0){
      MPI_Barrier(comm_shared); 
    }
    
  }
//
// Tell everyone what we just did.
//
  if (mpi_rank==0) {
    printf("  --> MAS(C) Couple: Time: %14.8e  Time+Dt: %14.8e  s_cor: %10.8f\n",
          t_global,
	        time_interp,s_cor);
    printf("  -->                idx0: %03d  idx1: %03d  masTime0: %14.8e  masTime1: %14.8e\n",
          masFileIndex0+1,
          masFileIndex1+1,
          masTime[masFileIndex0],
          masTime[masFileIndex1]);
    if (config.masHelCouple > 0){
      printf("  --> MAS(H) Couple: Time: %14.8e  Time+Dt: %14.8e  s_hel: %10.8f\n",
          t_global,
	        time_interp,s_hel);
      printf("  -->                idx0: %03d  idx1: %03d  masTime0: %14.8e  masTime1: %14.8e\n",
          masHelFileIndex0+1,
          masHelFileIndex1+1,
          masHelTime[masHelFileIndex0],
          masHelTime[masHelFileIndex1]);
    }      
  }

}
/*----------------- END masGetInterpData(dt)  ------------------------*/
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--*/     void                                                   /*--*/
/*--*/     masReadData(Index_t fileIndex,                         /*--*/
/*--*/            float *masBp[], float *masBt[], float *masBr[], /*--*/
/*--*/            float *masVp[], float *masVt[], float *masVr[], /*--*/
/*--*/            float *masD[])                                  /*--*/
/*--*/                                                            /*--*/
/*--                                                                --*/
/*--This function reads the MAS data from a HDF4 file.              --*/
/*--------------------------------------------------------------------*/
{/*-------------------------------------------------------------------*/

  int32 sd_id, sds_id;

  intn  status;

  int32 start[3] = {0,0,0};
  int32 edges[3];

  char fileNames[7][MAX_STRING_SIZE];

  double timer_tmp=0;

  timer_tmp = MPI_Wtime();

  if (config.masDigits == 3) {

    sprintf(fileNames[0], "%sbp%03d.hdf", config.masDirectory, fileIndex + 1);
    sprintf(fileNames[1], "%sbt%03d.hdf", config.masDirectory, fileIndex + 1);
    sprintf(fileNames[2], "%sbr%03d.hdf", config.masDirectory, fileIndex + 1);
    sprintf(fileNames[3], "%svp%03d.hdf",  config.masDirectory, fileIndex + 1);
    sprintf(fileNames[4], "%svt%03d.hdf",  config.masDirectory, fileIndex + 1);
    sprintf(fileNames[5], "%svr%03d.hdf",  config.masDirectory, fileIndex + 1);
    sprintf(fileNames[6], "%srho%03d.hdf",  config.masDirectory, fileIndex + 1);

  } else {

    sprintf(fileNames[0], "%sbp%06d.hdf", config.masDirectory, fileIndex + 1);
    sprintf(fileNames[1], "%sbt%06d.hdf", config.masDirectory, fileIndex + 1);
    sprintf(fileNames[2], "%sbr%06d.hdf", config.masDirectory, fileIndex + 1);
    sprintf(fileNames[3], "%svp%06d.hdf",  config.masDirectory, fileIndex + 1);
    sprintf(fileNames[4], "%svt%06d.hdf",  config.masDirectory, fileIndex + 1);
    sprintf(fileNames[5], "%svr%06d.hdf",  config.masDirectory, fileIndex + 1);
    sprintf(fileNames[6], "%srho%06d.hdf",  config.masDirectory, fileIndex + 1);

  }

  // reading in Bp
  edges[0] = masBppDimMax[0];
  edges[1] = masBptDimMax[0];
  edges[2] = masBprDimMax[0];

  sd_id = SDstart(fileNames[0], DFACC_READ);
  sds_id = SDselect (sd_id, 3);

  status = SDreaddata (sds_id, start, NULL, edges, (VOIDP)*masBp); ERR(status);
  status = SDendaccess (sds_id); ERR(status);
  status = SDend (sd_id); ERR(status);

  // reading in Bt
  edges[0] = masBtpDimMax[0];
  edges[1] = masBttDimMax[0];
  edges[2] = masBtrDimMax[0];

  sd_id = SDstart(fileNames[1], DFACC_READ);
  sds_id = SDselect (sd_id, 3);

  status = SDreaddata (sds_id, start, NULL, edges, (VOIDP)*masBt); ERR(status);
  status = SDendaccess (sds_id); ERR(status);
  status = SDend (sd_id); ERR(status);

  // reading in Br
  edges[0] = masBrpDimMax[0];
  edges[1] = masBrtDimMax[0];
  edges[2] = masBrrDimMax[0];

  sd_id = SDstart(fileNames[2], DFACC_READ);
  sds_id = SDselect (sd_id, 3);

  status = SDreaddata (sds_id, start, NULL, edges, (VOIDP)*masBr); ERR(status);
  status = SDendaccess (sds_id); ERR(status);
  status = SDend (sd_id); ERR(status);

  // reading in Vp
  edges[0] = masVppDimMax[0];
  edges[1] = masVptDimMax[0];
  edges[2] = masVprDimMax[0];

  sd_id = SDstart(fileNames[3], DFACC_READ);
  sds_id = SDselect (sd_id, 3);

  status = SDreaddata (sds_id, start, NULL, edges, (VOIDP)*masVp); ERR(status);
  status = SDendaccess (sds_id); ERR(status);
  status = SDend (sd_id); ERR(status);

  // reading in Vt
  edges[0] = masVtpDimMax[0];
  edges[1] = masVttDimMax[0];
  edges[2] = masVtrDimMax[0];

  sd_id = SDstart(fileNames[4], DFACC_READ);
  sds_id = SDselect (sd_id, 3);

  status = SDreaddata (sds_id, start, NULL, edges, (VOIDP)*masVt); ERR(status);
  status = SDendaccess (sds_id); ERR(status);
  status = SDend (sd_id); ERR(status);

  // reading in Vr
  edges[0] = masVrpDimMax[0];
  edges[1] = masVrtDimMax[0];
  edges[2] = masVrrDimMax[0];

  sd_id = SDstart(fileNames[5], DFACC_READ);
  sds_id = SDselect (sd_id, 3);

  status = SDreaddata (sds_id, start, NULL, edges, (VOIDP)*masVr); ERR(status);
  status = SDendaccess (sds_id); ERR(status);
  status = SDend (sd_id); ERR(status);

  // reading in D
  edges[0] = masDpDimMax[0];
  edges[1] = masDtDimMax[0];
  edges[2] = masDrDimMax[0];

  sd_id = SDstart(fileNames[6], DFACC_READ);
  sds_id = SDselect (sd_id, 3);

  status = SDreaddata (sds_id, start, NULL, edges, (VOIDP)*masD); ERR(status);
  status = SDendaccess (sds_id); ERR(status);
  status = SDend (sd_id); ERR(status);

  if (mpi_rank == 0) printf("  --> IO MAS: Coronal sequence %03d read.\n",fileIndex+1);

  timer_mas_io = timer_mas_io + (MPI_Wtime() - timer_tmp);

}/*-------- END masReadData()  -----------------------*/
/*------------------------------------------------------------------*/


/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--*/     void                                                   /*--*/
/*--*/     masHelReadData(Index_t fileIndex,                      /*--*/
/*--*/            float *masBp[], float *masBt[], float *masBr[], /*--*/
/*--*/            float *masVp[], float *masVt[], float *masVr[], /*--*/
/*--*/            float *masD[])                                  /*--*/
/*--*/                                                            /*--*/
/*--                                                                --*/
/*--This function checks if the simulation time is right to load    --*/
/*--another timestep file from the MAS data, and does if so.        --*/
/*--------------------------------------------------------------------*/
{/*-------------------------------------------------------------------*/

  int32 sd_id, sds_id;

  intn  status;

  int32 start[3] = {0,0,0};
  int32 edges[3];

  char fileNames[7][MAX_STRING_SIZE];

  double timer_tmp=0;

  timer_tmp = MPI_Wtime();

  if (config.masHelDigits == 3) {

    sprintf(fileNames[0], "%sbp%03d.hdf", config.masHelDirectory, fileIndex + 1);
    sprintf(fileNames[1], "%sbt%03d.hdf", config.masHelDirectory, fileIndex + 1);
    sprintf(fileNames[2], "%sbr%03d.hdf", config.masHelDirectory, fileIndex + 1);
    sprintf(fileNames[3], "%svp%03d.hdf",  config.masHelDirectory, fileIndex + 1);
    sprintf(fileNames[4], "%svt%03d.hdf",  config.masHelDirectory, fileIndex + 1);
    sprintf(fileNames[5], "%svr%03d.hdf",  config.masHelDirectory, fileIndex + 1);
    sprintf(fileNames[6], "%srho%03d.hdf",  config.masHelDirectory, fileIndex + 1);

  } else {

    sprintf(fileNames[0], "%sbp%06d.hdf", config.masHelDirectory, fileIndex + 1);
    sprintf(fileNames[1], "%sbt%06d.hdf", config.masHelDirectory, fileIndex + 1);
    sprintf(fileNames[2], "%sbr%06d.hdf", config.masHelDirectory, fileIndex + 1);
    sprintf(fileNames[3], "%svp%06d.hdf",  config.masHelDirectory, fileIndex + 1);
    sprintf(fileNames[4], "%svt%06d.hdf",  config.masHelDirectory, fileIndex + 1);
    sprintf(fileNames[5], "%svr%06d.hdf",  config.masHelDirectory, fileIndex + 1);
    sprintf(fileNames[6], "%srho%06d.hdf",  config.masHelDirectory, fileIndex + 1);

  }

  // reading in Bp
  edges[0] = masHelBppDimMax[0];
  edges[1] = masHelBptDimMax[0];
  edges[2] = masHelBprDimMax[0];

  sd_id = SDstart(fileNames[0], DFACC_READ);
  sds_id = SDselect (sd_id, 3);

  status = SDreaddata (sds_id, start, NULL, edges, (VOIDP)*masBp); ERR(status);
  status = SDendaccess (sds_id); ERR(status);
  status = SDend (sd_id); ERR(status);

  // reading in Bt
  edges[0] = masHelBtpDimMax[0];
  edges[1] = masHelBttDimMax[0];
  edges[2] = masHelBtrDimMax[0];

  sd_id = SDstart(fileNames[1], DFACC_READ);
  sds_id = SDselect (sd_id, 3);

  status = SDreaddata (sds_id, start, NULL, edges, (VOIDP)*masBt); ERR(status);
  status = SDendaccess (sds_id); ERR(status);
  status = SDend (sd_id); ERR(status);

  // reading in Br
  edges[0] = masHelBrpDimMax[0];
  edges[1] = masHelBrtDimMax[0];
  edges[2] = masHelBrrDimMax[0];

  sd_id = SDstart(fileNames[2], DFACC_READ);
  sds_id = SDselect (sd_id, 3);

  status = SDreaddata (sds_id, start, NULL, edges, (VOIDP)*masBr); ERR(status);
  status = SDendaccess (sds_id); ERR(status);
  status = SDend (sd_id); ERR(status);

  // reading in Vp
  edges[0] = masHelVppDimMax[0];
  edges[1] = masHelVptDimMax[0];
  edges[2] = masHelVprDimMax[0];

  sd_id = SDstart(fileNames[3], DFACC_READ);
  sds_id = SDselect (sd_id, 3);

  status = SDreaddata (sds_id, start, NULL, edges, (VOIDP)*masVp); ERR(status);
  status = SDendaccess (sds_id); ERR(status);
  status = SDend (sd_id); ERR(status);

  // reading in Vt
  edges[0] = masHelVtpDimMax[0];
  edges[1] = masHelVttDimMax[0];
  edges[2] = masHelVtrDimMax[0];

  sd_id = SDstart(fileNames[4], DFACC_READ);
  sds_id = SDselect (sd_id, 3);

  status = SDreaddata (sds_id, start, NULL, edges, (VOIDP)*masVt); ERR(status);
  status = SDendaccess (sds_id); ERR(status);
  status = SDend (sd_id); ERR(status);

  // reading in Vr
  edges[0] = masHelVrpDimMax[0];
  edges[1] = masHelVrtDimMax[0];
  edges[2] = masHelVrrDimMax[0];

  sd_id = SDstart(fileNames[5], DFACC_READ);
  sds_id = SDselect (sd_id, 3);

  status = SDreaddata (sds_id, start, NULL, edges, (VOIDP)*masVr); ERR(status);
  status = SDendaccess (sds_id); ERR(status);
  status = SDend (sd_id); ERR(status);

  // reading in D
  edges[0] = masHelDpDimMax[0];
  edges[1] = masHelDtDimMax[0];
  edges[2] = masHelDrDimMax[0];

  sd_id = SDstart(fileNames[6], DFACC_READ);
  sds_id = SDselect (sd_id, 3);

  status = SDreaddata (sds_id, start, NULL, edges, (VOIDP)*masD); ERR(status);
  status = SDendaccess (sds_id); ERR(status);
  status = SDend (sd_id); ERR(status);

  if (mpi_rank == 0) printf("  --> IO MAS: Helio sequence %03d read.\n",fileIndex+1);

  timer_mas_io = timer_mas_io + (MPI_Wtime() - timer_tmp);

}/*-------- END masHelReadData()  ----------------------------------*/
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     masReadFieldIndex(void)                              /*--*/
/*--                                                              --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  int32 sd_id, sds_id;

  intn status;

  int32 dim_sizes[H4_MAX_VAR_DIMS];
  int32 rank, data_type, n_attrs;
  char  name[H4_MAX_NC_NAME];

  char fileNames[7][MAX_STRING_SIZE];

  Scalar_t rMin, rMax, rTemp;

  double timer_tmp=0;

  timer_tmp = MPI_Wtime();

  if (config.masDigits == 3) {

    sprintf(fileNames[0], "%sbp%03d.hdf", config.masDirectory, 1);
    sprintf(fileNames[1], "%sbt%03d.hdf", config.masDirectory, 1);
    sprintf(fileNames[2], "%sbr%03d.hdf", config.masDirectory, 1);
    sprintf(fileNames[3], "%svp%03d.hdf", config.masDirectory, 1);
    sprintf(fileNames[4], "%svt%03d.hdf", config.masDirectory, 1);
    sprintf(fileNames[5], "%svr%03d.hdf", config.masDirectory, 1);
    sprintf(fileNames[6], "%srho%03d.hdf", config.masDirectory, 1);

  } else {

    sprintf(fileNames[0], "%sbp%06d.hdf", config.masDirectory, 1);
    sprintf(fileNames[1], "%sbt%06d.hdf", config.masDirectory, 1);
    sprintf(fileNames[2], "%sbr%06d.hdf", config.masDirectory, 1);
    sprintf(fileNames[3], "%svp%06d.hdf", config.masDirectory, 1);
    sprintf(fileNames[4], "%svt%06d.hdf", config.masDirectory, 1);
    sprintf(fileNames[5], "%svr%06d.hdf", config.masDirectory, 1);
    sprintf(fileNames[6], "%srho%06d.hdf", config.masDirectory, 1);

  }


  // masBp Dimensions
  sd_id = SDstart(fileNames[0], DFACC_READ);

  // masBppDim
  sds_id = SDselect (sd_id, 0);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masBppDimMax[0] = dim_sizes[0];
  masBppDim = (float *) malloc (sizeof(float) * (int)(masBppDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masBppDimMax, (VOIDP)masBppDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masBptDim
  sds_id = SDselect (sd_id, 1);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masBptDimMax[0] = dim_sizes[0];
  masBptDim = (float *) malloc (sizeof(float) * (int)(masBptDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masBptDimMax, (VOIDP)masBptDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masBprDim
  sds_id = SDselect (sd_id, 2);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masBprDimMax[0] = dim_sizes[0];
  masBprDim = (float *) malloc (sizeof(float) * (int)(masBprDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masBprDimMax, (VOIDP)masBprDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // grab the min and max for the r
  rMin = masBprDim[0];
  rMax = masBprDim[masBprDimMax[0] - 1];

  status = SDend (sd_id); ERR(status);

  // masBt Dimensions
  sd_id = SDstart(fileNames[1], DFACC_READ);

  // masBtpDim
  sds_id = SDselect (sd_id, 0);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masBtpDimMax[0] = dim_sizes[0];
  masBtpDim = (float *) malloc (sizeof(float) * (int)(masBtpDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masBtpDimMax, (VOIDP)masBtpDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masBttDim
  sds_id = SDselect (sd_id, 1);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masBttDimMax[0] = dim_sizes[0];
  masBttDim = (float *) malloc (sizeof(float) * (int)(masBttDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masBttDimMax, (VOIDP)masBttDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masBtrDim
  sds_id = SDselect (sd_id, 2);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masBtrDimMax[0] = dim_sizes[0];
  masBtrDim = (float *) malloc (sizeof(float) * (int)(masBtrDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masBtrDimMax, (VOIDP)masBtrDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // grab the min and max for the r
  rTemp = masBtrDim[0];
  if (rTemp > rMin)
    rMin = rTemp;

  rTemp = masBtrDim[masBtrDimMax[0] - 1];
  if (rTemp < rMax)
    rMax = rTemp;

  status = SDend (sd_id); ERR(status);

  // masBr Dimensions
  sd_id = SDstart(fileNames[2], DFACC_READ);

  // masBrpDim
  sds_id = SDselect (sd_id, 0);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masBrpDimMax[0] = dim_sizes[0];
  masBrpDim = (float *) malloc (sizeof(float) * (int)(masBrpDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masBrpDimMax, (VOIDP)masBrpDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masBrtDim
  sds_id = SDselect (sd_id, 1);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masBrtDimMax[0] = dim_sizes[0];
  masBrtDim = (float *) malloc (sizeof(float) * (int)(masBrtDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masBrtDimMax, (VOIDP)masBrtDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masBrrDim
  sds_id = SDselect (sd_id, 2);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masBrrDimMax[0] = dim_sizes[0];
  masBrrDim = (float *) malloc (sizeof(float) * (int)(masBrrDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masBrrDimMax, (VOIDP)masBrrDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // grab the min and max for the r
  rTemp = masBrrDim[0];
  if (rTemp > rMin)
    rMin = rTemp;

  rTemp = masBrrDim[masBrrDimMax[0] - 1];
  if (rTemp < rMax)
    rMax = rTemp;

  status = SDend (sd_id); ERR(status);


  // masVp Dimensions
  sd_id = SDstart(fileNames[3], DFACC_READ);

  // masVppDim
  sds_id = SDselect (sd_id, 0);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masVppDimMax[0] = dim_sizes[0];
  masVppDim = (float *) malloc (sizeof(float) * (int)(masVppDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masVppDimMax, (VOIDP)masVppDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masVptDim
  sds_id = SDselect (sd_id, 1);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masVptDimMax[0] = dim_sizes[0];
  masVptDim = (float *) malloc (sizeof(float) * (int)(masVptDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masVptDimMax, (VOIDP)masVptDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masVprDim
  sds_id = SDselect (sd_id, 2);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masVprDimMax[0] = dim_sizes[0];
  masVprDim = (float *) malloc (sizeof(float) * (int)(masVprDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masVprDimMax, (VOIDP)masVprDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // grab the min and max for the r
  rTemp = masVprDim[0];
  if (rTemp > rMin)
    rMin = rTemp;

  rTemp = masVprDim[masVprDimMax[0] - 1];
  if (rTemp < rMax)
    rMax = rTemp;

  status = SDend (sd_id); ERR(status);

  // masVt Dimensions
  sd_id = SDstart(fileNames[4], DFACC_READ);

  // masVtpDim
  sds_id = SDselect (sd_id, 0);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masVtpDimMax[0] = dim_sizes[0];
  masVtpDim = (float *) malloc (sizeof(float) * (int)(masVtpDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masVtpDimMax, (VOIDP)masVtpDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masVttDim
  sds_id = SDselect (sd_id, 1);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masVttDimMax[0] = dim_sizes[0];
  masVttDim = (float *) malloc (sizeof(float) * (int)(masVttDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masVttDimMax, (VOIDP)masVttDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masVtrDim
  sds_id = SDselect (sd_id, 2);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masVtrDimMax[0] = dim_sizes[0];
  masVtrDim = (float *) malloc (sizeof(float) * (int)(masVtrDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masVtrDimMax, (VOIDP)masVtrDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // grab the min and max for the r
  rTemp = masVtrDim[0];
  if (rTemp > rMin)
    rMin = rTemp;

  rTemp = masVtrDim[masVtrDimMax[0] - 1];
  if (rTemp < rMax)
    rMax = rTemp;

  status = SDend (sd_id); ERR(status);

  // masVr Dimensions
  sd_id = SDstart(fileNames[5], DFACC_READ);

  // masVrpDim
  sds_id = SDselect (sd_id, 0);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masVrpDimMax[0] = dim_sizes[0];
  masVrpDim = (float *) malloc (sizeof(float) * (int)(masVrpDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masVrpDimMax, (VOIDP)masVrpDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masVrtDim
  sds_id = SDselect (sd_id, 1);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masVrtDimMax[0] = dim_sizes[0];
  masVrtDim = (float *) malloc (sizeof(float) * (int)(masVrtDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masVrtDimMax, (VOIDP)masVrtDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masVrrDim
  sds_id = SDselect (sd_id, 2);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masVrrDimMax[0] = dim_sizes[0];
  masVrrDim = (float *) malloc (sizeof(float) * (int)(masVrrDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masVrrDimMax, (VOIDP)masVrrDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // grab the min and max for the r
  rTemp = masVrrDim[0];
  if (rTemp > rMin)
    rMin = rTemp;

  rTemp = masVrrDim[masVrrDimMax[0] - 1];
  if (rTemp < rMax)
    rMax = rTemp;

  status = SDend (sd_id); ERR(status);


  // masD Dimensions
  sd_id = SDstart (fileNames[6], DFACC_READ);

  // masDpDim
  sds_id = SDselect (sd_id, 0);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masDpDimMax[0] = dim_sizes[0];
  masDpDim = (float *) malloc (sizeof(float) * (int)(masDpDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masDpDimMax, (VOIDP)masDpDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masDtDim
  sds_id = SDselect (sd_id, 1);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masDtDimMax[0] = dim_sizes[0];
  masDtDim = (float *) malloc (sizeof(float) * (int)(masDtDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masDtDimMax, (VOIDP)masDtDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masDrDim
  sds_id = SDselect (sd_id, 2);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masDrDimMax[0] = dim_sizes[0];
  masDrDim = (float *) malloc (sizeof(float) * (int)(masDrDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masDrDimMax, (VOIDP)masDrDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // grab the min and max for the r
  rTemp = masDrDim[0];
  if (rTemp > rMin)
    rMin = rTemp;

  rTemp = masDrDim[masDrDimMax[0] - 1];
  if (rTemp < rMax)
    rMax = rTemp;

  status = SDend (sd_id); ERR(status);


  // set rScale, and masRadialMin/Max
  config.rScale = rMin * RSAU;
  config.masRadialMin = rMin;
  config.masRadialMax = rMax;

  timer_mas_io = timer_mas_io + (MPI_Wtime() - timer_tmp);

}/*-------- END masReadFieldIndex()  -----------------------*/
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     masHelReadFieldIndex(void)                              /*--*/
/*--                                                              --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  int32 sd_id, sds_id;

  intn status;

  int32 dim_sizes[H4_MAX_VAR_DIMS];
  int32 rank, data_type, n_attrs;
  char  name[H4_MAX_NC_NAME];

  char fileNames[7][MAX_STRING_SIZE];

  Scalar_t rMin, rMax, rTemp;

  double timer_tmp=0;

  timer_tmp = MPI_Wtime();

  if (config.masHelDigits == 3) {

    sprintf(fileNames[0], "%sbp%03d.hdf", config.masHelDirectory, 1);
    sprintf(fileNames[1], "%sbt%03d.hdf", config.masHelDirectory, 1);
    sprintf(fileNames[2], "%sbr%03d.hdf", config.masHelDirectory, 1);
    sprintf(fileNames[3], "%svp%03d.hdf", config.masHelDirectory, 1);
    sprintf(fileNames[4], "%svt%03d.hdf", config.masHelDirectory, 1);
    sprintf(fileNames[5], "%svr%03d.hdf", config.masHelDirectory, 1);
    sprintf(fileNames[6], "%srho%03d.hdf", config.masHelDirectory, 1);

  } else {

    sprintf(fileNames[0], "%sbp%06d.hdf", config.masHelDirectory, 1);
    sprintf(fileNames[1], "%sbt%06d.hdf", config.masHelDirectory, 1);
    sprintf(fileNames[2], "%sbr%06d.hdf", config.masHelDirectory, 1);
    sprintf(fileNames[3], "%svp%06d.hdf", config.masHelDirectory, 1);
    sprintf(fileNames[4], "%svt%06d.hdf", config.masHelDirectory, 1);
    sprintf(fileNames[5], "%svr%06d.hdf", config.masHelDirectory, 1);
    sprintf(fileNames[6], "%srho%06d.hdf", config.masHelDirectory, 1);

  }


  // masBp Dimensions
  sd_id = SDstart(fileNames[0], DFACC_READ);

  // masHelBppDim
  sds_id = SDselect (sd_id, 0);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelBppDimMax[0] = dim_sizes[0];
  masHelBppDim = (float *) malloc (sizeof(float) * (int)(masHelBppDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelBppDimMax, (VOIDP)masHelBppDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masHelBptDim
  sds_id = SDselect (sd_id, 1);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelBptDimMax[0] = dim_sizes[0];
  masHelBptDim = (float *) malloc (sizeof(float) * (int)(masHelBptDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelBptDimMax, (VOIDP)masHelBptDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masHelBprDim
  sds_id = SDselect (sd_id, 2);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelBprDimMax[0] = dim_sizes[0];
  masHelBprDim = (float *) malloc (sizeof(float) * (int)(masHelBprDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelBprDimMax, (VOIDP)masHelBprDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // grab the min and max for the r
  rMin = masHelBprDim[0];
  rMax = masHelBprDim[masHelBprDimMax[0] - 1];

  status = SDend (sd_id); ERR(status);

  // masHelBt Dimensions
  sd_id = SDstart(fileNames[1], DFACC_READ);

  // masHelBtpDim
  sds_id = SDselect (sd_id, 0);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelBtpDimMax[0] = dim_sizes[0];
  masHelBtpDim = (float *) malloc (sizeof(float) * (int)(masHelBtpDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelBtpDimMax, (VOIDP)masHelBtpDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masHelBttDim
  sds_id = SDselect (sd_id, 1);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelBttDimMax[0] = dim_sizes[0];
  masHelBttDim = (float *) malloc (sizeof(float) * (int)(masHelBttDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelBttDimMax, (VOIDP)masHelBttDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masHelBtrDim
  sds_id = SDselect (sd_id, 2);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelBtrDimMax[0] = dim_sizes[0];
  masHelBtrDim = (float *) malloc (sizeof(float) * (int)(masHelBtrDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelBtrDimMax, (VOIDP)masHelBtrDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // grab the min and max for the r
  rTemp = masHelBtrDim[0];
  if (rTemp > rMin)
    rMin = rTemp;

  rTemp = masHelBtrDim[masHelBtrDimMax[0] - 1];
  if (rTemp < rMax)
    rMax = rTemp;

  status = SDend (sd_id); ERR(status);

  // masHelBr Dimensions
  sd_id = SDstart(fileNames[2], DFACC_READ);

  // masHelBrpDim
  sds_id = SDselect (sd_id, 0);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelBrpDimMax[0] = dim_sizes[0];
  masHelBrpDim = (float *) malloc (sizeof(float) * (int)(masHelBrpDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelBrpDimMax, (VOIDP)masHelBrpDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masHelBrtDim
  sds_id = SDselect (sd_id, 1);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelBrtDimMax[0] = dim_sizes[0];
  masHelBrtDim = (float *) malloc (sizeof(float) * (int)(masHelBrtDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelBrtDimMax, (VOIDP)masHelBrtDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masHelBrrDim
  sds_id = SDselect (sd_id, 2);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelBrrDimMax[0] = dim_sizes[0];
  masHelBrrDim = (float *) malloc (sizeof(float) * (int)(masHelBrrDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelBrrDimMax, (VOIDP)masHelBrrDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // grab the min and max for the r
  rTemp = masHelBrrDim[0];
  if (rTemp > rMin)
    rMin = rTemp;

  rTemp = masHelBrrDim[masHelBrrDimMax[0] - 1];
  if (rTemp < rMax)
    rMax = rTemp;

  status = SDend (sd_id); ERR(status);


  // masHelVp Dimensions
  sd_id = SDstart(fileNames[3], DFACC_READ);

  // masHelVppDim
  sds_id = SDselect (sd_id, 0);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelVppDimMax[0] = dim_sizes[0];
  masHelVppDim = (float *) malloc (sizeof(float) * (int)(masHelVppDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelVppDimMax, (VOIDP)masHelVppDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masHelVptDim
  sds_id = SDselect (sd_id, 1);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelVptDimMax[0] = dim_sizes[0];
  masHelVptDim = (float *) malloc (sizeof(float) * (int)(masHelVptDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelVptDimMax, (VOIDP)masHelVptDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masHelVprDim
  sds_id = SDselect (sd_id, 2);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelVprDimMax[0] = dim_sizes[0];
  masHelVprDim = (float *) malloc (sizeof(float) * (int)(masHelVprDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelVprDimMax, (VOIDP)masHelVprDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // grab the min and max for the r
  rTemp = masHelVprDim[0];
  if (rTemp > rMin)
    rMin = rTemp;

  rTemp = masHelVprDim[masHelVprDimMax[0] - 1];
  if (rTemp < rMax)
    rMax = rTemp;

  status = SDend (sd_id); ERR(status);

  // masHelVt Dimensions
  sd_id = SDstart(fileNames[4], DFACC_READ);

  // masHelVtpDim
  sds_id = SDselect (sd_id, 0);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelVtpDimMax[0] = dim_sizes[0];
  masHelVtpDim = (float *) malloc (sizeof(float) * (int)(masHelVtpDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelVtpDimMax, (VOIDP)masHelVtpDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masHelVttDim
  sds_id = SDselect (sd_id, 1);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelVttDimMax[0] = dim_sizes[0];
  masHelVttDim = (float *) malloc (sizeof(float) * (int)(masHelVttDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelVttDimMax, (VOIDP)masHelVttDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masHelVtrDim
  sds_id = SDselect (sd_id, 2);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelVtrDimMax[0] = dim_sizes[0];
  masHelVtrDim = (float *) malloc (sizeof(float) * (int)(masHelVtrDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelVtrDimMax, (VOIDP)masHelVtrDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // grab the min and max for the r
  rTemp = masHelVtrDim[0];
  if (rTemp > rMin)
    rMin = rTemp;

  rTemp = masHelVtrDim[masHelVtrDimMax[0] - 1];
  if (rTemp < rMax)
    rMax = rTemp;

  status = SDend (sd_id); ERR(status);

  // masHelVr Dimensions
  sd_id = SDstart(fileNames[5], DFACC_READ);

  // masHelVrpDim
  sds_id = SDselect (sd_id, 0);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelVrpDimMax[0] = dim_sizes[0];
  masHelVrpDim = (float *) malloc (sizeof(float) * (int)(masHelVrpDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelVrpDimMax, (VOIDP)masHelVrpDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masHelVrtDim
  sds_id = SDselect (sd_id, 1);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelVrtDimMax[0] = dim_sizes[0];
  masHelVrtDim = (float *) malloc (sizeof(float) * (int)(masHelVrtDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelVrtDimMax, (VOIDP)masHelVrtDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masHelVrrDim
  sds_id = SDselect (sd_id, 2);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelVrrDimMax[0] = dim_sizes[0];
  masHelVrrDim = (float *) malloc (sizeof(float) * (int)(masHelVrrDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelVrrDimMax, (VOIDP)masHelVrrDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // grab the min and max for the r
  rTemp = masHelVrrDim[0];
  if (rTemp > rMin)
    rMin = rTemp;

  rTemp = masHelVrrDim[masHelVrrDimMax[0] - 1];
  if (rTemp < rMax)
    rMax = rTemp;

  status = SDend (sd_id); ERR(status);


  // masHelD Dimensions
  sd_id = SDstart (fileNames[6], DFACC_READ);

  // masHelDpDim
  sds_id = SDselect (sd_id, 0);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelDpDimMax[0] = dim_sizes[0];
  masHelDpDim = (float *) malloc (sizeof(float) * (int)(masHelDpDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelDpDimMax, (VOIDP)masHelDpDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masHelDtDim
  sds_id = SDselect (sd_id, 1);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelDtDimMax[0] = dim_sizes[0];
  masHelDtDim = (float *) malloc (sizeof(float) * (int)(masHelDtDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelDtDimMax, (VOIDP)masHelDtDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // masHelDrDim
  sds_id = SDselect (sd_id, 2);
  status = SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs);
  masHelDrDimMax[0] = dim_sizes[0];
  masHelDrDim = (float *) malloc (sizeof(float) * (int)(masHelDrDimMax[0]));
  status = SDreaddata (sds_id, masDimMin, NULL, masHelDrDimMax, (VOIDP)masHelDrDim); ERR(status);
  status = SDendaccess (sds_id); ERR(status);

  // grab the min and max for the r
  rTemp = masHelDrDim[0];
  if (rTemp > rMin)
    rMin = rTemp;

  rTemp = masHelDrDim[masHelDrDimMax[0] - 1];
  if (rTemp < rMax)
    rMax = rTemp;

  status = SDend (sd_id); ERR(status);


  // set rScale, and masHelRadialMin/Max
  config.masHelRadialMin = rMin;
  config.masHelRadialMax = rMax;

  timer_mas_io = timer_mas_io + (MPI_Wtime() - timer_tmp);

}/*-------- END masHelReadFieldIndex()  ----------------------------*/
/*------------------------------------------------------------------*/
