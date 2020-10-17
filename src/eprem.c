/*------------------------------------------------------------------------------
 ---- EPREM: main.c
 ----
 -- A spatial grid in 3-d space with grid points (nodes) associated with
 -- magnetic field lines anchored at the sun's surface and freely floating
 -- nodes. The spatial grid is housed in a data structure based on nested cubes
 -- whose surfaces are regularly subdivided into square arrays of square cells.
 -- Each cell defines a grid node (vertex) at its center. Each vertex has four
 -- surface neighbor vertices in the NEWS directions defined by its cell's
 -- edge adjacencies (which wrap at cube edges).
 -- A cube's collection of vertices and their surface neighbor relations
 -- (edges) defines a sub-graph of the entire grid of nested cubes.
 -- A cube houses a "shell" of grid nodes, the nested shells comprise the grid.
 -- Each node has a space position vector, r, defining its location
 -- in 3-d space. Initially, a shell sits at some fixed distance from the
 -- sun's surface. As time
 -- advances, the vertices' positions change with the locally defined flow
 -- field, moving outward from the sun's surface. But, graph adjacencies do
 -- not change.
 --
 -- Periodically, a new set of vertices is generated at the sun's surface.
 -- This causes the old vertex data (the shell) on the inner-most cube to move
 -- outward to the next enclosing cube in the data structure.
 -- That causes a ripple effect of shell data moving outward from cube to cube.
 -- For each proc, its inner-most cube serves as a communications buffer.
 -- Consequently, the inner-most shell on each proc resides on the next-to-
 -- inner-most cube (the inner-most cube just buffers a copy of that shell.)
 -- Each vertex then has two additional neighbor links besides its NEWS
 -- neighbors, streamIn and streamOut, identifying shell node neighbors in
 -- the next inner shell and next outer shell, respectively. (NB--These stream
 -- neighbor links connect the active grid shells. Thus, these stream links
 -- jump over the inner-most cube as that shell is not active in grid updates.)
 -- The "stream" notation indicates that these links are defined along
 -- magnetic field lines. The inner-most cube on the inner-most proc never
 -- has its grid node data updated by the main simulation loop as that shell
 -- serves as a template in spawning new shells.
 --
 -- NB-A good notation would be helpful for separating the notion of the cubic
 -- data structure's elements (vertices and graph edges) from the conceptual
 -- idea of the spatially located shell nodes and their associated data. The
 -- cube structure is static, while the shell nodes stream through it.
 -- Nested shells of nodes "live in" nested sub-divided cube sub-graphs, and
 -- shells migrate outward from cube to cube, eventually disappearing.
 -- Shells have a spatial distribution, initially on the surface of a sphere,
 -- while cubes have a logical 3-d structure consisting six orthogonal square
 -- grids (number of rows = number of columns) connected as if on the faces of
 -- a cube. Cube vertices form the data structure, while a shell node is a
 -- conceptual collection of data that resides temporarily in a cube vertex.
 -- Vertices have adjacent vertices, while shell nodes have neighbors. A
 -- shell node can find its neighbor by using graph edges to get to adjacent
 -- cube vertices, where it will find its neighbor node's information.
 -- Unfortunately, this notation is not adhered to in a consistant fashion,
 -- and a clean separation between the concepts of cubes and shells is not
 -- made.
 --
 -- ___________________CHANGE HISTORY_________________________
 --
 -- NOTE TO DEVS: BE SURE TO UPDATE VERSION NUMBER IN GLOBAL.H!
 --
 -- ## Version 1.8.0, 09/22/2020
 --       - Updated constants and MAS normalizations.
 --         Now, Omega shoudl match exactly.
 --       - Updated heliospheric coupling.
 --
 -- ## Version 1.7.9, 09/02/2020
 --       - Fixed bug in reading heliosphere MAS data.
 --
 -- ## Version 1.7.8, 08/13/2020
 --       - Refactored seed function to be more readable.
 --
 -- ## Version 1.7.7, 07/10/2020
 --       - Removed cross-field diffusion subcycling.
 --         Routine now back to original but with the proper
 --         scaling in ShellData.
 --
 -- ## Version 1.7.6, 06/23/2020
 --       - made some changes to driftVelocity so that it has
 --         consisten units (vd/c). I also added the charge/mass
 --         dependence of the driftVelocity.
 --       - Modified DriftShellData( ) to make consistent use of node.n.dl
 --         which have AU lengths as is the consistent scale for lengths
 --         in the code
 --
 -- ## Version 1.7.5, 06/23/2020
 --       - Fixed bug in energeticPartles.c, ShellData()
 --         This code used to produce lateral lengths
 --         such as grid[idx_frcs(face,row,col,shell)].n.dlPer
 --         and grid[idx_frcs(face,row,col,shell)].e.dl that were
 --         in code units. I have multipled r.x, r.y, r.z 
 --         and r1.x, r1.y, r1.z by config.rScale so that 
 --         the units for the neighbor distances (dl and dlPer)
 --         are in AU, the natural unit used for lengths throughout. 
 --         This fixes DiffuseShellData(shell, dt) so that the 
 --         cross-field diffusion is calculated appropriately. 
 --         Equations such as delN = dt_kper / (node.n.dlPer * node.n.dlPer + VERYSMALL)
 --         now have appropriate units (not mixed units, AU^2/(codunit)^2
 --         which essentially killed the cross-field diffusion
 --
 -- ### Version 1.7.4, 06/16/2020
 -- ### modified by NS
 --       - Added subcycling to cross-field diffusion at the 
 --         calculated stablility time step.
 --
 -- ### Version 1.7.3, 06/02/2020
 -- ### modified by RC
 --       - Changed the way the preEruptionDuration works.
 --         It now delays all particle advances until the CME 
 --         is about to erupt.  This makes the delayed runs behave the 
 --         same as the unstable intant runs.
 --
 -- ### Version 1.7.2, 06/01/2020
 -- ### modified by RC
 --       - Some bug fixes.
 --       - Reinitializes seed particles at start of epCalc.
 --
 -- ### Version 1.7.1, 05/25/2020
 -- ### modified by RC
 --       - Fixed where energectic particles are initialized.
 --       - Fixed seed test with above change.
 --       - Modified output a bit.
 --       - Commented out Monte-carlo seeding alg.  the updated RK4
 --         gets well within tolerance now so it is unnessesary.
 --       - Fixed issue that the mas coupling domain was not being rotated.
 --
 -- ### Version 1.7.0, 05/22/2020
 -- ### modified by RC,MY
 --       - Fixed BUG and did cleanup of RK4 advancing of nodes.
 --       - Cleaned up comments in main eprem file.
 --       - Overhauled MAS file reading and time step setting (in new setDt()).
 --       - Made gather and scatter operations asynchronous.
 --       - Fixed start and stop time to be the actual time (not 1970!).
 --       - Various renaming and cleanup.
 --
 -- ### Version 1.6.2, 04/30/2020
 -- ### modified by RC
 --       - Added heliospheric dshmin parameter.
 --       - BUG fixes for helio and file reading.
 --         Now testsuite works.
 --       - BUG fix for MAS IO. Removed incorrect pointer pointing.
 --         Now we read from file each time.  This changes testsuite data
 --         slightly.
 --       - Fixed bug in curl routine.  Phi was never set.
 --       - Updated MAS IO to only load files as needed, otherwise copy.
 --
 -- ### Version 1.6.1, 04/29/2020
 -- ### modified by MY
 --       - Updated some old MPI-1 to new standard.
 --         Code now compiles with GNU 8.
 --
 -- ### Version 1.6.0, 03/17/2020
 -- ### modified by MG,RC
 --       - Added heliospheric MAS coupling.
 -- 
 -- ### Version 1.5.0, 07/31/2019
 -- ### modified by RC
 --       - Added MPI shared arrays so now only one rank per node
 --         reads and stores mas data and all other can read it.
 --
 -- ### Version 1.4.1, 07/15/2019
 -- ### modified by MG
 --       - Fixed bug related to delayed IO output step.
 --
 -- ### Version 1.4.0, 06/25/2019
 -- ### modified by MG/RC
 --       - Focusing bug fixes:
 --         - Sign change (paper typo)
 --         - DT missing (alg bug)
 --       - Set default focusing limit to 1.0.
 --       - Modified deafults for trace-back radii.
 --       - Modified some input defaults to avoid issues.
 --       - Now allowed to use negative MFP exponent.
 --
 -- ### Version 1.3.1, 02/18/2019
 -- ### modified by MG
 --       - Bug fix - focusingLimit not being set.
 --
 -- ### Version 1.3.0, 02/13/2019
 -- ### modified by MG/RC
 --       - Added 4th order node advection.
 --       - Added subtime node advection.
 --       - Added backward node seeding.
 --       - Added flux limiter.
 --
 -- ### Version 1.2.0, 05/03/2018
 -- ### modified by MG
 --       - Major shock handeling overhall.
 --       - Changed mfp.
 --       - Changed seed scaling.
 --
 -- ### Version 1.1.2, 06/21/2017
 -- ### modified by RC,MG
 --       - Changed IO to output obs files named linearly
 --         and added streamMapping.txt output to give information
 --         about the stream identifications.
 --
 -- ### Version 1.1.1, 06/20/2017
 -- ### modified by RC,MG
 --       - Eliminated numShellsPerRank.
 --         Now the user specifies the resolution directly
 --         using the new numNodesPerStream input parameter.
 --         Now any number of processors up to the number of streams
 --         can be used for any run.
 --
 -- ### Version 1.1.0, 06/16/2017
 -- ### modified by RC,MG
 --       - Code no longer needs to be recompiled to use
 --         different paramters.
 --         Added numRows, numColums, numEnergySteps,
 --         numMuSteps, numShellsPerRank to input file.
 --         Removed Makefile.am.template and made universal
 --         Makefile.am.
 --      -  Added c99 standard flags to compilation options
 --         and made rpath options more portable.
 --      -  Updated MPI calls for ripple to be non-blocking.
 --
 -- ### Version 1.0.0, 04/19/2017
 -- ### modified by RC,MG
 --       - Fixed final bug that was changing solution based on
 --         processor count.  This problem was incorrectly
 --         calculating the curl for MAS runs.
 --       - Fixed major cross-physics bugs.
 --       - Added parallel IO for some EPREM output.
 --
 -- ### Version 0.9.0, 03/15/2017
 -- ### modified by RC
 --       - Fixed bug that was changing solution based on
 --         processor count.
 --
 -- ### Version 0.8.0, 02/20/2017
 -- ### modified by RC
 --       - Added version number to code.
 --
 -- ___________________END CHANGE HISTORY_____________________
 -----------------------------------------------------------------------*/

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

#include <time.h>
#include <math.h>
#include "global.h"
#include "configuration.h"
#include "energeticParticlesInit.h"
#include "energeticParticles.h"
#include "cubeShellInit.h"
#include "searchTypes.h"
#include "unifiedOutput.h"
#include "observerOutput.h"
#include "flow.h"
#include "simCore.h"
#include "timers.h"

/* Initialize all global timers. */
double timer_diffusestream=0;
double timer_adiabaticfocus=0;
double timer_adiabaticchange=0;
double timer_diffuseshell=0;
double timer_driftshell=0;
double timer_eptotal=0;
double timer_eprem_io=0;
double timer_init=0;
double timer_other=0;
double timer_step=0;
double timer_wall=0;
double timer_start=0;
double timer_MPIgatherscatter=0;
double timer_MPIsendrecv=0;

/*====================== MAIN =======================================*/
int main(int argc, char *argv[]) {

  int epInit;
  double timer_tmp;
  int rciter;
  
  // Initialize MPI
  initMPI(argc, argv);

  epInit=0;
  timer_tmp = 0;
  rciter=0;
  simStarted=0;

  // Record the starting MPI time and real time.
  if (mpi_rank == 0) time ( &start_time );
  timer_start = MPI_Wtime();

  // Read and set runtime parameters
  initGlobalParameters(argv[1]);

  // Initialize MPI Offsets for Gatherv and Scatterv
  initMPIOffsets();

  // Memory allocation for the main global variables
  allocateGlobalVariables();

  // Initialize MPI Types
  initMPITypes();

  // Initialize flags
  flagParamInit();

  // Initialize time variables
  timeInitialization();

  // Define the streams which will be computed
  defineComputeLines();

  // Set simulation parameters
  simCoreInit();

  // Intialize e-grid and ep attributes
  initEnergeticParticlesGrids();

 // Initialize MHD.
  updateMhd( config.tDel );

  // Initialize the cube / shell structure
  // and set node positions through backward integration.
  gridStructInit();

  // Create the names for the output files
  buildOutputNames();

  // Initialize unified netCDF output
  if (config.unifiedOutput > 0) initObserverDataNetCDF();

  // Initialize point observer netCDF output
  if (config.numObservers > 0) initPointObserverDataNetCDF();

  // Initialize netCDF file for domain output
  if (config.epremDomain > 0) initDomainDumpNetCDF();

  // Initialize netCDF file for unstructured domain output
  if (config.unstructuredDomain > 0) initUnstructuredDomainDumpNetCDF();

  // Write the simulation params out to the screen and file
  RunParamsOut();

  // Get time taken for initialization:
  timer_init = MPI_Wtime() - timer_start;

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ------------------  MAIN SIMULATION LOOP  ---------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  if(mpi_rank==0){
    printf(" \n");
    printf("*************************************************\n");
    printf("****SEEDING NODES********************************\n");
    printf("*************************************************\n");
  }

  do
  {
    timer_step = MPI_Wtime();
    
    rciter+=1;

    // The seeding of nodes takes TOTAL_NUM_SHELLS iterations.
    // Since the t_global might not be exactly simStartTimeDay, we use 
    // iteration count here and sync the t_global.
    
    if (rciter == TOTAL_NUM_SHELLS+1)
    {
      if (simStarted == 0){
        simStarted = 1;
        // Intialize seed population of energetic particles.
        // Since seed function is r-dependent, this needs to be done
        // after pushing all nodes out to their initial positions.
        initEnergeticParticles();
        if (mpi_rank == 0){
          printf(" \n");
          printf("*************************************************\n");
          printf("****STARTING SIMULATION**************************\n");
          printf("*************************************************\n");
          printf("NOTE:  Initializing seed population.\n");
          printf("NOTE:  Syncing t_global [%14.8e] to config.simStartTimeDay [%14.8e]\n", 
                 t_global, config.simStartTimeDay);
        }
        t_global = config.simStartTimeDay;
      }
    }
    
    if (t_global*DAY >= config.epCalcStartTime)
    {
      if (simStarted == 1){
        simStarted = 2;
        if (mpi_rank == 0){
          printf(" \n");
          printf("*************************************************\n");
          printf("****STARTING ENERGETIC PARTICLE CALCULATIONS*****\n");
          printf("*************************************************\n");
        }
      }
    }

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // ------ DISPLAY CURRENT TIME INFO TO SCREEN. -----------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    
    if (mpi_rank == 0){ 
      printf("Step: %06d  TIME: %14.8e  [JD %9.5f]  DTIME: %14.8e  [JD %7.5f]\n", 
             rciter, t_global, t_global*DAY, config.tDel, config.tDel*DAY);
    }

    // -------------------------------------------------------------------------    
    // -------------------------------------------------------------------------
    // -------------------  I/O OF CURRENT TIME MHD AND EPART ------------------
    // -------------------------------------------------------------------------    
    // -------------------------------------------------------------------------

    // NOTE!  This MUST be here and not below masGetInterpData.
    //        This is because the I/O (and moveNodes()) use getMasNode()
    //        which uses "s_xxx" and the MAS files from the previous step,
    //        which makes this work right because it makes it at the right time.

    if ( (num_loops % config.dumpFreq) == 0 )
    {
       dataDumpIO();
    }

    // -------------------------------------------------------------------------    
    // -------------------------------------------------------------------------
    // ---------------- MOVE NODES USING CURRENT TIME's MHD V ------------------
    // -------------------------------------------------------------------------    
    // -------------------------------------------------------------------------
      
    // Rotate the node seed positions and ripple the shells out
    rotSunAndSpawnShell( config.tDel );
    moveNodes( config.tDel );

    // -------------------------------------------------------------------------    
    // -------------------------------------------------------------------------
    // --------- UPDATE MHD FOR TIME+DT AND SAVE CURRENT MHD TO MHD-OLD --------
    // -------------------------------------------------------------------------    
    // -------------------------------------------------------------------------

    // NOTE:  updateMhd only uses tDel for computing div-V and it is never used.
    updateMhd( config.tDel );

    // -------------------------------------------------------------------------    
    // -------------------------------------------------------------------------
    // ---------------- UPDATE values used for perpendicular solvers -----------
    // -------------------------------------------------------------------------    
    // -------------------------------------------------------------------------   
    
    if ( (config.useDrift > 0) || (config.useShellDiffusion > 0) )
    {
      ShellData();
    }

    // -------------------------------------------------------------------------    
    // -------------------------------------------------------------------------
    // ----------------  PARTICLE COMPUTATION ----------------------------------
    // -------------------------------------------------------------------------    
    // -------------------------------------------------------------------------
    
    if ( (t_global*DAY >= config.epCalcStartTime) )
    {
      timer_tmp = MPI_Wtime();
      
      // Re-intialize seed population of energetic particles.
      // Since nodes are pushed before this point, this needs to be done
      // here.
      if (epInit == 0){
        epInit = 1;
        if (mpi_rank == 0) printf("NOTE:  Reinitializing seed population.\n");
        initEnergeticParticles();
      }
      
      updateEnergeticParticles();
      // For the seed test, re-init seed population.
      if (config.seedFunctionTest > 0){
        initEnergeticParticles();
        if (mpi_rank == 0) printf("NOTE:  Reinitializing seed population.\n");
      }
      
      timer_eptotal = timer_eptotal + (MPI_Wtime() - timer_tmp);
    }

    // -------------------------------------------------------------------------    
    // -------------------------------------------------------------------------
    // ----------------  INCREMENT TIME->TIME+DT--------------------------------
    // -------------------------------------------------------------------------    
    // -------------------------------------------------------------------------
    
    t_global += config.tDel;
    
    // IO LOOP counter.
    num_loops++;
    
    timer_tmp = MPI_Wtime();
    if (mpi_rank == 0) printf("  --> Compute time for step: %18.4f seconds.\n",timer_tmp-timer_step);
    
    
  } while(t_global <= config.simStopTimeDay);
  
  //
  // Dump final solution state.
  //
  //dataDumpIO();

  // ---------------------------------------------------------------------------    
  // ---------------------------------------------------------------------------
  // ------------------ CLEAN UP AND EXIT --------------------------------------
  // ---------------------------------------------------------------------------    
  // ---------------------------------------------------------------------------
  
  DumpRunTimes();
  config_destroy(&cfg);
  cleanupMPIWindows();
  MPI_Finalize();
  
  return(0);

}
//============================================================================
