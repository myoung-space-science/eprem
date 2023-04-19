/*-----------------------------------------------
 -- EMMREM: cubeShellInit.c
 --
 -- Routines to initialize the grid's cube-shell data structure and
 -- ancillaries.
 --
 -- ______________CHANGE HISTORY______________
 --
 -- ______________END CHANGE HISTORY______________
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
#include <math.h>
#include <time.h>
#include <float.h>

#include "global.h"
#include "cubeShellInit.h"
#include "cubeShellStruct.h"
#include "observerOutput.h"
#include "searchTypes.h"
#include "configuration.h"
#include "energeticParticlesBoundary.h"
#include "geometry.h"
#include "error.h"

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/          void                                                /*---*/
/*--*/    gridStructInit( void )                                    /*---*/
/*--*                                                                *---*/
/*--* Grid Structure Initialization. Data structure (grid) links to  *---*/
/*--* indicate neighbors of each grid cell are set to form nested    *---*/
/*--* gridded cube surfaces. Cube cells get their spatial positions  *---*/
/*--* set to points on the unit cube at the origin.                  *---*/
/*--* Cube cells then get there spatial positions projected onto     *---*/
/*--* the unit sphere. Cube cells are thus associated spatially      *---*/
/*--* with sphere nodes, which carry the simulation information.     *---*/
/*--* Initially, all spheres are identically located on the surface  *---*/
/*--* of the unit sphere.                                            *---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{

  /*-- Give each array element on the inner-most cube its cube --*/
  /*-- surface grid neighbor links.                            --*/

  initNEWS( INNER_SHELL );

  /*-- Give each array element of the inner-most cube its        --*/
  /*-- unit cube coordinates. Cube center is at origin (0,0,0).  --*/
  /*-- Array elements are aligned at centers of cells of the     --*/
  /*-- cube-face grids.                                          --*/

  initCubeCoords( INNER_SHELL );

  /*-- Give each array element of the inner-most cube its        --*/
  /*-- unit sphere coordinates in both (x,y,z) cartesian and     --*/
  /*-- (r, azimuth, zenith) polar coordinates.                   --*/

  initSphereCoords( INNER_SHELL );

  initCopyAll();
  initStreamNeighbors();

  checkNEWS( INNER_SHELL );

  if (mpi_rank == 0) initInnerShellValues();

} /*------ END  GridStructInit( ) --------------------*/
/*----------------------------------------------------*/

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/          void                                                /*---*/
/*--*/    initNEWS( Index_t shell )                                 /*---*/
/*--*                                                                *---*/
/*--*  Give each node on this shell its cube-surface grid links to   *---*/
/*--*  its NEWS neighbors on the cube surface.                       *---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{
  Index_t face, row, col;
  Index_t face2, row2, col2;
  Index_t col2p;
  Index_t i;
  Scalar_t dh, dw;

  dw = 1.0/(1.0*FACE_COLS);
  dh = 1.0/(1.0*FACE_ROWS);

  /*-- Do all the NEWS neighbors for all interior nodes.--*/
  /*-- These nodes are not on cube edges, thus they do  --*/
  /*-- not need any special treatment as their NEWS     --*/
  /*-- neigbors are just their array neigbors.          --*/
  /*-- We also take care of the boundary cases in part, --*/
  /*-- but we will have to fix up the connections that  --*/
  /*-- wrap around to other faces. For now we ignore    --*/
  /*-- that problem and assign all neighbors naively.   --*/

  for (face = 0; face < NUM_FACES; face++){
    for (row = 0; row < FACE_ROWS; row++){
      for (col = 0; col < FACE_COLS; col++){

        grid[idx_frcs(face,row,col,shell)].n.face  = face;
        grid[idx_frcs(face,row,col,shell)].n.row   = row - 1;
        grid[idx_frcs(face,row,col,shell)].n.col   = col;
        grid[idx_frcs(face,row,col,shell)].n.shell = shell;
        grid[idx_frcs(face,row,col,shell)].n.rank  = mpi_rank ;

        grid[idx_frcs(face,row,col,shell)].e.face  = face;
        grid[idx_frcs(face,row,col,shell)].e.row   = row;
        grid[idx_frcs(face,row,col,shell)].e.col   = col + 1;
        grid[idx_frcs(face,row,col,shell)].e.shell = shell;
        grid[idx_frcs(face,row,col,shell)].e.rank  = mpi_rank ;

        grid[idx_frcs(face,row,col,shell)].w.face  = face;
        grid[idx_frcs(face,row,col,shell)].w.row   = row;
        grid[idx_frcs(face,row,col,shell)].w.col   = col - 1;
        grid[idx_frcs(face,row,col,shell)].w.shell = shell;
        grid[idx_frcs(face,row,col,shell)].w.rank  = mpi_rank ;

        grid[idx_frcs(face,row,col,shell)].s.face  = face;
        grid[idx_frcs(face,row,col,shell)].s.row   = row + 1;
        grid[idx_frcs(face,row,col,shell)].s.col   = col;
        grid[idx_frcs(face,row,col,shell)].s.shell = shell;
        grid[idx_frcs(face,row,col,shell)].s.rank  = mpi_rank ;

	    }/*for col*/
    }/*for row*/
  }/*for face*/
  /*--------------------------------------------------------- */
  /*-- Now for the 12 special cases of cube edges.          --*/
  /*-- Consider the problem of unfolding a cardboard box:   --*/
  /*--                                                      --*/
  /*--              LEFT                                    --*/
  /*--               ::                                     --*/
  /*--    BACK  ::   TOP   ::  FRONT  ::   BOTTOM           --*/
  /*--               ::                                     --*/
  /*--              RIGHT                                   --*/
  /*--                                                      --*/
  /*-- Each outward face is face up in the above layout.    --*/
  /*-- The "::" show the face connections that are not cut. --*/
  /*-- These would be cube edges. Assuming the words were   --*/
  /*-- printed to appear upright when the box is assembled, --*/
  /*-- BACK should be shown above rotated 90 degrees clock- --*/
  /*-- wise, LEFT should be rotated by 180, and FRONT       --*/
  /*-- should be be rotated by 90 counter-clockwise.        --*/
  /*-- This corresponds to the array convention we adopted  --*/
  /*-- above. That is, looking at any face above from the   --*/
  /*-- point of view where the lettering appears upright,   --*/
  /*-- the upper-left corner of every face is cell [0,0]    --*/
  /*-- (assuming the indicated rotations have been done).   --*/
  /*--------------------------------------------------------- */

  /*-- connect bottom/S of TOP to top/N of RIGHT --*/
  face  = TOP_FACE;
  row   = FACE_ROWS - 1;      /*-- bottom of TOP. --*/
  face2 = RIGHT_FACE;
  row2  = 0;                  /*-- top of RIGHT. --*/
  for (col = 0; col < FACE_COLS; col++){
    col2 = col;  /*-- same x alignment. --*/
    grid[idx_frcs(face,row,col,shell)].s.face  = face2;
    grid[idx_frcs(face,row,col,shell)].s.row   = row2;
    grid[idx_frcs(face,row,col,shell)].s.col   = col2;
    grid[idx_frcs(face,row,col,shell)].s.rank   = mpi_rank ;

    grid[idx_frcs(face2,row2,col2,shell)].n.face  = face;
    grid[idx_frcs(face2,row2,col2,shell)].n.row   = row;
    grid[idx_frcs(face2,row2,col2,shell)].n.col   = col;
    grid[idx_frcs(face2,row2,col2,shell)].n.rank   = mpi_rank ;
  }

  /*-- connect top/N of TOP to top/N of LEFT --*/
  face  = TOP_FACE;
  row   = 0;                  /*-- top of TOP. --*/
  face2 = LEFT_FACE;
  row2  = 0;                  /*-- top of LEFT. --*/
  for (col = 0; col < FACE_COLS; col++){
    col2 = (FACE_COLS - 1) - col;  /*-- counter x alignment. --*/
    grid[idx_frcs(face,row,col,shell)].n.face  = face2;
    grid[idx_frcs(face,row,col,shell)].n.row   = row2;
    grid[idx_frcs(face,row,col,shell)].n.col   = col2;
    grid[idx_frcs(face,row,col,shell)].n.rank   = mpi_rank;

    grid[idx_frcs(face2,row2,col2,shell)].n.face  = face;
    grid[idx_frcs(face2,row2,col2,shell)].n.row   = row;
    grid[idx_frcs(face2,row2,col2,shell)].n.col   = col;
    grid[idx_frcs(face2,row2,col2,shell)].n.rank   = mpi_rank;
  }

  /*-- connect left/W of TOP to top/N of BACK --*/
  face  = TOP_FACE;
  col   = 0;                  /*-- left of TOP. --*/
  face2 = BACK_FACE;
  row2  = 0;                  /*-- top of BACK. --*/

  for (row = 0; row < FACE_ROWS; row++){
    col2 = (Index_t) ((0.5*dh/dw) - 0.5 + (row*dh/dw) + 0.499) ;
    /*-- y alignment: BACK col to TOP row. --*/
    grid[idx_frcs(face,row,col,shell)].w.face  = face2;
    grid[idx_frcs(face,row,col,shell)].w.row   = row2;
    grid[idx_frcs(face,row,col,shell)].w.col   = col2;
    grid[idx_frcs(face,row,col,shell)].w.rank   = mpi_rank;
  }

  for (col2 = 0; col2 < FACE_COLS; col2++) {
    row = (Index_t) ((0.5*dw/dh) - 0.5 + (col2*dw/dh) + 0.499) ;
    grid[idx_frcs(face2,row2,col2,shell)].n.face  = face;
    grid[idx_frcs(face2,row2,col2,shell)].n.row   = row;
    grid[idx_frcs(face2,row2,col2,shell)].n.col   = col;
    grid[idx_frcs(face2,row2,col2,shell)].n.rank   = mpi_rank;
  }

  /*-- connect right/E of TOP to top/N of FRONT --*/
  face  = TOP_FACE;
  col   = FACE_COLS - 1;      /*-- right of TOP. --*/
  face2 = FRONT_FACE;
  row2  = 0;                  /*-- top of FRONT. --*/
  for (row = 0; row < FACE_ROWS; row++){
    col2p = (Index_t) ((0.5*dh/dw) - 0.5 + (row*dh/dw) + 0.499) ;
    col2 = (FACE_COLS - 1) - col2p;  /*-- y: FRONT -col to TOP row. --*/
    grid[idx_frcs(face,row,col,shell)].e.face  = face2;
    grid[idx_frcs(face,row,col,shell)].e.row   = row2;
    grid[idx_frcs(face,row,col,shell)].e.col   = col2;
    grid[idx_frcs(face,row,col,shell)].e.rank   = mpi_rank;
  }

  for (col2p = 0; col2p < FACE_COLS; col2p++){
    col2 = (FACE_COLS - 1) - col2p; /* col2 runs from high to low number */
    row = (Index_t) ((0.5*dw/dh) - 0.5 + (col2p*dw/dh) + 0.499) ; /* row runs from low to high */
    grid[idx_frcs(face2,row2,col2,shell)].n.face  = face;
    grid[idx_frcs(face2,row2,col2,shell)].n.row   = row;
    grid[idx_frcs(face2,row2,col2,shell)].n.col   = col;
    grid[idx_frcs(face2,row2,col2,shell)].n.rank   = mpi_rank;
  }

  /*-- connect left/W of BOTTOM to bottom/S of FRONT --*/
  face  = BOTTOM_FACE;
  row   = 0;                  /*-- row 0 => row +i =  -y direction. --*/
  col   = 0;                  /*-- left of BOTTOM. --*/

  face2 = FRONT_FACE;
  col2  = FACE_COLS - 1;      /*-- col k => col k-i =  -y direction. --*/
  row2  = FACE_ROWS - 1;      /*-- bottom of FRONT. --*/

  for (row = 0; row < FACE_ROWS; row++){
    col2p = (Index_t) ((0.5*dh/dw) - 0.5 + (row*dh/dw) + 0.499) ;
    col2 = (FACE_COLS - 1) - col2p;
    grid[idx_frcs(face,row,col,shell)].w.face  = face2;
    grid[idx_frcs(face,row,col,shell)].w.row   = row2;
    grid[idx_frcs(face,row,col,shell)].w.col   = col2;
    grid[idx_frcs(face,row,col,shell)].w.rank   = mpi_rank;
  }

  /* run through columns on Front Face */
  for (col2p=0; col2p < FACE_COLS; col2p++){
    col2 = (FACE_COLS - 1) - col2p;  /* col 2 runs high to low */
    row = (Index_t) ((0.5*dw/dh) - 0.5 + (col2p*dw/dh) + 0.499) ; /* row runs from low to high */
    grid[idx_frcs(face2,row2,col2,shell)].s.face  = face;
    grid[idx_frcs(face2,row2,col2,shell)].s.row   = row;
    grid[idx_frcs(face2,row2,col2,shell)].s.col   = col;
    grid[idx_frcs(face2,row2,col2,shell)].s.rank  = mpi_rank;
  }

  /*-- connect right/E of BOTTOM to bottom/S of BACK --*/
  face  = BOTTOM_FACE;
  row   = 0;                  /*-- row 0 => row +i =  -y direction. --*/
  col   = FACE_COLS - 1;      /*-- right of BOTTOM. --*/

  face2 = BACK_FACE;
  col2  = 0;                  /*-- col 0 => col +i =  -y direction. --*/
  row2  = FACE_ROWS - 1;      /*-- bottom of BACK. --*/

  for (row = 0; row < FACE_ROWS; row++){
    col2 = (Index_t) ((0.5*dh/dw) - 0.5 + (row*dh/dw) + 0.499) ;
    grid[idx_frcs(face,row,col,shell)].e.face  = face2;
    grid[idx_frcs(face,row,col,shell)].e.row   = row2;
    grid[idx_frcs(face,row,col,shell)].e.col   = col2;
    grid[idx_frcs(face,row,col,shell)].e.rank   = mpi_rank;
  }

  for (col2 = 0; col2 < FACE_COLS; col2++){
    row = (Index_t) ((0.5*dw/dh) - 0.5 + (col2*dw/dh) + 0.499) ;
    grid[idx_frcs(face2,row2,col2,shell)].s.face  = face;
    grid[idx_frcs(face2,row2,col2,shell)].s.row   = row;
    grid[idx_frcs(face2,row2,col2,shell)].s.col   = col;
    grid[idx_frcs(face2,row2,col2,shell)].s.rank   = mpi_rank;
  }

  /*-- connect top/N of BOTTOM to bottom/S of LEFT --*/
  face  = BOTTOM_FACE;
  row   = 0;                  /*-- top of BOTTOM.                   --*/
  col   = 0;                  /*-- col 0 => col +i =  -x direction. --*/

  face2 = LEFT_FACE;
  row2  = FACE_ROWS - 1;      /*-- bottom of LEFT. --*/
  col2  = 0;                  /*-- col 0 => col +i =  -x direction. --*/

  for (i = 0; i < FACE_COLS; i++){
    grid[idx_frcs(face,row,col,shell)].n.face  = face2;
    grid[idx_frcs(face,row,col,shell)].n.row   = row2;
    grid[idx_frcs(face,row,col,shell)].n.col   = col2;
    grid[idx_frcs(face,row,col,shell)].n.rank   = mpi_rank;

    grid[idx_frcs(face2,row2,col2,shell)].s.face  = face;
    grid[idx_frcs(face2,row2,col2,shell)].s.row   = row;
    grid[idx_frcs(face2,row2,col2,shell)].s.col   = col;
    grid[idx_frcs(face2,row2,col2,shell)].s.rank   = mpi_rank;

    col++;   /*-- move -y direction for BOTTOM. --*/
    col2++;  /*-- move -y direction for LEFT. --*/
  }

  /*-- connect bottom/S of BOTTOM to bottom/S of RIGHT --*/
  face  = BOTTOM_FACE;
  row   = FACE_ROWS - 1;      /*-- bottom of BOTTOM.               --*/
  col   = 0;                  /*-- col 0 => col +i  =  -x direction. --*/

  face2 = RIGHT_FACE;
  row2  = FACE_ROWS - 1;      /*-- bottom of RIGHT. --*/
  col2  = FACE_COLS - 1;      /*-- col k => col k-i =  -x direction. --*/

  for (i = 0; i < FACE_COLS; i++){
    grid[idx_frcs(face,row,col,shell)].s.face  = face2;
    grid[idx_frcs(face,row,col,shell)].s.row   = row2;
    grid[idx_frcs(face,row,col,shell)].s.col   = col2;
    grid[idx_frcs(face,row,col,shell)].s.rank  = mpi_rank;

    grid[idx_frcs(face2,row2,col2,shell)].s.face  = face;
    grid[idx_frcs(face2,row2,col2,shell)].s.row   = row;
    grid[idx_frcs(face2,row2,col2,shell)].s.col   = col;
    grid[idx_frcs(face2,row2,col2,shell)].s.rank  = mpi_rank;

    col++;   /*-- move -y direction for BOTTOM. --*/
    col2--;  /*-- move -y direction for RIGHT. --*/
  }

  /*-- connect right/E of FRONT to left/W of LEFT --*/
  face  = FRONT_FACE;
  row   = 0;                  /*-- row 0 => row +i  = -z direction. --*/
  col   = FACE_COLS - 1;      /*-- right of FRONT. --*/

  face2 = LEFT_FACE;
  row2  = 0;                  /*-- row 0 => row +i  = -z direction. --*/
  col2  = 0;                  /*-- left of LEFT. --*/

  for (i = 0; i < FACE_ROWS; i++){
    grid[idx_frcs(face,row,col,shell)].e.face  = face2;
    grid[idx_frcs(face,row,col,shell)].e.row   = row2;
    grid[idx_frcs(face,row,col,shell)].e.col   = col2;
    grid[idx_frcs(face,row,col,shell)].e.rank   = mpi_rank;

    grid[idx_frcs(face2,row2,col2,shell)].w.face  = face;
    grid[idx_frcs(face2,row2,col2,shell)].w.row   = row;
    grid[idx_frcs(face2,row2,col2,shell)].w.col   = col;
    grid[idx_frcs(face2,row2,col2,shell)].w.rank  = mpi_rank;

    row++;   /*-- move -z direction for FRONT. --*/
    row2++;  /*-- move -z direction for LEFT. --*/
  }

  /*-- connect right/E of LEFT to left/W of BACK --*/
  face  = LEFT_FACE;
  row   = 0;                  /*-- row 0 => row +i  = -z direction. --*/
  col   = FACE_COLS - 1;      /*-- right. --*/

  face2 = BACK_FACE;
  row2  = 0;                  /*-- row 0 => row +i  = -z direction. --*/
  col2  = 0;                  /*-- left. --*/

  for (i = 0; i < FACE_ROWS; i++){
    grid[idx_frcs(face,row,col,shell)].e.face  = face2;
    grid[idx_frcs(face,row,col,shell)].e.row   = row2;
    grid[idx_frcs(face,row,col,shell)].e.col   = col2;
    grid[idx_frcs(face,row,col,shell)].e.rank   = mpi_rank;

    grid[idx_frcs(face2,row2,col2,shell)].w.face  = face;
    grid[idx_frcs(face2,row2,col2,shell)].w.row   = row;
    grid[idx_frcs(face2,row2,col2,shell)].w.col   = col;
    grid[idx_frcs(face2,row2,col2,shell)].w.rank   = mpi_rank;

    row++;   /*-- move -z direction. --*/
    row2++;  /*-- move -z direction. --*/
  }

  /*-- connect right/E of BACK to left/W of RIGHT --*/
  face  = BACK_FACE;
  row   = 0;                  /*-- row 0 => row +i  = -z direction. --*/
  col   = FACE_COLS - 1;      /*-- right. --*/

  face2 = RIGHT_FACE;
  row2  = 0;                  /*-- row 0 => row +i  = -z direction. --*/
  col2  = 0;                  /*-- left. --*/

  for (i = 0; i < FACE_ROWS; i++){
    grid[idx_frcs(face,row,col,shell)].e.face  = face2;
    grid[idx_frcs(face,row,col,shell)].e.row   = row2;
    grid[idx_frcs(face,row,col,shell)].e.col   = col2;
    grid[idx_frcs(face,row,col,shell)].e.rank   = mpi_rank;

    grid[idx_frcs(face2,row2,col2,shell)].w.face  = face;
    grid[idx_frcs(face2,row2,col2,shell)].w.row   = row;
    grid[idx_frcs(face2,row2,col2,shell)].w.col   = col;
    grid[idx_frcs(face2,row2,col2,shell)].w.rank   = mpi_rank;

    row++;   /*-- move -z direction. --*/
    row2++;  /*-- move -z direction. --*/
  }

  /*-- connect right/E of RIGHT to left/W of FRONT --*/
  face  = RIGHT_FACE;
  row   = 0;                  /*-- row 0 => row +i  = -z direction. --*/
  col   = FACE_COLS - 1;      /*-- right. --*/

  face2 = FRONT_FACE;
  row2  = 0;                  /*-- row 0 => row +i  = -z direction. --*/
  col2  = 0;                  /*-- left. --*/

  for (i = 0; i < FACE_ROWS; i++){
    grid[idx_frcs(face,row,col,shell)].e.face  = face2;
    grid[idx_frcs(face,row,col,shell)].e.row   = row2;
    grid[idx_frcs(face,row,col,shell)].e.col   = col2;
    grid[idx_frcs(face,row,col,shell)].e.rank   = mpi_rank;

    grid[idx_frcs(face2,row2,col2,shell)].w.face  = face;
    grid[idx_frcs(face2,row2,col2,shell)].w.row   = row;
    grid[idx_frcs(face2,row2,col2,shell)].w.col   = col;
    grid[idx_frcs(face2,row2,col2,shell)].w.rank   = mpi_rank;

    row++;   /*-- move -z direction. --*/
    row2++;  /*-- move -z direction. --*/
  }

} /*------ END  initNEWS( ) --------------------*/
/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/          void                                                /*---*/
/*--*/    checkNEWS( Index_t shell )                                 /*---*/
/*--*                                                                *---*/
/* Check to make sure initialization worked */
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{
  Index_t face, row, col;

  Neighbor_t neighbor;
  for (face=0; face<NUM_FACES; face++){
    for (row=0; row<FACE_ROWS; row++) {
      for (col=0; col<FACE_COLS; col++){

        neighbor = grid[idx_frcs(face,row,col,shell)].n;
        if ( ( ( neighbor.face >= NUM_FACES )
              ||   (neighbor.row >= FACE_ROWS))
            || (neighbor.col >= FACE_COLS) ) {
          printf("face=%d, row=%d, col=%d \n",(int) face, (int) row, (int) col);
          printf("    Problem Neighbor North \n");
          printf("       face=%d, row=%d, col=%d\n",
                 (int) grid[idx_frcs(face,row,col,shell)].n.face,
                 (int) grid[idx_frcs(face,row,col,shell)].n.row,
                 (int) grid[idx_frcs(face,row,col,shell)].n.col ); }

        neighbor = grid[idx_frcs(face,row,col,shell)].e;
        if ( ( ( neighbor.face >= NUM_FACES )
              ||   (neighbor.row >= FACE_ROWS))
            || (neighbor.col >= FACE_COLS) ) {
          printf("face=%d, row=%d, col=%d \n",(int) face, (int) row, (int) col);
          printf("    Problem Neighbor East \n");
          printf("       face=%d, row=%d, col=%d\n",
                 (int) grid[idx_frcs(face,row,col,shell)].e.face,
                 (int) grid[idx_frcs(face,row,col,shell)].e.row,
                 (int) grid[idx_frcs(face,row,col,shell)].e.col ); }

        neighbor = grid[idx_frcs(face,row,col,shell)].w;
        if ( ( ( neighbor.face >= NUM_FACES )
              ||   (neighbor.row >= FACE_ROWS))
            || (neighbor.col >= FACE_COLS) ) {
          printf("face=%d, row=%d, col=%d \n",(int) face, (int) row, (int) col);
          printf("    Problem Neighbor West \n");
          printf("       face=%d, row=%d, col=%d\n",
                 (int) grid[idx_frcs(face,row,col,shell)].w.face,
                 (int) grid[idx_frcs(face,row,col,shell)].w.row,
                 (int) grid[idx_frcs(face,row,col,shell)].w.col ); }

        neighbor = grid[idx_frcs(face,row,col,shell)].s;
        if ( ( ( neighbor.face >= NUM_FACES )
              ||   (neighbor.row >= FACE_ROWS))
            || (neighbor.col >= FACE_COLS) ) {
          printf("face=%d, row=%d, col=%d \n",(int) face, (int) row, (int) col);
          printf("    Problem Neighbor South \n");
          printf("       face=%d, row=%d, col=%d\n",
                 (int) grid[idx_frcs(face,row,col,shell)].s.face,
                 (int) grid[idx_frcs(face,row,col,shell)].s.row,
                 (int) grid[idx_frcs(face,row,col,shell)].s.col );
          printf("Grid faces not initialized correctly!\nCheck observer faces!\n");
          panic("Problem initializing EPREM grid!");
        }
      }
    }
  }
} /* END check NEWS */
/*------------------*/

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/          void                                                /*---*/
/*--*/    initCubeCoords( Index_t shell )                           /*---*/
/*--*                                                                *---*/
/*--*  Give this shell its initial x-y-z coordinates.                *---*/
/*--* The nodes on the inner shell are initially positioned          *---*/
/*--* on the unit cube. Later they will be projected to the sphere.  *---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{
  Index_t row, col, i;
  Coord_t cellWidth  = 1.0/(FACE_COLS);
  Coord_t cellHeight = 1.0/(FACE_ROWS);
  Coord_t x;
  Coord_t y;
  Coord_t z;

  /*-- TOP face ---------------------------------
   *--* Array orientation: viewing from outside the cube, origin at cube
   *--* center, top face in plane z = +0.5, front face in plane x = +0.5.
   *--* Viewing down on top face from +z direction, feet at -y, head at +y,
   *--* right arm in +x direction.
   *--* [row, col] index goes as,
   *--*                 [0, 0] => (-x, +y)
   *--*         [0, FACE_COLS] => (+x, +y)
   *--*         [FACE_ROWS, 0] => (-x, -y)
   *--* [FACE_ROWS, FACE_COLS] => (+x, -y)
   *----------------------------------------------------*/
  /*---- 2-d loop over every cell in TOP face -----------------*/
  /*---- Assign (x,y,z) coordinates for each cell. ----------- */
  /*---- Proceeds in row-major order from [0,0]. ------------- */
  x =  -0.5;                  /*-- -x cube extreme. --*/
  x += (0.5)*cellWidth;       /*-- move to center of cell. --*/
  y =   0.5;                  /*-- +y cube extreme. --*/
  y -= (0.5)*cellHeight;      /*-- move to center of cell. --*/
  row = col = 0;
  for ( i = 0; i < FACE_SIZE; i++) {
    if (col == FACE_COLS){   /*-- The end of the row? --*/
	    col = 0;               /*-- reset to -x extreme --*/
      x =  -0.5;             /*-- -x cube extreme. --*/
      x += (0.5)*cellWidth;  /*-- move to center of cell. --*/
	    row++;                 /*-- move one row in -y direction. --*/
	    y -= cellHeight;
    }
    grid[idx_frcs(TOP_FACE,row,col,shell)].r.x = x;
    grid[idx_frcs(TOP_FACE,row,col,shell)].r.y = y;
    grid[idx_frcs(TOP_FACE,row,col,shell)].r.z = +0.5;

    col++;                /*-- move in +x direction to next cell. --*/
    x += cellWidth;
  }

  /*-- BOTTOM face ---------------------------------
   *--* Face in plane z = -0.5, front face in plane x = +0.5.
   *--* Viewing up from -z direction, feet at -y, head at +y,
   *--* right arm in -x direction.
   *--*                 [0, 0] => (+x, +y)
   *--*         [0, FACE_COLS] => (-x, +y)
   *--*         [FACE_ROWS, 0] => (+x, -y)
   *--* [FACE_ROWS, FACE_COLS] => (-x, -y)
   *----------------------------------------------------*/
  /*---- 2-d loop over every cell in BOTTOM face -----------------*/
  /*---- Proceeds in row-major order from [0,0] => (+x, +y). -----*/
  x =  0.5;                   /*-- +x cube extreme. --*/
  x -= (0.5)*cellWidth;       /*-- move to center of cell. --*/
  y =   0.5;                  /*-- +y cube extreme. --*/
  y -= (0.5)*cellHeight;      /*-- move to center of cell. --*/
  row = col = 0;
  for ( i = 0; i < FACE_SIZE; i++) {
    if (col == FACE_COLS){  /*-- The end of the row? --*/
	    col = 0;              /*-- reset to +x extreme --*/
      x =  0.5;             /*-- +x cube extreme. --*/
      x -= (0.5)*cellWidth; /*-- move to center of cell. --*/
	    row++;                /*-- move one row in -y direction. --*/
	    y -= cellHeight;
    }
    grid[idx_frcs(BOTTOM_FACE,row,col,shell)].r.x = x;
    grid[idx_frcs(BOTTOM_FACE,row,col,shell)].r.y = y;
    grid[idx_frcs(BOTTOM_FACE,row,col,shell)].r.z = -0.5;

    col++;                /*-- move in -x direction to next cell. --*/
    x -= cellWidth;
  }

  /*-- FRONT face ---------------------------------
   *--* Face in plane x = +0.5.
   *--* Viewing up from +x direction, feet at -z, head at +z,
   *--* right arm in +y direction.
   *--*                 [0, 0] => (-y, +z)
   *--*         [0, FACE_COLS] => (+y, +z)
   *--*         [FACE_ROWS, 0] => (-y, -z)
   *--* [FACE_ROWS, FACE_COLS] => (+y, -z)
   *----------------------------------------------------*/
  /*---- 2-d loop over every cell in FRONT face -----------------*/
  /*---- Proceeds in row-major order from [0,0] => (-y, +z). -----*/
  y  = -0.5;                  /*-- -y cube extreme. --*/
  y += (0.5)*cellWidth;       /*-- move to center of cell. --*/
  z  =  0.5;                  /*-- +z cube extreme. --*/
  z -= (0.5)*cellHeight;      /*-- move to center of cell. --*/
  row = col = 0;
  for ( i = 0; i < FACE_SIZE; i++) {
    if (col == FACE_COLS){  /*-- The end of the row? --*/
	    col = 0;               /*-- reset to +y extreme --*/
      y   = -0.5;            /*-- -y cube extreme. --*/
      y  += (0.5)*cellWidth; /*-- move to center of cell. --*/
	    row++;                 /*-- move one row in -z direction. --*/
	    z -= cellHeight;
    }
    grid[idx_frcs(FRONT_FACE,row,col,shell)].r.x = +0.5;
    grid[idx_frcs(FRONT_FACE,row,col,shell)].r.y = y;
    grid[idx_frcs(FRONT_FACE,row,col,shell)].r.z = z;

    col++;                /*-- move in +y direction to next cell. --*/
    y += cellWidth;
  }

  /*-- BACK face ---------------------------------
   *--* Face in plane x = -0.5.
   *--* Viewing from -x direction, feet at -z, head at +z,
   *--* right arm in -y direction.
   *--*                 [0, 0] => (+y, +z)
   *--*         [0, FACE_COLS] => (-y, +z)
   *--*         [FACE_ROWS, 0] => (+y, -z)
   *--* [FACE_ROWS, FACE_COLS] => (-y, -z)
   *----------------------------------------------------*/
  /*---- 2-d loop over every cell in BACK face -----------------*/
  /*---- Proceeds in row-major order from [0,0] => (+y, +z). -----*/
  y =  0.5;                  /*-- +y cube extreme. --*/
  y -= (0.5)*cellWidth;       /*-- move to center of cell. --*/
  z =   0.5;                  /*-- +z cube extreme. --*/
  z -= (0.5)*cellHeight;      /*-- move to center of cell. --*/
  row = col = 0;
  for ( i = 0; i < FACE_SIZE; i++) {
    if (col == FACE_COLS){  /*-- The end of the row? --*/
	    col = 0;              /*-- reset to +y extreme --*/
      y =  0.5;             /*-- +y cube extreme. --*/
      y -= (0.5)*cellWidth; /*-- move to center of cell. --*/
	    row++;                /*-- move one row in -z direction. --*/
	    z -= cellHeight;
    }
    grid[idx_frcs(BACK_FACE,row,col,shell)].r.x = -0.5;
    grid[idx_frcs(BACK_FACE,row,col,shell)].r.y = y;
    grid[idx_frcs(BACK_FACE,row,col,shell)].r.z = z;

    col++;                /*-- move in -y direction to next cell. --*/
    y -= cellWidth;
  }

  /*-- RIGHT face ---------------------------------
   *--* Face in plane y = -0.5.
   *--* Viewing from -y direction, feet at -z, head at +z,
   *--* right arm in +x direction.
   *--*                 [0, 0] => (-x, +z)
   *--*         [0, FACE_COLS] => (+x, +z)
   *--*         [FACE_ROWS, 0] => (-x, -z)
   *--* [FACE_ROWS, FACE_COLS] => (+x, -z)
   *----------------------------------------------------*/
  /*---- 2-d loop over every cell in RIGHT face -----------------*/
  /*---- Proceeds in row-major order from [0,0] => (-x, +z). -----*/
  x  = -0.5;                  /*-- -x cube extreme. --*/
  x += (0.5)*cellWidth;       /*-- move to center of cell. --*/
  z  =  0.5;                  /*-- +z cube extreme. --*/
  z -= (0.5)*cellHeight;      /*-- move to center of cell. --*/
  row = col = 0;
  for ( i = 0; i < FACE_SIZE; i++) {
    if (col == FACE_COLS){  /*-- The end of the row? --*/
	    col = 0;              /*-- reset to -x extreme --*/
      x  = -0.5;            /*-- -x cube extreme. --*/
      x += (0.5)*cellWidth; /*-- move to center of cell. --*/
	    row++;                /*-- move one row in -z direction. --*/
	    z -= cellHeight;
    }
    grid[idx_frcs(RIGHT_FACE,row,col,shell)].r.x = x;
    grid[idx_frcs(RIGHT_FACE,row,col,shell)].r.y = -0.5;
    grid[idx_frcs(RIGHT_FACE,row,col,shell)].r.z = z;

    col++;                /*-- move in +x direction to next cell. --*/
    x += cellWidth;
  }

  /*-- LEFT face ---------------------------------
   *--* Face in plane y = +0.5.
   *--* Viewing from +y direction, feet at -z, head at +z,
   *--* right arm in -x direction.
   *--*                 [0, 0] => (+x, +z)
   *--*         [0, FACE_COLS] => (-x, +z)
   *--*         [FACE_ROWS, 0] => (+x, -z)
   *--* [FACE_ROWS, FACE_COLS] => (-x, -z)
   *----------------------------------------------------*/
  /*---- 2-d loop over every cell in LEFT face -----------------*/
  /*---- Proceeds in row-major order from [0,0] => (+x, +z). -----*/
  x =  0.5;                  /*-- +x cube extreme. --*/
  x -= (0.5)*cellWidth;       /*-- move to center of cell. --*/
  z =   0.5;                  /*-- +z cube extreme. --*/
  z -= (0.5)*cellHeight;      /*-- move to center of cell. --*/
  row = col = 0;
  for ( i = 0; i < FACE_SIZE; i++) {
    if (col == FACE_COLS){  /*-- The end of the row? --*/
	    col = 0;              /*-- reset to +x extreme --*/
      x =  0.5;             /*-- +x cube extreme. --*/
      x -= (0.5)*cellWidth; /*-- move to center of cell. --*/
	    row++;                /*-- move one row in -z direction. --*/
	    z -= cellHeight;
    }
    grid[idx_frcs(LEFT_FACE,row,col,shell)].r.x = x;
    grid[idx_frcs(LEFT_FACE,row,col,shell)].r.y = +0.5;
    grid[idx_frcs(LEFT_FACE,row,col,shell)].r.z = z;

    col++;                /*-- move in -x direction to next cell. --*/
    x -= cellWidth;
  }

} /*------ END  initCubeCoords( ) --------------------*/
/*----------------------------------------------------*/


/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/          void                                                /*---*/
/*--*/    initSphereCoords ( Index_t shell )                        /*---*/
/*--*                                                                *---*/
/*--*  Scale each node's position vector, r, to the unit sphere.     *---*/
/*--*  Only the inner shell is done, other shells get copied later.  *---*/
/*--*  Scaling is r = ( r / ||r|| ) for sphere of radius 1.          *---*/
/*--*  Assign rmag = ||r||. Also "azi" and "zen", the azimuth and    *---*/
/*--*  zenith/co-latitude of each node, get initialized.             *---*/
/*--*                                                                *---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{
  Index_t face, row, col;
  Vec_t r;
  Scalar_t  rmag;
  Radian_t  zen;
  Radian_t  azi;

  /*---- Loop over cube surface, rescale each node's r  --*/
  /*---- to unit by scaling (x,y,z) by ||r|| ( -> 1!).  --*/
  /*---- Assign azi = atan2( r.y, r.x ). We take        --*/
  /*---- azi to be in [-pi, pi].                        --*/
  /*---- Assign zen = acos(r.z / ||r||). We take        --*/
  /*---- zen to be in [0, pi].                          --*/
  for (face = 0; face < NUM_FACES; face++)
  {
    for (row  = 0; row  < FACE_ROWS; row++ )
    {
      for (col  = 0; col  < FACE_COLS; col++ )
      {

        r = grid[idx_frcs(face,row,col,shell)].r;
        rmag  = sqrt( (r.x * r.x) + (r.y * r.y) + (r.z * r.z) );
        r.x   = r.x / rmag;
        r.y   = r.y / rmag;
        r.z   = r.z / rmag;
        rmag  = 1.0;

        azi = atan2(r.y, r.x);
        zen = acos(r.z/rmag);

        grid[idx_frcs(face,row,col,shell)].rmag = 1.0;

        grid[idx_frcs(face,row,col,shell)].azi  = azi;
        grid[idx_frcs(face,row,col,shell)].zen  = zen;

        grid[idx_frcs(face,row,col,shell)].r.z =
        grid[idx_frcs(face,row,col,shell)].rmag * cos(zen) ;

        grid[idx_frcs(face,row,col,shell)].r.x =
        grid[idx_frcs(face,row,col,shell)].rmag * sin(zen)*cos(azi) ;

        grid[idx_frcs(face,row,col,shell)].r.y =
        grid[idx_frcs(face,row,col,shell)].rmag * sin(zen)*sin(azi) ;

      }
    }
  }

} /*------ END  initSphereCoords ( ) --------------------*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/


/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/          void                                                /*---*/
/*--*/    initCopyAll(void)                                         /*---*/
/*--*                                                                *---*/
/*--* Copy inner cube data to all shells.                            *---*/
/*--* Copy the NEWS links and position vector, r, from the inner     *---*/
/*--* shell to all the other shells. All shells then are at the      *---*/
/*--* unit sphere's surface.                                         *---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{
  Index_t face, row, col, shell;
  Radian_t zen, azi;
  Scalar_t rmag0;

  Neighbor_t neighborN;
  Neighbor_t neighborE;
  Neighbor_t neighborW;
  Neighbor_t neighborS;


  for (face = 0; face < NUM_FACES; face++ )
  {
    for (row = 0; row < FACE_ROWS; row++ )
    {
      for (col = 0; col < FACE_COLS; col++ )
      {

        neighborN = grid[idx_frcs(face,row,col,INNER_SHELL)].n;
        neighborE = grid[idx_frcs(face,row,col,INNER_SHELL)].e;
        neighborW = grid[idx_frcs(face,row,col,INNER_SHELL)].w;
        neighborS = grid[idx_frcs(face,row,col,INNER_SHELL)].s;
        azi       = grid[idx_frcs(face,row,col,INNER_SHELL)].azi;
        zen       = grid[idx_frcs(face,row,col,INNER_SHELL)].zen;
        rmag0     = grid[idx_frcs(face,row,col,INNER_SHELL)].rmag;

        for (shell = INNER_ACTIVE_SHELL; shell < LOCAL_NUM_SHELLS; shell++)
        {

          grid[idx_frcs(face,row,col,shell)].n = neighborN;
          grid[idx_frcs(face,row,col,shell)].e = neighborE;
          grid[idx_frcs(face,row,col,shell)].w = neighborW;
          grid[idx_frcs(face,row,col,shell)].s = neighborS;

          grid[idx_frcs(face,row,col,shell)].rmag = rmag0 + displGrid[mpi_rank] + shell;
        //  Original code was the following but I think it should have been LOCAL_NUM_SHELLS-1 in which case
        //  the above "+ mpi_rank +" should be removed.
        //  grid[idx_frcs(face,row,col,shell)].rmag = config.rScale + shell + LOCAL_NUM_SHELLS * mpi_rank;

          grid[idx_frcs(face,row,col,shell)].azi  = azi;
          grid[idx_frcs(face,row,col,shell)].zen  = zen;

          grid[idx_frcs(face,row,col,shell)].r.z =
          grid[idx_frcs(face,row,col,shell)].rmag * cos(zen);

          grid[idx_frcs(face,row,col,shell)].r.x =
          grid[idx_frcs(face,row,col,shell)].rmag * sin(zen)*cos(azi);

          grid[idx_frcs(face,row,col,shell)].r.y =
          grid[idx_frcs(face,row,col,shell)].rmag * sin(zen)*sin(azi);

          /*-- fix up shell numbers copied from inner shell --*/
          grid[idx_frcs(face,row,col,shell)].n.shell = shell;
          grid[idx_frcs(face,row,col,shell)].e.shell = shell;
          grid[idx_frcs(face,row,col,shell)].w.shell = shell;
          grid[idx_frcs(face,row,col,shell)].s.shell = shell;

        }
      }
    }
  }

}
/*------ END  initCopyAll ( ) --------------------*/
/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/          void                                                /*---*/
/*--*/    initStreamNeighbors(void)                                 /*---*/
/*--*                                                                *---*/
/*--* Give all shells stream neighbors. For the inner-most proc,     *---*/
/*--* its inner-most shell has no streamIn neighbor. (This is also   *---*/
/*--* true for the next shell, INNER_ACTIVE_SHELL.) A symmetric      *---*/
/*--* situation holds for the outter-most proc's OUTER_SHELL.       *---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{
  Index_t face, row, col, shell;

  for (shell = 0; shell < LOCAL_NUM_SHELLS; shell++){
    for (face  = 0; face  < NUM_FACES;  face++ ){
      for (row   = 0; row   < FACE_ROWS;  row++  ){
        for (col   = 0; col   < FACE_COLS;  col++  ){

          grid[idx_frcs(face,row,col,shell)].streamIn.shell  = shell - 1;
          grid[idx_frcs(face,row,col,shell)].streamOut.shell = shell + 1;
          grid[idx_frcs(face,row,col,shell)].streamIn.face   = face;
          grid[idx_frcs(face,row,col,shell)].streamOut.face  = face;
          grid[idx_frcs(face,row,col,shell)].streamIn.row    = row;
          grid[idx_frcs(face,row,col,shell)].streamOut.row   = row;
          grid[idx_frcs(face,row,col,shell)].streamIn.col    = col;
          grid[idx_frcs(face,row,col,shell)].streamOut.col   = col;
          grid[idx_frcs(face,row,col,shell)].streamIn.rank   = mpi_rank;
          grid[idx_frcs(face,row,col,shell)].streamOut.rank  = mpi_rank;

          /*-- Fix up rank/shell of inner/outer shells' stream neighbors. --*/
          /*-- Inner two shells usually link to the next inner proc.      --*/
          /*-- However, when we are on the inner-most proc itself, we     --*/
          /*-- fix this up by signaling there is no neighbor link in the  --*/
          /*-- streamIn direction, "NO_STREAM_NEIGHBOR". Likewise, a      --*/
          /*-- symmetric situation occurs at the outter-most proc for     --*/
          /*-- the streamOut link. The shell index also needs to be fixed --*/
          /*-- to accommodate for the inner-most shell being inactive.    --*/
          /*-- That is, the outter-most cube's streamOut neighbors are on --*/
          /*-- the second-inner-most cube on the next proc out (which     --*/
          /*-- holds the INNER_ACTIVE_SHELL on the next proc out).        --*/

          if ( (shell == INNER_SHELL) || (shell == INNER_ACTIVE_SHELL) ){

            grid[idx_frcs(face,row,col,shell)].streamIn.rank
            = PREV_PROC( mpi_rank );
            grid[idx_frcs(face,row,col,shell)].streamIn.shell
            = OUTER_SHELL;

            if (mpi_rank == INNER_PROC){
              grid[idx_frcs(face,row,col,shell)].streamIn.rank
              = NO_STREAM_NEIGHBOR;
            }
          }

          if (shell == OUTER_SHELL){
            grid[idx_frcs(face,row,col,shell)].streamOut.rank
            = NEXT_PROC( mpi_rank );
            grid[idx_frcs(face,row,col,shell)].streamOut.shell
            = INNER_ACTIVE_SHELL;

            if (mpi_rank == OUTER_PROC){
              grid[idx_frcs(face,row,col,shell)].streamOut.rank
              = NO_STREAM_NEIGHBOR;
            }
          }

        }}}}
} /*------ END  initStreamNeighbors ( ) ---------------------------------*/
/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/    void                                                      /*---*/
/*--*/    initInnerShellValues( void )                              /*---*/
/*--*                                                                *---*/
/*--* Initializing the values on the inner most shell on proc 0.     *---*/
/*--* This will be used as the template for when a new shell spawns. *---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{

  Index_t face, row, col;
  Index_t shell = 0;
  Index_t species, energy, mu;

  for (face = 0; face < NUM_FACES; face++) {
    for (row = 0; row < FACE_ROWS; row++) {
      for (col = 0; col < FACE_COLS; col++) {

        grid[idx_frcs(face,row,col,shell)].rOld.x = 0.0;
        grid[idx_frcs(face,row,col,shell)].rOld.y = 0.0;
        grid[idx_frcs(face,row,col,shell)].rOld.z = 0.0;
        grid[idx_frcs(face,row,col,shell)].rOlder.z = 0.0;
        grid[idx_frcs(face,row,col,shell)].rOlder.z = 0.0;
        grid[idx_frcs(face,row,col,shell)].rOlder.z = 0.0;
        grid[idx_frcs(face,row,col,shell)].ds = 0.0;
        grid[idx_frcs(face,row,col,shell)].dsOld = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdDensity = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdDensityOld = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdDivV = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdBr = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdBphi = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdBtheta = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdBmag = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdBmagPlus = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdBmagMinus = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdBvec.x = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdBvec.y = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdBvec.z = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdVvec.x = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdVvec.y = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdVvec.z = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdVsphOld.r = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdVsphOld.theta = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdVsphOld.phi = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdVr = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdVtheta = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdVphi = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdVmag = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdDlnB = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdDlnN = 0.0;
        grid[idx_frcs(face,row,col,shell)].mhdDuPar = 0.0;
        grid[idx_frcs(face,row,col,shell)].curlBoverB2.r = 0.0;
        grid[idx_frcs(face,row,col,shell)].curlBoverB2.theta = 0.0;
        grid[idx_frcs(face,row,col,shell)].curlBoverB2.phi = 0.0;

        for (species = 0; species < NUM_SPECIES; species++) {
          for (energy = 0; energy < NUM_ESTEPS; energy++) {
            for (mu = 0; mu < NUM_MUSTEPS; mu++) {

              eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)] = DBL_MIN;

            }
          }
        }

      }
    }
  }

} /*------ END  initStreamNeighbors ( ) ---------------------------------*/
/*-----------------------------------------------------------------------*/



