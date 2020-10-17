/*-------------------------------------------------
 -- EMMREM: energeticParticlesInit.c
 --
 -- Utilities for energetic particle functionality.
 --
 --     initEnergeticParticles( void ) : set initial distributions.
 --
 -- _________________CHANGE HISTORY_____________
 -- _____________END_CHANGE HISTORY_____________
 --------------------------------------------------*/

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
#include "global.h"
#include "configuration.h"
#include "energeticParticlesInit.h"
#include "energeticParticlesBoundary.h"
#include "error.h"

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/             void                                         /*--*/
/*--*/     initEnergeticParticlesGrids( void )                  /*--*/
/*--                                                              --*/
/*--   Define the energy grid and EP attributes                   --*/
/*--                                                              --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  Index_t species, energy, mu;
  Scalar_t kin, p;
  Index_t face, row, col, shell;
  const double one   = 1.0;
  const double two   = 2.0;

  /*  T = Kin. Energy/( A * mp c^2 )                          -- */
  /*  p = momentum / ( A * mp * c )                           -- */
  /*  lnp(i) = lnpmin + dlnp*energy                           -- */
  for (species =0;              species < NUM_SPECIES; species++  ){

    /*    kin = config.charge[species]*config.eMin*MEV/(config.mass[species]*MP*C*C); */
    /* reset the config.eMin and config.eMax for Energy/nucleon */
    kin = config.eMin * MEV/(MP*C*C);
    lnpmin[species] = log(sqrt(kin*kin+two*kin ));

    kin = config.eMax * MEV/(MP*C*C);
    lnpmax[species] = log(sqrt(kin*kin+two*kin ));

    dlnp[species]   = (lnpmax[species]- lnpmin[species])/ (NUM_ESTEPS - 1.0);
  }

  for (species =0;              species < NUM_SPECIES; species++  ){
    for (energy = 0;            energy  < NUM_ESTEPS;  energy++   ){
      p = exp(lnpmin[species] + energy * dlnp[species]);
      pgrid[idx_se(species,energy)]    = p;
      egrid[idx_se(species,energy)]    = sqrt(one+p*p) - one;
      /*    printf("energy %d %e\n",energy, egrid[idx_se(species,energy)]*MP*C*C/MEV); */
      vgrid[idx_se(species,energy)]    = p/sqrt(one+p*p);
      rigidity[idx_se(species,energy)] = (config.mass[species]/config.charge[species]) * p * MZERO;
      rigidity[idx_se(species,energy)] = pow(rigidity[idx_se(species,energy)], config.rigidityPower);
      /*      printf("energy %d %e rigidity %e M0 %e\n",energy, egrid[idx_se(species,energy)]*MP*C*C/MEV, rigidity[idx_se(species,energy)],MZERO); */

    }
  }

  for (mu = 0; mu < NUM_MUSTEPS ; mu++) {
    mugrid[mu] = -one + ( mu * two / (NUM_MUSTEPS - one) );

  }

  /*-- 4d loop for every node in every shell. --*/
  for (face  = 0;              face  < NUM_FACES;  face++  ) {
    for (row   = 0;              row   < FACE_ROWS;  row++   ) {
      for (col   = 0;              col   < FACE_COLS;  col++   ) {
        for (shell = INNER_SHELL ;   shell < LOCAL_NUM_SHELLS; shell++ ) {
          for (species = 0;       species < NUM_SPECIES;  species++  ){
            for (energy  = 0;       energy  < NUM_ESTEPS;   energy++   ){
              for (mu      = 0;       mu      < NUM_MUSTEPS;  mu++   ){
                eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)] = 0.0;

              }}}/*-- endfor ---*/
        }}}}/*-- endfor() -*/

}/*-------- END initEnergeticParticlesGrids()  ---------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/             void                                         /*--*/
/*--*/     initEnergeticParticles( void )                       /*--*/
/*--                                                              --*/
/*--   Define EP attributes on grid                               --*/
/*--                                                              --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  Index_t species, energy, mu;
  Index_t face, row, col, shell;

  /* initializes all node points to the VS distribution */

  /*-- 4d loop for every node in every shell. --*/
  for (face  = 0;              face  < NUM_FACES;  face++  ) {
    for (row   = 0;              row   < FACE_ROWS;  row++   ) {
      for (col   = 0;              col   < FACE_COLS;  col++   ) {
        for (shell = INNER_SHELL ;   shell < LOCAL_NUM_SHELLS; shell++ ) {

          /*-- 3d loop for every species/energy/mu --*/
          for (species = 0;       species < NUM_SPECIES;  species++  ){
            for (energy  = 0;       energy  < NUM_ESTEPS;   energy++   ){
              for (mu      = 0;       mu      < NUM_MUSTEPS;  mu++   ){

                if ((config.useBoundaryFunction > 0) &&
                    (config.boundaryFunctionInitDomain > 0)) {

                  eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)] =
                  sepSeedFunction(egrid[idx_se(species,energy)],
                                  grid[idx_frcs(face,row,col,shell)].rmag);

                }
              }
            }
          }
        }
      }
    }
  }



}/*-------- END initEnergeticParticles()   -------------------------*/
/*------------------------------------------------------------------*/

