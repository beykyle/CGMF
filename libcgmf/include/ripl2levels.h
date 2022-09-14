/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file ripl2levels.h

  \brief Read RIPL-2 discrete level files

*/

#ifndef __RIPL2LEVELS_H__
#define __RIPL2LEVELS_H__

#define LEVELDIRECTORY "levels/"
#define RIPLDATA "ripl3-levels.dat"
#define DISCRETEDATA "discreteLevels.dat" // used to be called 'cgmfDiscreteLevels-2015.dat'
//#define DISCRETEDATA "cgmf-discreteLevels-add10.dat" // used to be called 'cgmfDiscreteLevels-2015.dat'

#include "structur.h"

/**************************************/
/*      Selection for Nmax            */
/**************************************/
typedef enum {normal=0, extended=1, reassign=2, all=3} MaxLevelCtl;

/**************************************/
/*      ripl2levels.cpp               */
/**************************************/
#ifdef DEVUTIL
void	riplReadDiscreteLevelData ();
int     riplReadDiscreteLevels (ZAnumber *, Level *, MaxLevelCtl);
// void    riplReadDiscreteLevels (Nucleus *, MaxLevelCtl , int );
#endif

void	getRiplDiscreteLevels (Nucleus *, int);
void	getRiplDiscreteLevels (Nucleus *);
void	fixLevels(Nucleus *, bool, MaxLevelCtl);

void    riplDiscreteLevelsCleanup ();

void    readDiscreteLevelData ();
void	getDiscreteLevels (Nucleus *, int);
void	getDiscreteLevels (Nucleus *);

#endif //__RIPL2LEVELS_H__
