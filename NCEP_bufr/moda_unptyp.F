	MODULE MODA_UNPTYP

#ifndef NFILES_H
#define NFILES_H
	  USE MODV_NFILES
#endif

#ifdef DYNAMIC_ALLOCATION
	  INTEGER, ALLOCATABLE :: MSGUNP(:)
#else
	  INTEGER :: MSGUNP(NFILES)
#endif

	END MODULE
