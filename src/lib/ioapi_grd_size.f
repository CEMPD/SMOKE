
       INTEGER FUNCTION IOAPI_GRD_SIZE( NCOLS, NROWS, NLAYS, 
     &                                  NVARS, NSTEPS )

C***********************************************************************
C  function body starts at line
C
C  DESCRIPTION:
C      This function returns the approximate size in megabytes of a 
C      gridded I/O API file based on the number of columns, rows, 
C      layers, variables, and time steps.
C
C  PRECONDITIONS REQUIRED:

C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created by C. Seppanen 3/03
C
C**************************************************************************
C
C Project Title: EDSS Tools Library
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

       IMPLICIT NONE
       
C........  Function arguments
       INTEGER, INTENT (IN) :: NCOLS   ! number of columns
       INTEGER, INTENT (IN) :: NROWS   ! number of rows
       INTEGER, INTENT (IN) :: NLAYS   ! number of layers
       INTEGER, INTENT (IN) :: NVARS   ! number of variables
       INTEGER, INTENT (IN) :: NSTEPS  ! number of time steps

C........  Other local variables
       INTEGER              NCELLS     ! number of grid cells
       INTEGER              HDRSIZE    ! size of header in bytes
       INTEGER              RECSIZE    ! size of single record in bytes

       CHARACTER(16) :: PROGNAME = 'IOAPI_GRD_SIZE'  ! program name

C***********************************************************************
C   begin body of function IOAPI_GRD_SIZE

C........  Calculate number of grid cells
       NCELLS = NCOLS * NROWS
       
C........  Calculate size of header
       HDRSIZE = ( 9860 + 116 * NVARS + 4 * NLAYS ) / 1000000
       
C........  Calculate size of individual records
       RECSIZE = ( 8 * NVARS + 4 * NLAYS * NVARS * NCELLS ) / 1000000

C........  Calculate total number of bytes in file       
       IOAPI_GRD_SIZE = HDRSIZE + RECSIZE * NSTEPS

       RETURN
      
       END FUNCTION IOAPI_GRD_SIZE
