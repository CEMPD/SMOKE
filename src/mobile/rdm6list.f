
        SUBROUTINE RDM6LIST( LDEV )

C***********************************************************************
C  subroutine body starts at line 64
C
C  DESCRIPTION:
C       Reads the list of MOBILE6 input scenarios
C
C  PRECONDITIONS REQUIRED:
C       LDEV must be opened
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     10/01: Created by C. Seppanen
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
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
C***********************************************************************

C.........  MODULES for public variables
                
C...........   This module contains emission factor tables and related
        USE MODEMFAC, ONLY: M6LIST
                
        IMPLICIT NONE

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER        GETFLINE
        
        EXTERNAL  GETFLINE

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: LDEV     ! M6LIST file unit no.

C...........   Other local variables               

        INTEGER IOS                       ! I/O status        
        INTEGER NLINES                    ! number of lines

        CHARACTER*16 :: PROGNAME = 'RDM6LIST'   ! program name
        
C***********************************************************************
C   begin body of subroutine RDM6LIST

C.........  Get number of lines in file
        NLINES = GETFLINE( LDEV, 'M6LIST file' )
        
C.........  Allocate memory for storing the file
        ALLOCATE( M6LIST( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'M6LIST', PROGNAME )
        
C.........  Store lines of M6LIST file
        CALL RDLINES( LDEV, 'M6LIST file', NLINES, M6LIST )

        RETURN
        
        END SUBROUTINE RDM6LIST
        