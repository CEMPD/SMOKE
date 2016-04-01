
        SUBROUTINE BLDPROCIDX

C***********************************************************************
C  subroutine BLDPROCIDX body starts at line
C
C  DESCRIPTION:
C      This subroutine builds a mapping from the list of MOVES emission
C      processes to the pollutant-species combinations.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 4/10 by C. Seppanen
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
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: RPDFLAG, RPVFLAG, EMPOLIDX

C.........  This module contains the major data structure and control flags
        USE MODMERGE, ONLY: NSMATV, TSVDESC

        IMPLICIT NONE

C.........  INCLUDES:
        INCLUDE 'EMCNST3.EXT'    !  emissions constant parameters
        INCLUDE 'MVSCNST3.EXT'   !  MOVES constants

C.........  EXTERNAL FUNCTIONS and their descriptions:
        INTEGER   INDEX1
        
        EXTERNAL  INDEX1

C.........  Other local variables
        INTEGER   J, K, L1, V  ! counters and indexes
        INTEGER   IOS       ! error status

        CHARACTER(IOVLEN3) :: CPROC ! tmp process buffer
        CHARACTER(PLSLEN3) :: SVBUF ! tmp speciation name buffer
        CHARACTER(300) :: MESG      ! message buffer

        CHARACTER(16) :: PROGNAME = 'BLDPROCIDX' ! program name

C***********************************************************************
C   begin body of subroutine BLDPROCIDX

        ALLOCATE( EMPOLIDX( NSMATV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMPOLIDX', PROGNAME )
        EMPOLIDX = 0   ! array

C.........  Loop through pollutant-species combos
        DO V = 1, NSMATV
        
            SVBUF = TSVDESC( V )
            L1 = INDEX( SVBUF, ETJOIN )
            CPROC = SVBUF( 1:L1-1 )

C.............  Find emission process in master MOVES list
            IF( RPDFLAG ) THEN
                J = INDEX1( CPROC, MXMVSDPROCS, MVSDPROCS )
            ELSE IF( RPVFLAG ) THEN
                J = INDEX1( CPROC, MXMVSVPROCS, MVSVPROCS )
            ELSE
                J = INDEX1( CPROC, MXMVSPPROCS, MVSPPROCS )
            END IF

            IF( J .LE. 0 ) THEN

C.................  Check if process is handled by other modes
                IF( RPDFLAG ) THEN
                    J = INDEX1( CPROC, MXMVSVPROCS, MVSVPROCS )
                    K = INDEX1( CPROC, MXMVSPPROCS, MVSPPROCS )
                ELSE IF( RPVFLAG ) THEN
                    J = INDEX1( CPROC, MXMVSDPROCS, MVSDPROCS )
                    K = INDEX1( CPROC, MXMVSPPROCS, MVSPPROCS )
                ELSE
                    J = INDEX1( CPROC, MXMVSDPROCS, MVSDPROCS )
                    K = INDEX1( CPROC, MXMVSVPROCS, MVSVPROCS )
                END IF
                
                IF( J .LE. 0 .AND. K .LE. 0 ) THEN
                    MESG = 'ERROR: Requested emission process ' // 
     &                TRIM( CPROC ) // ' is not output by MOVES.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
            END IF

        END DO

        END SUBROUTINE BLDPROCIDX
