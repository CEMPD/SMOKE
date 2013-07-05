
        SUBROUTINE WRSRCGRPS( VNAME, JDATE, JTIME )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C     This subroutine writes the emissions for each source
C     group and grid cell into the CMAQ inline emissions file.
C     The first time the subroutine is called, it also writes
C     the group metadata (group number, number of sources, dummy
C     stack parameters, etc.) into the stack groups file.
C
C  PRECONDITIONS REQUIRED:
C     Stack groups and inline emissions files opened for output
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 7/2013 by C. Seppanen
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2013, Environmental Modeling for Policy Development
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
C**************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE, ONLY: EMGGRD, NSRCGRP, NSGOUTPUT, GRPCNT,
     &                      IGRPNUM, INLINENAME, SRCGRPNAME

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NGRID, NCOLS, NROWS, 
     &                     XORIG, YORIG, XCELL, YCELL
        
        IMPLICIT NONE

C.........  INCLUDES:
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: VNAME   ! variable name to output
        INTEGER     , INTENT (IN) :: JDATE   ! Julian date to output (YYYYDDD)
        INTEGER     , INTENT (IN) :: JTIME   ! time to output (HHMMSS)

C...........   Local allocatable arrays
        INTEGER,      ALLOCATABLE :: ISTACK ( : ) ! group number
        INTEGER,      ALLOCATABLE :: STKCNT ( : ) ! num. srcs per group
        INTEGER,      ALLOCATABLE :: ROW    ( : ) ! row number
        INTEGER,      ALLOCATABLE :: COL    ( : ) ! column number
        REAL,         ALLOCATABLE :: XLOCA  ( : ) ! x-location at center of grid cell
        REAL,         ALLOCATABLE :: YLOCA  ( : ) ! y-location at center of grid cell
        REAL,         ALLOCATABLE :: OUTEMIS( : ) ! output emissions

C...........   Other local variables
        INTEGER          C, G, K       ! counters and indices
        INTEGER          IOS           ! i/o status
        INTEGER          ROWNUM        ! grid cell row
        INTEGER          COLNUM        ! grid cell column
        
        REAL             XLOCACELL     ! x-location for grid cell
        REAL             YLOCACELL     ! y-location for grid cell

        LOGICAL, SAVE :: FIRSTTIME = .TRUE. ! true: first time routine called

        CHARACTER(300)   MESG         ! message buffer

        CHARACTER(16) :: PROGNAME = 'WRSRCGRPS' ! program name

***********************************************************************
C   begin body of subroutine WRSRCGRPS

        IF( FIRSTTIME ) THEN
        
C.............  Output stack groups file
            ALLOCATE( ISTACK( NSGOUTPUT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ISTACK', PROGNAME )
            ALLOCATE( STKCNT( NSGOUTPUT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'STKCNT', PROGNAME )
            ALLOCATE( ROW( NSGOUTPUT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ROW', PROGNAME )
            ALLOCATE( COL( NSGOUTPUT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'COL', PROGNAME )
            ALLOCATE( XLOCA( NSGOUTPUT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'XLOCA', PROGNAME )
            ALLOCATE( YLOCA( NSGOUTPUT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'YLOCA', PROGNAME )
            ISTACK = 0  ! array
            STKCNT = 0
            ROW = 0
            COL = 0
            XLOCA = BADVAL3
            YLOCA = BADVAL3
            
            K = 0
            DO C = 1, NGRID

C.................  Determine row and column for current grid cell
                ROWNUM = C / NCOLS   ! integer math
                IF( MOD( C, NCOLS ) .GT. 0 ) ROWNUM = ROWNUM + 1
                COLNUM = C - ( ROWNUM-1 ) * NCOLS

C.................  Calculate x and y-position at center of grid cell
                XLOCACELL = XORIG + ( COLNUM-1 ) * XCELL + 0.5 * XCELL
                YLOCACELL = YORIG + ( ROWNUM-1 ) * YCELL + 0.5 * YCELL

                DO G = 1, NSRCGRP

C.....................  Skip missing values
                    IF( GRPCNT( C, G ) == 0 ) CYCLE
                    
                    K = K + 1
                    ISTACK( K ) = IGRPNUM( G )
                    STKCNT( K ) = GRPCNT( C, G )
                    ROW   ( K ) = ROWNUM
                    COL   ( K ) = COLNUM
                    XLOCA ( K ) = XLOCACELL
                    YLOCA ( K ) = YLOCACELL
                END DO
            END DO

            MESG = 'Error writing to output file "' // SRCGRPNAME // '"'
        
            IF( .NOT. WRITE3( SRCGRPNAME, 'ISTACK', JDATE, JTIME, ISTACK )) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        
            IF( .NOT. WRITE3( SRCGRPNAME, 'STKCNT', JDATE, JTIME, STKCNT )) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        
            IF( .NOT. WRITE3( SRCGRPNAME, 'ROW', JDATE, JTIME, ROW )) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        
            IF( .NOT. WRITE3( SRCGRPNAME, 'COL', JDATE, JTIME, COL )) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        
            IF( .NOT. WRITE3( SRCGRPNAME, 'XLOCA', JDATE, JTIME, XLOCA )) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        
            IF( .NOT. WRITE3( SRCGRPNAME, 'YLOCA', JDATE, JTIME, YLOCA )) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            
            DEALLOCATE( ISTACK, STKCNT, ROW, COL, XLOCA, YLOCA )

C.............  Allocate space for output emissions
            ALLOCATE( OUTEMIS( NSGOUTPUT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'OUTEMIS', PROGNAME )

            FIRSTTIME = .FALSE.
        END IF
        
        OUTEMIS = 0.  ! array
        
        K = 0
        DO C = 1, NGRID
            DO G = 1, NSRCGRP

C.................  Skip missing values
                IF( GRPCNT( C, G ) == 0 ) CYCLE

                K = K + 1
                OUTEMIS( K ) = EMGGRD( C, G )
            END DO
        END DO

        IF( .NOT. WRITESET( INLINENAME, VNAME, ALLFILES,
     &                      JDATE, JTIME, OUTEMIS ) ) THEN

             MESG = 'Could not write "' // VNAME //
     &              '" to file "' // INLINENAME // '"'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

        END IF
        
        RETURN

        END SUBROUTINE WRSRCGRPS
