
        SUBROUTINE WRSRCGRPS( VNAME, JDATE, JTIME, INPUTFLAG, INPUTEMIS )

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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C**************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE M3UTILIO

        USE MODMERGE, ONLY: EMGGRD, NSRCGRP, NSGOUTPUT, GRPCNT,
     &                      IGRPNUM, SGINLNNAME, SRCGRPNAME,
     &                      PFLAG, PVNAME, PVSDATE, PVSTIME, ISRCGRP

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NGRID, NCOLS, NROWS,
     &                     GDTYP, GRDNM, P_ALP, P_BET, P_GAM,
     &                     XCENT, YCENT, XORIG, YORIG, XCELL, YCELL

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV, ONLY: NGROUP, NELEVGRPS, EMELEVGRP,
     &                     ELEVSTKGRP, ELEVSRCGRP, ELEVSTKCNT, SGFIREFLAG

        USE MODGRDLIB

        IMPLICIT NONE

C.........  INCLUDES:
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file descriptions
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C.........  EXTERNAL FUNCTIONS and their descriptions:
C       REAL     ENVREAL
C       INTEGER  FIND1, INDEX1

C        EXTERNAL ENVREAL, FIND1, INDEX1

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: VNAME   ! variable name to output
        INTEGER     , INTENT (IN) :: JDATE   ! Julian date to output (YYYYDDD)
        INTEGER     , INTENT (IN) :: JTIME   ! time to output (HHMMSS)
        LOGICAL     , INTENT (IN) :: INPUTFLAG      ! indicate if input emissions need to be used
        REAL        , INTENT (IN) :: INPUTEMIS( * ) ! input emissions (optional)

C...........   Local allocatable arrays
        INTEGER,      ALLOCATABLE :: ISTACK ( : ) ! group number
        INTEGER,      ALLOCATABLE :: STKCNT ( : ) ! num. srcs per group
        INTEGER,      ALLOCATABLE :: ROW    ( : ) ! row number
        INTEGER,      ALLOCATABLE :: COL    ( : ) ! column number
        INTEGER,      ALLOCATABLE :: LMAJOR ( : ) ! major source flag
        INTEGER,      ALLOCATABLE :: LPING  ( : ) ! PinG source flag
        INTEGER,      ALLOCATABLE :: INTDATA( : ) ! generic integer data
        REAL,         ALLOCATABLE :: XLOCA  ( : ) ! x-location at center of grid cell
        REAL,         ALLOCATABLE :: YLOCA  ( : ) ! y-location at center of grid cell
        REAL,         ALLOCATABLE :: LAT    ( : ) ! latitude of YLOCA
        REAL,         ALLOCATABLE :: LONG   ( : ) ! longitude of XLOCA
        REAL,         ALLOCATABLE :: STKDM  ( : ) ! inside stack diameter
        REAL,         ALLOCATABLE :: STKHT  ( : ) ! stack height
        REAL,         ALLOCATABLE :: STKTK  ( : ) ! stack exit temperature
        REAL,         ALLOCATABLE :: STKVE  ( : ) ! stack exit velocity
        REAL,         ALLOCATABLE :: STKFLW ( : ) ! stack exit flow rate
        REAL,         ALLOCATABLE :: ACRES  ( : ) ! acres burned for a fire
        REAL, SAVE,   ALLOCATABLE :: OUTEMIS( : ) ! output emissions
        REAL,         ALLOCATABLE :: REALDATA( : )! generic real data

C...........   Other local variables
        INTEGER          C, G, K, IDX  ! counters and indices
        INTEGER          IOS           ! i/o status
        INTEGER          ROWNUM        ! grid cell row
        INTEGER          COLNUM        ! grid cell column
        INTEGER          ELEVIDX       ! starting index for elevated sources

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
            ALLOCATE( LMAJOR( NSGOUTPUT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LMAJOR', PROGNAME )
            ALLOCATE( LPING( NSGOUTPUT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LPING', PROGNAME )
            ALLOCATE( XLOCA( NSGOUTPUT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'XLOCA', PROGNAME )
            ALLOCATE( YLOCA( NSGOUTPUT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'YLOCA', PROGNAME )
            ALLOCATE( LAT( NSGOUTPUT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LAT', PROGNAME )
            ALLOCATE( LONG( NSGOUTPUT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LONG', PROGNAME )
            ALLOCATE( STKDM( NSGOUTPUT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'STKDM', PROGNAME )
            ALLOCATE( STKHT( NSGOUTPUT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'STKHT', PROGNAME )
            ALLOCATE( STKTK( NSGOUTPUT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'STKTK', PROGNAME )
            ALLOCATE( STKVE( NSGOUTPUT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'STKVE', PROGNAME )
            ALLOCATE( STKFLW( NSGOUTPUT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'STKFLW', PROGNAME )
            ISTACK = 0  ! array
            STKCNT = 0
            ROW = 0
            COL = 0
            LMAJOR = 0
            LPING = 0
            XLOCA = BADVAL3
            YLOCA = BADVAL3
            LAT = BADVAL3
            LONG = BADVAL3

C.............  Set dummy stack parameter arrays based on environment settings
            STKDM  = ENVREAL( 'SRCGRP_STKDM',  'Stack diameter', 0.1, IOS )
            STKHT  = ENVREAL( 'SRCGRP_STKHT',  'Stack height', 0.1, IOS )
            STKTK  = ENVREAL( 'SRCGRP_STKTK',  'Stack exit temperature', 273.0, IOS )
            STKVE  = ENVREAL( 'SRCGRP_STKVE',  'Stack exit velocity', 0.1, IOS )
            STKFLW = ENVREAL( 'SRCGRP_STKFLW', 'Stack exit flow rate', 0.1, IOS )

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

C.............  Convert x and y grid cell locations to lat/lon
            LAT = YLOCA
            LONG = XLOCA
            CALL CONVRTLL( K, GDTYP, GRDNM,
     &                     P_ALP, P_BET, P_GAM,
     &                     XCENT, YCENT, LONG, LAT )

C.............  Append data for elevated source groups
            IF( PFLAG ) THEN
                ELEVIDX = K + 1

                DO G = 1, NELEVGRPS
                    IDX = ELEVIDX + G - 1
                    ISTACK( IDX ) = IGRPNUM( ELEVSRCGRP( G ) )
                END DO

                STKCNT( ELEVIDX:NSGOUTPUT ) = ELEVSTKCNT

                ALLOCATE( INTDATA( NGROUP ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INTDATA', PROGNAME )

                CALL INT_READ3( PVNAME, 'ROW', 1, PVSDATE, PVSTIME, INTDATA )
                DO G = 1, NELEVGRPS
                    IDX = ELEVIDX + G - 1
                    ROW( IDX ) = INTDATA( ELEVSTKGRP( G ) )
                END DO

                CALL INT_READ3( PVNAME, 'COL', 1, PVSDATE, PVSTIME, INTDATA )
                DO G = 1, NELEVGRPS
                    IDX = ELEVIDX + G - 1
                    COL( IDX ) = INTDATA( ELEVSTKGRP( G ) )
                END DO

                CALL INT_READ3( PVNAME, 'LMAJOR', 1, PVSDATE, PVSTIME, INTDATA )
                DO G = 1, NELEVGRPS
                    IDX = ELEVIDX + G - 1
                    LMAJOR( IDX ) = INTDATA( ELEVSTKGRP( G ) )
                END DO

                CALL INT_READ3( PVNAME, 'LPING', 1, PVSDATE, PVSTIME, INTDATA )
                DO G = 1, NELEVGRPS
                    IDX = ELEVIDX + G - 1
                    LPING( IDX ) = INTDATA( ELEVSTKGRP( G ) )
                END DO

                DEALLOCATE( INTDATA )

                ALLOCATE( REALDATA( NGROUP ), STAT=IOS )
                CALL CHECKMEM( IOS, 'REALDATA', PROGNAME )

                CALL REAL_READ3( PVNAME, 'XLOCA', 1, PVSDATE, PVSTIME, REALDATA )
                DO G = 1, NELEVGRPS
                    IDX = ELEVIDX + G - 1
                    XLOCA( IDX ) = REALDATA( ELEVSTKGRP( G ) )
                END DO

                CALL REAL_READ3( PVNAME, 'YLOCA', 1, PVSDATE, PVSTIME, REALDATA )
                DO G = 1, NELEVGRPS
                    IDX = ELEVIDX + G - 1
                    YLOCA( IDX ) = REALDATA( ELEVSTKGRP( G ) )
                END DO

                CALL REAL_READ3( PVNAME, 'STKDM', 1, PVSDATE, PVSTIME, REALDATA )
                DO G = 1, NELEVGRPS
                    IDX = ELEVIDX + G - 1
                    STKDM( IDX ) = REALDATA( ELEVSTKGRP( G ) )
                END DO

                CALL REAL_READ3( PVNAME, 'STKHT', 1, PVSDATE, PVSTIME, REALDATA )
                DO G = 1, NELEVGRPS
                    IDX = ELEVIDX + G - 1
                    STKHT( IDX ) = REALDATA( ELEVSTKGRP( G ) )
                END DO

                CALL REAL_READ3( PVNAME, 'STKTK', 1, PVSDATE, PVSTIME, REALDATA )
                DO G = 1, NELEVGRPS
                    IDX = ELEVIDX + G - 1
                    STKTK( IDX ) = REALDATA( ELEVSTKGRP( G ) )
                END DO

                CALL REAL_READ3( PVNAME, 'STKVE', 1, PVSDATE, PVSTIME, REALDATA )
                DO G = 1, NELEVGRPS
                    IDX = ELEVIDX + G - 1
                    STKVE( IDX ) = REALDATA( ELEVSTKGRP( G ) )
                END DO

                CALL REAL_READ3( PVNAME, 'STKFLW', 1, PVSDATE, PVSTIME, REALDATA )
                DO G = 1, NELEVGRPS
                    IDX = ELEVIDX + G - 1
                    STKFLW( IDX ) = REALDATA( ELEVSTKGRP( G ) )
                END DO

C.................  If lat/lon is in existing stack groups file, append it
                IF( READ3( PVNAME, 'LATITUDE', 1, PVSDATE, PVSTIME, REALDATA ) ) THEN

                    DO G = 1, NELEVGRPS
                        IDX = ELEVIDX + G - 1
                        LAT( IDX ) = REALDATA( ELEVSTKGRP( G ) )
                    END DO

C.....................  Assume longitude is available if latitude was
                    CALL REAL_READ3( PVNAME, 'LONGITUDE', 1, PVSDATE, PVSTIME, REALDATA )
                    DO G = 1, NELEVGRPS
                        IDX = ELEVIDX + G - 1
                        LONG( IDX ) = REALDATA( ELEVSTKGRP( G ) )
                    END DO

                ELSE

C.....................  Otherwise, convert x and y grid cell locations
                    LAT( ELEVIDX:NSGOUTPUT )  = YLOCA( ELEVIDX:NSGOUTPUT )
                    LONG( ELEVIDX:NSGOUTPUT ) = XLOCA( ELEVIDX:NSGOUTPUT )
                    CALL CONVRTLL( NSGOUTPUT-ELEVIDX+1, GDTYP, GRDNM,
     &                             P_ALP, P_BET, P_GAM, XCENT, YCENT,
     &                             LONG( ELEVIDX:NSGOUTPUT ),
     &                             LAT(  ELEVIDX:NSGOUTPUT ) )

                END IF

C.................  If acres burned is in existing stack groups file, add to output
                IF( SGFIREFLAG ) THEN

                    ALLOCATE( ACRES( NSGOUTPUT ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'ACRES', PROGNAME )
                    ACRES = 0.

                    CALL REAL_READ3( PVNAME, 'ACRESBURNED', 1, PVSDATE, PVSTIME, REALDATA )
                    DO G = 1, NELEVGRPS
                        IDX = ELEVIDX + G - 1
                        ACRES( IDX ) = REALDATA( ELEVSTKGRP( G ) )
                    END DO

                END IF

                DEALLOCATE( REALDATA )
            END IF

            CALL INT_WRITE3( SRCGRPNAME, 'IGROUP', JDATE, JTIME, ISTACK )
            CALL INT_WRITE3( SRCGRPNAME, 'GRPCNT', JDATE, JTIME, STKCNT )
            CALL INT_WRITE3( SRCGRPNAME, 'ROW',    JDATE, JTIME, ROW )
            CALL INT_WRITE3( SRCGRPNAME, 'COL',    JDATE, JTIME, COL )
            CALL INT_WRITE3( SRCGRPNAME, 'LMAJOR', JDATE, JTIME, LMAJOR )
            CALL INT_WRITE3( SRCGRPNAME, 'LPING',  JDATE, JTIME, LPING )

            CALL REAL_WRITE3( SRCGRPNAME, 'XLOCA', JDATE, JTIME, XLOCA )
            CALL REAL_WRITE3( SRCGRPNAME, 'YLOCA', JDATE, JTIME, YLOCA )
            CALL REAL_WRITE3( SRCGRPNAME, 'LATITUDE',  JDATE, JTIME, LAT )
            CALL REAL_WRITE3( SRCGRPNAME, 'LONGITUDE', JDATE, JTIME, LONG )
            CALL REAL_WRITE3( SRCGRPNAME, 'STKDM', JDATE, JTIME, STKDM )
            CALL REAL_WRITE3( SRCGRPNAME, 'STKHT', JDATE, JTIME, STKHT )
            CALL REAL_WRITE3( SRCGRPNAME, 'STKTK', JDATE, JTIME, STKTK )
            CALL REAL_WRITE3( SRCGRPNAME, 'STKVE', JDATE, JTIME, STKVE )
            CALL REAL_WRITE3( SRCGRPNAME, 'STKFLW', JDATE, JTIME, STKFLW )

            IF( SGFIREFLAG ) THEN

                CALL REAL_WRITE3( SRCGRPNAME, 'ACRESBURNED', JDATE, JTIME, ACRES )
                DEALLOCATE( ACRES )

            END IF

            DEALLOCATE( ISTACK, STKCNT, ROW, COL, LMAJOR, LPING )
            DEALLOCATE( XLOCA, YLOCA, LAT, LONG )
            DEALLOCATE( STKDM, STKHT, STKTK, STKVE, STKFLW )

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

                IF( INPUTFLAG ) THEN
                    OUTEMIS( K ) = OUTEMIS( K ) + INPUTEMIS( K )
                END IF
            END DO
        END DO

C.........  Append emissions for elevated source groups
        IF( PFLAG ) THEN
            DO G = 1, NELEVGRPS

                K = K + 1
                OUTEMIS( K ) = EMELEVGRP( G )

            END DO
        END IF

        IF( .NOT. WRITESET( SGINLNNAME, VNAME, ALLFILES,
     &                      JDATE, JTIME, OUTEMIS ) ) THEN

             MESG = 'Could not write "' // VNAME //
     &              '" to file "' // SGINLNNAME // '"'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

        END IF

        RETURN

C*****************  INTERNAL SUBPROGRAMS  ******************************

        CONTAINS

C.............  This internal subprogram reads real data from an
C               I/O API file, and aborts if not successful.
            SUBROUTINE REAL_READ3( FILNAM, VARNAM, LAYER,
     &                             RDATE, RTIME, REALBUF )

C.............  Subprogram arguments
            CHARACTER(*) FILNAM       ! logical file name
            CHARACTER(*) VARNAM       ! variable name
            INTEGER      LAYER        ! layer number
            INTEGER      RDATE        ! read Julian date
            INTEGER      RTIME        ! read time
            REAL         REALBUF(*)   ! real data buffer

C.............  Local variables
            INTEGER L1, L2

C----------------------------------------------------------------------

            IF ( GET_VTYPE( FILNAM, VARNAM, RDATE, RTIME ) .NE. M3REAL ) THEN
                    MESG = 'Type for variable "' // TRIM( VARNAM ) //
     &                     '" in "'// TRIM( FILNAM ) // '" not M3REAL'
                    CALL M3EXIT( PROGNAME, RDATE, RTIME, MESG, 2 )
            ELSE IF ( .NOT. READ3( FILNAM, VARNAM, LAYER,
     &                             RDATE, RTIME, REALBUF ) ) THEN

                L1 = LEN_TRIM( VARNAM )
                L2 = LEN_TRIM( FILNAM )
                MESG = 'Could not read "' // TRIM( VARNAM ) //
     &                 '" from file "' // FILNAM( 1:L2 ) // '"'
                CALL M3EXIT( PROGNAME, RDATE, RTIME, MESG, 2 )

            END IF

            END SUBROUTINE REAL_READ3

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram reads integer data from an
C               I/O API file, and aborts if not successful.
            SUBROUTINE INT_READ3( FILNAM, VARNAM, LAYER,
     &                            RDATE, RTIME, INTBUF )

C.............  Subprogram arguments
            CHARACTER(*) FILNAM       ! logical file name
            CHARACTER(*) VARNAM       ! variable name
            INTEGER      LAYER        ! layer number
            INTEGER      RDATE        ! read Julian date
            INTEGER      RTIME        ! read time
            INTEGER      INTBUF(*)    ! integer data buffer

C----------------------------------------------------------------------

            IF ( GET_VTYPE( FILNAM, VARNAM, RDATE, RTIME ) .NE. M3INT ) THEN

                    MESG = 'Type for variable "' // TRIM( VARNAM ) //
     &                     '" in "'              // TRIM( FILNAM ) //
     &                     '" not M3INT'
                    CALL M3EXIT( PROGNAME, RDATE, RTIME, MESG, 2 )

            ELSE IF ( .NOT. READ3( FILNAM, VARNAM, LAYER,
     &                             RDATE, RTIME, INTBUF ) ) THEN

                MESG = 'Could not read "' // TRIM( VARNAM ) //
     &                 '" from file "'    // TRIM( FILNAM ) // '"'
                CALL M3EXIT( PROGNAME, RDATE, RTIME, MESG, 2 )

            END IF

            END SUBROUTINE INT_READ3

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram writes real data to an
C               I/O API file, and aborts if not successful.
            SUBROUTINE REAL_WRITE3( FILNAM, VARNAM,
     &                              WDATE, WTIME, REALBUF )

C.............  Subprogram arguments
            CHARACTER(*) FILNAM       ! logical file name
            CHARACTER(*) VARNAM       ! variable name
            INTEGER      WDATE        ! write Julian date
            INTEGER      WTIME        ! write time
            REAL         REALBUF(*)   ! real data buffer

C----------------------------------------------------------------------

            IF ( GET_VTYPE( FILNAM, VARNAM, WDATE, WTIME ) .NE. M3REAL ) THEN

                    MESG = 'Type for variable "' // TRIM( VARNAM ) //
     &                     '" in "'              // TRIM( FILNAM ) //
     &                     '" not M3REAL'
                    CALL M3EXIT( PROGNAME, WDATE, WTIME, MESG, 2 )

            ELSE IF ( .NOT. WRITE3( FILNAM, VARNAM,
     &                         WDATE, WTIME, REALBUF ) ) THEN

                MESG = 'Could not write "' // TRIM( VARNAM ) //
     &                 '" to file "'       // TRIM( FILNAM ) // '"'
                CALL M3EXIT( PROGNAME, WDATE, WTIME, MESG, 2 )

            END IF

            END SUBROUTINE REAL_WRITE3

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram writes integer data to an
C               I/O API file, and aborts if not successful.
            SUBROUTINE INT_WRITE3( FILNAM, VARNAM,
     &                             WDATE, WTIME, INTBUF )

C.............  Subprogram arguments
            CHARACTER(*) FILNAM       ! logical file name
            CHARACTER(*) VARNAM       ! variable name
            INTEGER      WDATE        ! write Julian date
            INTEGER      WTIME        ! write time
            INTEGER      INTBUF(*)    ! integer data buffer

C----------------------------------------------------------------------

            IF ( GET_VTYPE( FILNAM, VARNAM, WDATE, WTIME ) .NE. M3INT ) THEN

                    MESG = 'Type for variable "' // TRIM( VARNAM ) //
     &                     '" in "'              // TRIM( FILNAM ) //
     &                     '" not M3INT'
                    CALL M3EXIT( PROGNAME, WDATE, WTIME, MESG, 2 )

            ELSE IF ( .NOT. WRITE3( FILNAM, VARNAM,
     &                              WDATE, WTIME, INTBUF ) ) THEN

                MESG = 'Could not write "' // TRIM( VARNAM ) //
     &                 '" to file "'       // TRIM( FILNAM ) // '"'
                CALL M3EXIT( PROGNAME, WDATE, WTIME, MESG, 2 )

            END IF

            END SUBROUTINE INT_WRITE3

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram writes integer data to an
C               I/O API file, and aborts if not successful.
            INTEGER FUNCTION GET_VTYPE( FILNAM, VARNAM, JDATE, JTIME )

                CHARACTER(LEN=*), INTENT(IN   ) :: FILNAM, VARNAM
                INTEGER         , INTENT(IN   ) :: JDATE, JTIME

                INTEGER     K

                IF ( .NOT.DESC3( FILNAM ) ) THEN
                    MESG = 'File "' // TRIM( FILNAM ) // '" not available'
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF

                K = INDEX1( VARNAM, NVARS3D, VNAME3D )
                IF ( K .LE. 0 ) THEN
                    MESG = 'Variable "'          // TRIM( VARNAM ) //
     &                     '" not available in "'// TRIM( FILNAM ) // '"'
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF

                GET_VTYPE = VTYPE3D( K )
                RETURN

            END FUNCTION GET_VTYPE

        END SUBROUTINE WRSRCGRPS
