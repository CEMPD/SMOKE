
        SUBROUTINE OPENMRGOUT( NGRP )

C***********************************************************************
C  subroutine OPENMRGOUT body starts at line 83
C
C  DESCRIPTION:
C      The purpose of this subroutine is to open all of the necessary
C      output files for the merge routine (both I/O API and ASCII files)
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 2/99 by M. Houyoux
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV

        IMPLICIT NONE

C.........  INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C.........  EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*16    MULTUNIT
        INTEGER         PROMPTFFILE  
        CHARACTER*16    PROMPTMFILE  

        EXTERNAL  MULTUNIT, PROMPTFFILE, PROMPTMFILE

C...........  SUBROUTINE ARGUMENTS
       INTEGER, INTENT (IN) :: NGRP     ! Actual number of groups

C.........  Other local variables

        INTEGER         I, J, K, L, L1, L2, LJ, N, V

        CHARACTER*50    RTYPNAM      ! name for type of state/county report file
        CHARACTER*50    UNIT         ! tmp units buffer
        CHARACTER*300   BUFFER       ! tmp buffer
        CHARACTER*300   MESG         ! message buffer

        CHARACTER*16 :: PROGNAME = 'OPENMRGOUT' ! program name

C***********************************************************************
C   begin body of subroutine OPENMRGOUT

C.........  Set default output file names
        CALL MRGONAMS

C.........  Set up header for I/O API output files
        FTYPE3D = GRDDED3
        P_ALP3D = DBLE( P_ALP )
        P_BET3D = DBLE( P_BET )
        P_GAM3D = DBLE( P_GAM )
        XCENT3D = DBLE( XCENT )
        YCENT3D = DBLE( YCENT )
        XORIG3D = DBLE( XORIG )
        YORIG3D = DBLE( YORIG )
        XCELL3D = DBLE( XCELL )
        YCELL3D = DBLE( YCELL )
        SDATE3D = SDATE
        STIME3D = STIME
        TSTEP3D = TSTEP
        NCOLS3D = NCOLS
        NROWS3D = NROWS
        NLAYS3D = EMLAYS
        NTHIK3D = 1
        GDTYP3D = GDTYP
        VGTYP3D = VGTYP
        VGTOP3D = VGTOP
        GDNAM3D = GRDNM

C.........  Set constants number and values for variables
C.........  Do this regardless of whether we have outputs or not
        IF( SFLAG ) THEN
            NVARS3D = NMSPC

            K = 0
            LJ = -1
            DO N = 1, NGRP
                DO V = 1, VGRPCNT( N )

                    I = SIINDEX( V,N )
                    J = SPINDEX( V,N )
                    IF( J .EQ. LJ ) CYCLE    ! Do not repeat species

                    K = K + 1
                    VNAME3D( K ) = EMNAM  ( J )
                    UNITS3D( K ) = GRDUNIT( I )
                    VDESC3D( K ) = 'Model species ' // EMNAM( J )
                    VTYPE3D( K ) = M3REAL

                    LJ = J

                END DO
            END DO

        ELSE
            NVARS3D = NIPPA

            DO V = 1, NVARS3D
                VNAME3D( V ) = EANAM  ( V )
                UNITS3D( V ) = GRDUNIT( V )
                VDESC3D( V ) = 'Pollutant or activity ' // EANAM( V )
                VTYPE3D( V ) = M3REAL
            END DO

        ENDIF

C.........  Set up and open I/O API output file
        IF( LGRDOUT ) THEN
          
C.............  Prompt for and gridded open file(s)
            IF( AFLAG ) THEN
                AONAME = PROMPTMFILE(  
     &            'Enter name for AREA-SOURCE GRIDDED OUTPUT file',
     &            FSUNKN3, AONAME, PROGNAME )
            END IF 

            IF( MFLAG ) THEN
                MONAME = PROMPTMFILE(  
     &            'Enter name for MOBILE-SOURCE GRIDDED OUTPUT file',
     &            FSUNKN3, MONAME, PROGNAME )
            END IF 

            IF( PFLAG ) THEN
                PONAME = PROMPTMFILE(  
     &            'Enter name for POINT-SOURCE GRIDDED OUTPUT file',
     &            FSUNKN3, PONAME, PROGNAME )
            END IF 

            IF( XFLAG ) THEN
                TONAME = PROMPTMFILE(  
     &            'Enter name for MULTI-SOURCE GRIDDED OUTPUT file',
     &            FSUNKN3, TONAME, PROGNAME )
            END IF 

        END IF  ! End of gridded output

C.........  Open plume-in-grid output
        IF( PINGFLAG ) THEN

C.............  Override gridded file settings
            NCOLS3D = 1
            NROWS3D = NGROUP
            NLAYS3D = 1
            GDTYP3D = GDTYP
            VGTYP3D = BADVAL3
            VGTOP3D = BADVAL3

            FDESC3D = ' '   ! array
            
            PINGNAME = PROMPTMFILE( 
     &                     'Enter name for PING EMISSIONS file',
     &                     FSUNKN3, PINGNAME, PROGNAME )
        END IF

C.........  Open report file(s)

        IF( LREPSTA .OR. LREPCNY ) THEN

            IF( LREPSTA .AND. LREPCNY ) THEN
                RTYPNAM = 'STATE AND COUNTY'
            ELSE IF ( LREPSTA ) THEN
                RTYPNAM = 'STATE'
            ELSE IF ( LREPCNY ) THEN
                RTYPNAM = 'COUNTY'
            END IF

            L = LEN_TRIM( RTYPNAM )

            IF( AFLAG ) THEN
                ARDEV  = PROMPTFFILE(
     &                  'Enter name for AREA-SOURCE ' // 
     &                  RTYPNAM( 1:L ) // ' REPORT', 
     &                  .FALSE., .TRUE., AREPNAME, PROGNAME )
            END IF 

            IF( BFLAG ) THEN
                BRDEV  = PROMPTFFILE(
     &                  'Enter name for BIOGENIC ' // 
     &                  RTYPNAM( 1:L ) // ' REPORT', 
     &                  .FALSE., .TRUE., BREPNAME, PROGNAME )
            END IF 

            IF( MFLAG ) THEN
                MRDEV  = PROMPTFFILE(
     &                  'Enter name for MOBILE-SOURCE ' // 
     &                  RTYPNAM( 1:L ) // ' REPORT', 
     &                  .FALSE., .TRUE., MREPNAME, PROGNAME )
            END IF 

            IF( PFLAG ) THEN
                PRDEV  = PROMPTFFILE(
     &                  'Enter name for POINT-SOURCE ' // 
     &                  RTYPNAM( 1:L ) // ' REPORT', 
     &                  .FALSE., .TRUE., PREPNAME, PROGNAME )
            END IF 

            IF( XFLAG ) THEN
                TRDEV  = PROMPTFFILE(
     &                  'Enter name for TOTAL ' // 
     &                  RTYPNAM( 1:L ) // ' REPORT', 
     &                  .FALSE., .TRUE., TREPNAME, PROGNAME )
            END IF 

        END IF  ! End of state and/or county output

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94050   FORMAT( A, 1X, I2.2, A, 1X, A, 1X, I6.6, 1X,
     &          A, 1X, I3.3, 1X, A, 1X, I3.3, 1X, A   )

        END SUBROUTINE OPENMRGOUT
