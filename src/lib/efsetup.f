
        SUBROUTINE EFSETUP( FNAME, MODELNAM, TYPNAM, MXVAR, NVAR, 
     &                      VNAMES, VUNITS, VDESCS )
   
C***********************************************************************
C  subroutine EFSETUP body starts at line < >
C
C  DESCRIPTION:
C      Get the names of the emission factors given the emission factor
C      model name and locally obtained environment variable settings
C      NOTE - the NVAR argument is not used, but it is there for a future,
C      more flexible, version of the routine.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION HISTORY:
C
C***************************************************************************
C 
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C 
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C 
C See file COPYRIGHT for conditions of use.
C 
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C****************************************************************************

C...........   MODULES for public variables
C.........  This module contains emission factor tables and related
        USE MODEMFAC

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'M5CNST3.EXT'   !  Mobile5a/b constants
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data stru

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER     INDEX1
        EXTERNAL    INDEX1

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: FNAME           ! logical file or 'NONE'
        CHARACTER(*), INTENT (IN) :: MODELNAM        ! name of EF model
        CHARACTER(*), INTENT (IN) :: TYPNAM          ! factor type from EF model
        INTEGER     , INTENT (IN) :: MXVAR           ! max no of variables
        INTEGER     , INTENT(OUT) :: NVAR            ! actual no of variables
        CHARACTER(*), INTENT(OUT) :: VNAMES( MXVAR ) ! variable names
        CHARACTER(*), INTENT(OUT) :: VUNITS( MXVAR ) ! variable units
        CHARACTER(*), INTENT(OUT) :: VDESCS( MXVAR ) ! variable descriptions

C...........   Local variables
        INTEGER         I, K, L, L2 ! counters and indices

        INTEGER         IOS     ! status from retrieving E.V.s
        INTEGER         LJ      ! length of emission type joiner
        INTEGER         ML      ! length of MODELNAM buffer
        INTEGER         PINDX   ! polluntant index

        LOGICAL       :: DFLAG    = .FALSE.  ! true: process for diurnal emis
        LOGICAL       :: EFLAG    = .FALSE.  ! true: processing error found
        LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! true: first time routine called
        LOGICAL       :: IFLAG    = .FALSE.  ! true: input file available

        CHARACTER(LEN=IOVLEN3), SAVE :: VOLNAM = ' ' ! name of vol pollutant
        CHARACTER*300                   MESG         ! message buffer

        CHARACTER*16 :: PROGNAME = 'EFSETUP' ! program name

C***********************************************************************
C   begin body of subroutine EFSETUP

C.........  Set flag for diurnal or non-diurnal
        SELECT CASE( TYPNAM )
        CASE( 'DIURNAL' )
            DFLAG = .TRUE.

        CASE( 'NONDIURNAL' )
            DFLAG = .FALSE.

        CASE DEFAULT
            MESG = 'INTERNAL ERROR: TYPNAM must by "DIURNAL" or ' //
     &             '"NONDIURNAL" in program ' // PROGNAME
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END SELECT

C.........  Set flag for no input file available
        IFLAG = ( FNAME .EQ. 'NONE' ) 

C.........  Get length of model name
        ML = LEN_TRIM( MODELNAM )

C.........  Process for MOBILE5 model
        IF ( MODELNAM .EQ. 'MOBILE5' ) THEN

C.............  For new file, get environment variable for the volatile pol 
C               for mobile sources. For now, this routine only knows MOBILE5
            IF( IFLAG ) THEN
                IF( FIRSTIME ) THEN
                    MESG = 'Volatile pollutant type'
                    CALL ENVSTR( 'MB_HC_TYPE', MESG, 'VOC', VOLNAM, IOS)
                END IF

C.................  Create message about the volatile pollutant that is being 
C                   used
                L = LEN_TRIM( VOLNAM )
                IF ( IOS .LT. 0 ) THEN
                    MESG = 'WARNING: Default mobile volatile ' //
     &                     'pollutant "'// VOLNAM( 1:L ) // '" used.'

                ELSE
                    MESG = 'NOTE: Using mobile volatile pollutant "' //
     &                     VOLNAM( 1:L ) // '".'      
                END IF

C.............  For existing file, confirm MOBILE5 factors
            ELSE
c note: add here

            END IF

C.........  Abort if emission factor model name is not known
        ELSE
            MESG = 'Model name "' // MODELNAM( 1:ML ) // 
     &             '" not recognized by program ' // PROGNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Confirm that the pollutant of interest is valid for the model 
C           selected.
        IF( IFLAG ) THEN
            SELECT CASE( MODELNAM )

            CASE( 'MOBILE5' )

                L = LEN_TRIM( VOLNAM )
                PINDX = INDEX1( VOLNAM( 1:L ), NM5VPOL, M5VPOLS ) 

                IF( PINDX .LE. 0 ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                 'ERROR: Volatile pollutant type "' //
     &                 VOLNAM( 1:L ) // '" is invalid for the ' //
     &                 MODELNAM( 1:L ) // ' model'
                END IF
                  
            END SELECT
        END IF

C.........  Write out previously prepared message about the pollutant of 
C           interest.  This could have been set 3 sections above.
        IF( IFLAG ) CALL M3MSG2( MESG )

C.........  Abort if there has been a fatal problem up to this point
        IF( EFLAG ) THEN
            MESG = 'Problem configuring for ' // MODELNAM( 1:ML ) //
     &             ' emission factor model.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Store variable names, units, and descriptions, depending on the
C           model being used.
        SELECT CASE( MODELNAM )

        CASE( 'MOBILE5' )

C.................  Set maximum number of PSIs per scenario group
            MXPPGRP = MXM5SCEN

            IF( IFLAG ) THEN

C.................  Set variable names, units, description for CO and NOX from
C                   arrays defined in the MOBILE5 include file
                K = 0
                IF ( .NOT. DFLAG ) THEN
                    VNAMES( 1 ) = M5EFLST( 1 )   ! For CO
                    VNAMES( 2 ) = M5EFLST( 2 )   ! For NOX

                    VUNITS( 1 ) = M5EFUNT( 1 )
                    VUNITS( 2 ) = M5EFUNT( 2 )

                    VDESCS( 1 ) = M5EFDSC( 1 )
                    VDESCS( 2 ) = M5EFDSC( 2 )
                    K = 2
                END IF

C.................  Find volatile pol name in list of Mobile5 emission factor 
C                   names
C.................  Assign for diurnal OR non-diurnal (not both)
                LJ = LEN_TRIM( ETJOIN )
                DO I = 3, MXM5ALL

                    L  = INDEX( M5EFLST( I ), ETJOIN )
                    L2 = LEN_TRIM( M5EFLST( I ) )
                    IF( M5EFLST( I )( L+LJ:L2 ) .EQ. VOLNAM .AND.
     &                ( ( .NOT. DFLAG .AND. .NOT. M5DIURNL(I) ) .OR.
     &                  ( DFLAG .AND. M5DIURNL( I ) ) ) ) THEN
                        K = K + 1
                        VNAMES( K ) = M5EFLST( I )
                        VUNITS( K ) = M5EFUNT( I )
                        VDESCS( K ) = M5EFDSC( I )
                    END IF

                END DO

                IF( DFLAG ) THEN
                    NDIU = MXM5DIU
                ELSE
                    NNDI = MXM5NDI
                END IF

C.............  If file available
            ELSE 

                IF( .NOT. DESC3( FNAME ) ) THEN
                    MESG = 'Could not get description of ' //
     &                     TYPNAM // ' emission factors input file.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 ) 
                END IF

                VNAMES( 1:NVARS3D ) = VNAME3D( 1:NVARS3D )
                VUNITS( 1:NVARS3D ) = VNAME3D( 1:NVARS3D )
                VDESCS( 1:NVARS3D ) = VNAME3D( 1:NVARS3D )

                IF( DFLAG ) THEN
                    NDIU = NVARS3D
                ELSE
                    NNDI = NVARS3D
                END IF

            END IF

C.............  Allocate memory for and set the public variable for the 
C               non-diurnal and diurnal EF names...
C.............  For diurnal
            IF( DFLAG ) THEN
                ALLOCATE( DIUNAM( NDIU ), STAT=IOS )
                CALL CHECKMEM( IOS, 'DIUNAM', PROGNAME )
                ALLOCATE( DIUUNT( NDIU ), STAT=IOS )
                CALL CHECKMEM( IOS, 'NDIUNT', PROGNAME )
                ALLOCATE( DIUDSC( NDIU ), STAT=IOS )
                CALL CHECKMEM( IOS, 'NDIDSC', PROGNAME )

                DIUNAM( 1:NDIU ) = VNAMES( 1:NDIU )
                DIUUNT( 1:NDIU ) = VUNITS( 1:NDIU )
                DIUDSC( 1:NDIU ) = VDESCS( 1:NDIU )

C.............  Otherwise, for non-diurnal
            ELSE 
                ALLOCATE( NDINAM( NNDI ), STAT=IOS )
                CALL CHECKMEM( IOS, 'NDINAM', PROGNAME )
                ALLOCATE( NDIUNT( NNDI ), STAT=IOS )
                CALL CHECKMEM( IOS, 'NDIUNT', PROGNAME )
                ALLOCATE( NDIDSC( NNDI ), STAT=IOS )
                CALL CHECKMEM( IOS, 'NDIDSC', PROGNAME )

                NDINAM( 1:NNDI ) = VNAMES( 1:NNDI )
                NDIUNT( 1:NNDI ) = VUNITS( 1:NNDI )
                NDIDSC( 1:NNDI ) = VDESCS( 1:NNDI )

            END IF

        END SELECT

        FIRSTIME = .FALSE.

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

        END SUBROUTINE EFSETUP
