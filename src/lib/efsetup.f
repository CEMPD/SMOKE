
        SUBROUTINE EFSETUP( FNAME, MODELNAM, MXVAR, NVAR, 
     &                      VNAMES, VUNITS, VDESCS, VOLNAM )
   
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
        INCLUDE 'M6CNST3.EXT'   !  Mobile6 constants
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data stru

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER     INDEX1
        EXTERNAL    INDEX1

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: FNAME           ! logical file or 'NONE'
        CHARACTER(*), INTENT (IN) :: MODELNAM        ! name of EF model
        INTEGER     , INTENT (IN) :: MXVAR           ! max no of variables
        INTEGER     , INTENT(OUT) :: NVAR            ! actual no of variables
        CHARACTER(*), INTENT(OUT) :: VNAMES( MXVAR ) ! variable names
        CHARACTER(*), INTENT(OUT) :: VUNITS( MXVAR ) ! variable units
        CHARACTER(*), INTENT(OUT) :: VDESCS( MXVAR ) ! variable descriptions
        CHARACTER(*), INTENT(OUT) :: VOLNAM          ! volatile pollutant name

C...........   Local variables
        INTEGER         I, J, K, L, L2, N ! counters and indices

        INTEGER         IOS     ! status from retrieving E.V.s
        INTEGER         LJ      ! length of emission type joiner
        INTEGER         ML      ! length of MODELNAM buffer
        INTEGER         PINDX   ! polluntant index

        LOGICAL       :: EFLAG    = .FALSE.  ! true: processing error found
        LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! true: first time routine called
        LOGICAL       :: IFLAG    = .FALSE.  ! true: input file available

        CHARACTER*300                   MESG         ! message buffer

        CHARACTER*16 :: PROGNAME = 'EFSETUP' ! program name

C***********************************************************************
C   begin body of subroutine EFSETUP

C.........  Set flag for no input file available
        IFLAG = ( FNAME .EQ. 'NONE' ) 

C.........  Get length of model name
        ML = LEN_TRIM( MODELNAM )

C.........  Process for MOBILE6 model
        IF ( MODELNAM .EQ. 'MOBILE6' ) THEN

C.............  For new file, get environment variable for the volatile pol 
C               for mobile sources. For now, this routine only knows MOBILE6
            IF( IFLAG ) THEN
                IF( FIRSTIME ) THEN
                    MESG = 'Volatile pollutant type'
                    CALL ENVSTR( 'MB_HC_TYPE', MESG, 'TOG', VOLNAM, IOS)
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

C.............  For existing file, confirm MOBILE6 names
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

            CASE( 'MOBILE6' )

                L = LEN_TRIM( VOLNAM )
                PINDX = INDEX1( VOLNAM( 1:4 ), NM6VPOL, M6VPOLS ) 

                IF( PINDX .LE. 0 ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                 'ERROR: Volatile pollutant type "' //
     &                 VOLNAM( 1:L ) // '" is invalid for the ' //
     &                 MODELNAM( 1:ML ) // ' model'
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

        CASE( 'MOBILE6' )

            IF( IFLAG ) THEN

C.................  Set variable names, units, and descriptions from
C                   arrays defined in the MOBILE6 include file
  
                K = 0

C.................  Loop over all emission processes
                DO I = 1, MXM6EPR
                
C.....................  Loop over all pollutants                
                    DO J = 1, MXM6POLS
  
C.........................  Check if this is a valid pollutant/process combo
                        IF( M6POL2EF( I,J ) == -1 ) CYCLE
                        
                        K = K + 1
    
C.........................  Pull process name and description from include file                        
                        VNAMES( K ) = M6PROCS( I ) // ETJOIN
                        VDESCS( K ) = 'EFs for ' // M6PRCDSC( I )
    
C.........................  If pollutant is HC, append specified volatile pollutant
C                           name; otherwise, use name and desc from include file
                        IF( TRIM( M6POLS( J ) ) == 'HC' ) THEN
                            VNAMES( K ) = TRIM( VNAMES( K ) ) // VOLNAM
                            VDESCS( K ) = TRIM( VDESCS( K ) ) // ' ' // 
     &                                    VOLNAM
                        ELSE
                            VNAMES( K ) = TRIM( VNAMES( K ) ) //
     &                                    M6POLS( J )
                            VDESCS( K ) = TRIM( VDESCS( K ) ) // ' ' //
     &                                    M6POLDSC( J )
                        END IF
    
C.........................  Store units from include file                        
                        VUNITS( K ) = M6UNIT
                        
                    END DO  ! pollutant loop
                END DO  ! emission process loop

C.................  Store maximum number of emission factor types
                NEFS = MXM6EFS

C.............  If file available
            ELSE 

                IF( .NOT. DESC3( FNAME ) ) THEN
                    MESG = 'Could not get description of ' //
     &                     'emission factors input file.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 ) 
                END IF

                VNAMES( 1:NVARS3D ) = VNAME3D( 1:NVARS3D )
                VUNITS( 1:NVARS3D ) = VNAME3D( 1:NVARS3D )
                VDESCS( 1:NVARS3D ) = VNAME3D( 1:NVARS3D )

                NEFS = NVARS3D

            END IF

C.............  Allocate memory for and set the public variables for the 
C               EF names, units, and descriptions...

            ALLOCATE( EFSNAM( NEFS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EFSNAM', PROGNAME )
            ALLOCATE( EFSUNT( NEFS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EFSUNT', PROGNAME )
            ALLOCATE( EFSDSC( NEFS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EFSDSC', PROGNAME )

            EFSNAM( 1:NEFS ) = VNAMES( 1:NEFS )
            EFSUNT( 1:NEFS ) = VUNITS( 1:NEFS )
            EFSDSC( 1:NEFS ) = VDESCS( 1:NEFS )

        END SELECT

        FIRSTIME = .FALSE.

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

        END SUBROUTINE EFSETUP
