
        SUBROUTINE WRTSUP( FDEV, NSRC, NVAR, VARNAM )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine writes the temporal supplemental ASCII file
C
C  PRECONDITIONS REQUIRED:
C      Outfile file is opened
C      Used modules are populated for use by temporal
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines
C      Functions: I/O API functions
C
C  REVISION  HISTORY:
C      Created by M. Houyoux 10/2001
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
C***************************************************************************

C...........   MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC, ONLY: TPFLAG

C...........   This module contains the cross-reference tables
        USE MODXREF, ONLY: MDEX, WDEX, DDEX

C...........   This module contains the temporal profile tables
        USE MODTMPRL, ONLY: MONREF, WEKREF, HRLREF

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:

C...........   SUBROUTINE ARGUMENTS
        INTEGER         FDEV             ! output file unit number
        INTEGER         NSRC             ! number of sources
        INTEGER         NVAR             ! number of variables
        CHARACTER(*)    VARNAM( NVAR )   ! names of polltants/emis procs

C.........  Local allocatable arrays
        INTEGER, ALLOCATABLE :: MONPROF( : )  ! tmp monthly profiles by variable
        INTEGER, ALLOCATABLE :: WEKPROF( : )  ! tmp weekly profiles by variable
        INTEGER, ALLOCATABLE :: HRLPROF( : )  ! tmp diurnal profiles by variable

C.........  Local varables
        INTEGER       I, L1, L2, S, V      ! indices and counters

        INTEGER       IOS               ! i/o status
        INTEGER       PMON              ! monthly profile from previous iteration
        INTEGER       PWEK              ! weekly profile from previous iteration
        INTEGER       PHRL              ! hourly profile from previous iteration

        LOGICAL       MFLAG             ! true: monthly same for all pols
        LOGICAL       WFLAG             ! true: weekly same for all pols
        LOGICAL       HFLAG             ! true: hourly same for all pols

        CHARACTER(100) :: OUTFMT = ' '    ! output format
        CHARACTER(512) :: BUFFER = ' '    ! output variables buffer

        CHARACTER(16) :: PROGNAME = 'WRTSUP' !  program name

C***********************************************************************
C   begin body of subroutine WRTSUP

C.........  Allocate local temporary variables
        ALLOCATE( MONPROF( NVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MONPROF', PROGNAME )
        ALLOCATE( WEKPROF( NVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WEKPROF', PROGNAME )
        ALLOCATE( HRLPROF( NVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'HRLPROF', PROGNAME )

C.........  Write header with current variables
	L1 = LEN_TRIM( VARNAM( 1 ) )
        BUFFER = '"' // VARNAM( 1 )( 1:L1 ) // '"'
        DO V = 2, NVAR
            L1 = LEN_TRIM( VARNAM( V ) )
            L2 = LEN_TRIM( BUFFER )
            BUFFER = BUFFER(1:L2) // ', "' // VARNAM( V )(1:L1) // '"'
        END DO

        L2 = LEN_TRIM( BUFFER )
        WRITE( FDEV, '(A)' ) BUFFER( 1:L2 )

C.........  Create output format
        OUTFMT = '(A,I8,","'
        DO I = 1, NVAR, 50
            IF ( I .EQ. 1 ) THEN
                OUTFMT = TRIM( OUTFMT ) // ',50(I8,",")'
            ELSE
                OUTFMT = TRIM( OUTFMT ) // ',14X,/,50(I8,",")'
            END IF
        END DO 
        OUTFMT = TRIM( OUTFMT ) // ',I8)'

C.........  Loop through sources to output temporal profile info
        DO S = 1, NSRC

C.............  Retrieve profile numbers for all pollutants
            MFLAG = .TRUE.
            WFLAG = .TRUE.
            HFLAG = .TRUE.
            PMON = 0
            PWEK = 0
            PHRL = 0
            DO V = 1, NVAR
                MONPROF( V ) = MONREF( MDEX( S,V ) )
		WEKPROF( V ) = WEKREF( WDEX( S,V ) )
		HRLPROF( V ) = HRLREF( DDEX( S,V ) )

                IF( V .NE. 1 .AND.  
     &              MONPROF( V ) .NE. PMON ) MFLAG = .FALSE.
                IF( V .NE. 1 .AND.
     &              WEKPROF( V ) .NE. PWEK ) WFLAG = .FALSE.
                IF( V .NE. 1 .AND.
     &              HRLPROF( V ) .NE. PHRL ) HFLAG = .FALSE.

                PMON = MONPROF( V )
                PWEK = WEKPROF( V )
                PHRL = HRLPROF( V )

C.................  If source and pollutant is hour-specific, change hourly 
C                   profile to negative
c note: to do this, will need to have a list of all sources that are hour-specific
c    n: across all hours.
c                IF( LHSPOA( V ) .AND. ???

C.................  If source and pollutant is day-specific, change weekly 
C                   profile to negative
c note: same note as above.

            END DO

C.............  If source does not use monthly profile, change monthly 
C               profile to negative
            IF( MOD( TPFLAG( S ), MTPRFAC ) .NE. 0      ) THEN

                MONPROF = -MONPROF    ! array

            END IF

C.............  Write profile information by pollutant
            IF( MFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"M"', 1, MONPROF( 1 )
            ELSE
                WRITE( FDEV,OUTFMT ) '"M"', NVAR, (MONPROF(V),V=1,NVAR)
            END IF

            IF( MFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"W"', 1, WEKPROF( 1 )
            ELSE
                WRITE( FDEV,OUTFMT ) '"W"', NVAR, (WEKPROF(V),V=1,NVAR)
            END IF

            IF( MFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"H"', 1, HRLPROF( 1 )
            ELSE
                WRITE( FDEV,OUTFMT ) '"H"', NVAR, (HRLPROF(V),V=1,NVAR)
            END IF

        END DO  ! end source loop

C.........  Deallocate local temporary variables
	DEALLOCATE( MONPROF, WEKPROF, HRLPROF )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94100   FORMAT( 9( A, I2.2 ) )

        END SUBROUTINE WRTSUP
