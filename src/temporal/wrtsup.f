
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
C      Modified by C. Coats 07/2014 
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
C.........  MODINFO  contains the information about the source category
C.........  MODSOURC contains the inventory arrays
C........   MODTMPRL contains the temporal profile tables

        USE MODINFO,  ONLY: CATEGORY
        USE MODSOURC, ONLY: CSOURC, CSCC, TPFLAG
        USE MODTMPRL, ONLY: METPROF,  MTHPROF,  WEKPROF,  DOMPROF, HRLPROF,
     &                                MTHIDP,   WEKIDP,   DOMIDP,  HRLIDP,
     &                      METCOUNT, MTHCOUNT, WEKCOUNT, DOMCOUNT,
     &                      MONCOUNT, TUECOUNT, WEDCOUNT, THUCOUNT,
     &                      FRICOUNT, SATCOUNT, SUNCOUNT

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   SUBROUTINE ARGUMENTS

        INTEGER     , INTENT(IN   ) :: FDEV             ! output file unit number
        INTEGER     , INTENT(IN   ) :: NSRC             ! number of sources
        INTEGER     , INTENT(IN   ) :: NVAR             ! number of variables
        CHARACTER(*), INTENT(IN   ) :: VARNAM( NVAR )   ! names of polltants/emis procs

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER, EXTERNAL :: INDEX1

C...........   PARAMETERs and their descriptions:

        INTEGER, PARAMETER  :: ACAT = 1
        INTEGER, PARAMETER  :: MCAT = 2
        INTEGER, PARAMETER  :: PCAT = 3

        CHARACTER( 6), PARAMETER :: LOCCATS( 3 ) = ( / 'AREA  ', 'MOBILE', 'POINT ' / )

        CHARACTER(16), PARAMETER :: PROGNAME = 'WRTSUP' !  program name

C.........  Local varables

        INTEGER       I, S, V, D        ! indices and counters
        INTEGER       ICAT
        INTEGER       PMTH              ! monthly profile from previous iteration
        INTEGER       PWEK              ! weekly  profile from previous iteration
        INTEGER       PDOM              ! weekly  profile from previous iteration
        INTEGER       PMON              ! hourly  profile from previous iteration
        INTEGER       PTUE              ! hourly  profile from previous iteration
        INTEGER       PWED              ! hourly  profile from previous iteration
        INTEGER       PTHU              ! hourly  profile from previous iteration
        INTEGER       PFRI              ! hourly  profile from previous iteration
        INTEGER       PSAT              ! hourly  profile from previous iteration
        INTEGER       PSUN              ! hourly  profile from previous iteration

        LOGICAL       MTHFLAG           ! true: monthly       same for all pols
        LOGICAL       WEKFLAG           ! true: weekly        same for all pols
        LOGICAL       DOMFLAG           ! true: day-of-month  same for all pols
        LOGICAL       MONFLAG           ! true: Monday hourly same for all pols
        LOGICAL       TUEFLAG           ! true: Monday hourly same for all pols
        LOGICAL       WEDFLAG           ! true: Monday hourly same for all pols
        LOGICAL       THUFLAG           ! true: Monday hourly same for all pols
        LOGICAL       FRIFLAG           ! true: Monday hourly same for all pols
        LOGICAL       SATFLAG           ! true: Monday hourly same for all pols
        LOGICAL       SUNFLAG           ! true: Monday hourly same for all pols

        CHARACTER(ALLLEN3+2)    SOURCE  ! "<fip><scc>..."
        CHARACTER(100)          OUTFMT  ! output format
        CHARACTER(512)          BUFFER  ! output variables buffer
        CHARACTER(128)          MESG

C***********************************************************************
C   begin body of subroutine WRTSUP

C.............  Ensure that the CATEGORY is valid

        ICAT = INDEX1( CATEGORY, 3, LOCCATS )

        IF ( ICAT .LE. 0 ) THEN
            MESG = 'INTERNAL ERROR: category "' // TRIM( CATEGORY ) //
     &             '" is not valid'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Write header with current variables

        BUFFER = '"' // TRIM( VARNAM( 1 ) ) // '"'
        DO V = 2, NVAR
            BUFFER = TRIM( BUFFER ) // ', "' // TRIM( VARNAM( V ) ) // '"'
        END DO

        WRITE( FDEV, '(A)' ) TRIM( BUFFER )

C.........  Create CSV output format

        OUTFMT = '(A, ",", I8, ",", 2X, 2049( A, :, "," ) )'   !  supports up through MXVARS3=2048 vbles

C.........  Loop through sources to output temporal profile info

        MTHFLAG = ( MTHCOUNT .GT. 0 )
        WEKFLAG = ( WEKCOUNT .GT. 0 )
        DOMFLAG = ( DOMCOUNT .GT. 0 )
        MONFLAG = ( MONCOUNT .GT. 0 )
        TUEFLAG = ( TUECOUNT .GT. 0 )
        WEDFLAG = ( WEDCOUNT .GT. 0 )
        THUFLAG = ( THUCOUNT .GT. 0 )
        FRIFLAG = ( FRICOUNT .GT. 0 )
        SATFLAG = ( SATCOUNT .GT. 0 )
        SUNFLAG = ( SUNCOUNT .GT. 0 )

        DO S = 1, NSRC

C.............  Retrieve profile numbers for all pollutants

            IF ( MTHFLAG )  THEN
                PMTH = MTHPROF( 1,1 )
                DO V = 2, NVAR
                    IF( MTHPROF( S,V )   .NE. PMTH ) MTHFLAG = .FALSE.
                END DO
            END IF
            IF ( WEKFLAG )  THEN
                PWEK = WEKPROF( 1,1 )
                DO V = 2, NVAR
                    IF( WEKPROF( S,V )   .NE. PWEK ) WEKFLAG = .FALSE.
                END DO
            END IF
            IF ( DOMFLAG )  THEN
                PDOM = DOMPROF( 1,1 )
                DO V = 2, NVAR
                    IF( DOMPROF( S,V )   .NE. PDOM ) DOMFLAG = .FALSE.
                END DO
            END IF
            IF ( MONFLAG )  THEN
                PMON = HRLPROF( 1,1,1 )
                DO V = 2, NVAR
                    IF( HRLPROF( S,1,V ) .NE. PMON ) MONFLAG = .FALSE.
                END DO
            END IF
            IF ( TUEFLAG )  THEN
                PTUE = HRLPROF( 1,2,1 )
                DO V = 2, NVAR
                    IF( HRLPROF( S,2,V ) .NE. PTUE ) TUEFLAG = .FALSE.
                END DO
            END IF
            IF ( WEDFLAG )  THEN
                PWED = HRLPROF( 1,3,1 )
                DO V = 2, NVAR
                    IF( HRLPROF( S,3,V ) .NE. PWED ) WEDFLAG = .FALSE.
                END DO
            END IF
            IF ( THUFLAG )  THEN
                PTHU = HRLPROF( 1,4,1 )
                DO V = 2, NVAR
                    IF( HRLPROF( S,4,V ) .NE. PTHU ) THUFLAG = .FALSE.
                END DO
            END IF
            IF ( FRIFLAG )  THEN
                PFRI = HRLPROF( 1,5,1 )
                DO V = 2, NVAR
                    IF( HRLPROF( S,5,V ) .NE. PFRI ) FRIFLAG = .FALSE.
                END DO
            END IF
            IF ( SATFLAG )  THEN
                PSAT = HRLPROF( 1,6,1 )
                DO V = 2, NVAR
                    IF( HRLPROF( S,6,V ) .NE. PSAT ) SATFLAG = .FALSE.
                END DO
            END IF
            IF ( SUNFLAG )  THEN
                PSUN = HRLPROF( 1,7,1 )
                DO V = 2, NVAR
                    IF( HRLPROF( S,7,V ) .NE. PSUN ) SUNFLAG = .FALSE.
                END DO
            END IF

C.............  Write profile information by pollutant

            CALL MKSOURCE( ICAT, S, SOURCE )

            IF( MTHCOUNT .EQ. 0 ) THEN
                CONTINUE
            ELSE IF( MTHFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"MTH"',    1, TRIM( SOURCE ), MTHIDP( MTHPROF(S,1 ) )
            ELSE
                WRITE( FDEV,OUTFMT ) '"MTH"', NVAR, TRIM( SOURCE ), ( TRIM( MTHIDP( MTHPROF(S,V) ) ),V=1,NVAR )
            END IF

            IF( WEKCOUNT .EQ. 0 ) THEN
                CONTINUE
            ELSE IF( WEKFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"WEK"',    1, TRIM( SOURCE ), WEKIDP( WEKPROF(S,1 ) )
            ELSE
                WRITE( FDEV,OUTFMT ) '"WEK"', NVAR, TRIM( SOURCE ), ( TRIM( WEKIDP( WEKPROF(S,V) ) ),V=1,NVAR )
            END IF

            IF( DOMCOUNT .EQ. 0 ) THEN
                CONTINUE
            ELSE IF( DOMFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"DOM"',    1, TRIM( SOURCE ), DOMIDP( DOMPROF(S,1) )
            ELSE
                WRITE( FDEV,OUTFMT ) '"DOM"', NVAR, TRIM( SOURCE ), ( TRIM( DOMIDP( DOMPROF(S,V) ) ),V=1,NVAR )
            END IF

            IF( MONCOUNT .EQ. 0 ) THEN
                CONTINUE
            ELSE IF( MONFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"MON"',    1, TRIM( SOURCE ), HRLIDP( HRLPROF(S,1,1) )
            ELSE
                WRITE( FDEV,OUTFMT ) '"MON"', NVAR, TRIM( SOURCE ), ( TRIM( HRLIDP( HRLPROF(S,1,V) ) ),V=1,NVAR )
            END IF

            IF( TUECOUNT .EQ. 0 ) THEN
                CONTINUE
            ELSE IF( TUEFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"TUE"',    1, TRIM( SOURCE ), HRLIDP( HRLPROF(S,2,1) )
            ELSE
                WRITE( FDEV,OUTFMT ) '"TUE"', NVAR, TRIM( SOURCE ), ( TRIM( HRLIDP( HRLPROF(S,2,V) ) ),V=1,NVAR )
            END IF

            IF( WEDCOUNT .EQ. 0 ) THEN
                CONTINUE
            ELSE IF( WEDFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"WED"',    1, TRIM( SOURCE ), HRLIDP( HRLPROF(S,3,1) )
            ELSE
                WRITE( FDEV,OUTFMT ) '"WED"', NVAR, TRIM( SOURCE ), ( TRIM( HRLIDP( HRLPROF(S,3,V) ) ),V=1,NVAR )
            END IF

            IF( THUCOUNT .EQ. 0 ) THEN
                CONTINUE
            ELSE IF( THUFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"THU"',    1, TRIM( SOURCE ), HRLIDP( HRLPROF(S,4,1) )
            ELSE
                WRITE( FDEV,OUTFMT ) '"THU"', NVAR, TRIM( SOURCE ), ( TRIM( HRLIDP( HRLPROF(S,4,V) ) ),V=1,NVAR )
            END IF

            IF( FRICOUNT .EQ. 0 ) THEN
                CONTINUE
            ELSE IF( FRIFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"FRI"',    1, TRIM( SOURCE ), HRLIDP( HRLPROF(S,5,1) )
            ELSE
                WRITE( FDEV,OUTFMT ) '"FRI"', NVAR, TRIM( SOURCE ), ( TRIM( HRLIDP( HRLPROF(S,5,V) ) ),V=1,NVAR )
            END IF

            IF( SATCOUNT .EQ. 0 ) THEN
                CONTINUE
            ELSE IF( SATFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"SAT"',    1, TRIM( SOURCE ), HRLIDP( HRLPROF(S,6,1) )
            ELSE
                WRITE( FDEV,OUTFMT ) '"SAT"', NVAR, TRIM( SOURCE ), ( TRIM( HRLIDP( HRLPROF(S,6,V) ) ),V=1,NVAR )
            END IF

            IF( SUNCOUNT .EQ. 0 ) THEN
                CONTINUE
            ELSE IF( SUNFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"SUN"',    1, TRIM( SOURCE ), HRLIDP( HRLPROF(S,7,1) )
            ELSE
                WRITE( FDEV,OUTFMT ) '"SUN"', NVAR, TRIM( SOURCE ), ( TRIM( HRLIDP( HRLPROF(S,7,V) ) ),V=1,NVAR )
            END IF

        END DO  ! end source loop

        RETURN

        CLOSE( FDEV )

C******************  INTERNAL SUBPROGRAMS  *****************************

      CONTAINS


        SUBROUTINE MKSOURCE( ICAT, SRC, SOURCE )

            INTEGER,          INTENT(IN   ) :: ICAT, SRC
            CHARACTER(LEN=*), INTENT(  OUT) :: SOURCE

            CHARACTER(1), PARAMETER :: QUOTES = '"'
            CHARACTER(1), PARAMETER :: SLASH  = '/'

            CHARACTER(FIPLEN3)  CFIP            !  tmp (character) FIPS code
            CHARACTER(SCCLEN3)  TSCC            !  tmp 10-digit SCC
            CHARACTER(CHRLEN3)  CPNT            !  tmp point ID
            CHARACTER(CHRLEN3)  CPLT            !  tmp plant ID
            CHARACTER(CHRLEN3)  CSTK            !  tmp stack ID
            CHARACTER(CHRLEN3)  CSEG            !  tmp segment ID
            CHARACTER(CHRLEN3)  CPL5            !  tmp plt char 5
            INTEGER             L

            TSCC    = CSCC( S )
            CFIP    = CSOURC( S )

            IF ( ICAT .EQ. ACAT ) THEN

                SOURCE = QUOTES // TRIM( ADJUSTL( CFIP ) ) //
     &                    SLASH // TRIM( ADJUSTL( TSCC ) ) // QUOTES

            ELSE IF ( ICAT .EQ. MCAT ) THEN

                SOURCE = QUOTES // TRIM( ADJUSTL( CFIP ) ) //
     &                    SLASH // TRIM( ADJUSTL( TSCC ) ) // QUOTES

            ELSE IF ( ICAT .EQ. PCAT ) THEN

                CPLT = CSOURC( S )( PTBEGL3(2):PTENDL3(2) )
                CPNT = CSOURC( S )( PTBEGL3(3):PTENDL3(3) )
                CSTK = CSOURC( S )( PTBEGL3(4):PTENDL3(4) )
                CSEG = CSOURC( S )( PTBEGL3(5):PTENDL3(5) )
                CPL5 = CSOURC( S )( PTBEGL3(6):PTENDL3(6) )

                SOURCE = QUOTES // TRIM( ADJUSTL( CFIP ) ) //
     &                    SLASH // TRIM( ADJUSTL( TSCC ) ) //
     &                    SLASH // TRIM( ADJUSTL( CPLT ) ) //
     &                    SLASH // TRIM( ADJUSTL( CPNT ) ) //
     &                    SLASH // TRIM( ADJUSTL( CSTK ) ) //
     &                    SLASH // TRIM( ADJUSTL( CSEG ) ) //
     &                    SLASH // TRIM( ADJUSTL( CPL5 ) ) // QUOTES

            ELSE

                CALL M3EXIT( 'WRTSUP/MKSOURCE', 0, 0, 'Illegal source category', 2 )

            END IF

            RETURN

        END SUBROUTINE MKSOURCE


      END SUBROUTINE WRTSUP
