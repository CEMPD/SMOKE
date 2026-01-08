
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
C      Modified by C. Coats 07/2014:  new GENTPRO CSV profiles and cross-references
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***************************************************************************

C...........   MODULES for public variables
C.........  MODINFO  contains the information about the source category
C.........  MODSOURC contains the inventory arrays
C........   MODTMPRL contains the temporal profile tables

        USE M3UTILIO

        USE MODINFO,  ONLY: CATEGORY
        USE MODSOURC, ONLY: CSOURC, CSCC, TPFLAG
        USE MODDAYHR, ONLY: INDXH, NHRSRC, INDXD, NDYSRC
        USE MODTMPRL, ONLY: METPROF,  MTHPROF,  WEKPROF,  DOMPROF, HRLPROF,
     &                                MTHIDP,   WEKIDP,   DOMIDP,  HRLIDP,
     &                      METCOUNT, MTHCOUNT, WEKCOUNT, DOMCOUNT,
     &                      MONCOUNT, TUECOUNT, WEDCOUNT, THUCOUNT,
     &                      FRICOUNT, SATCOUNT, SUNCOUNT

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   SUBROUTINE ARGUMENTS

        INTEGER     , INTENT(IN   ) :: FDEV             ! output file unit number
        INTEGER     , INTENT(IN   ) :: NSRC             ! number of sources
        INTEGER     , INTENT(IN   ) :: NVAR             ! number of variables
        CHARACTER(*), INTENT(IN   ) :: VARNAM( NVAR )   ! names of polltants/emis procs

C...........   EXTERNAL FUNCTIONS and their descriptions:

C       INTEGER, EXTERNAL :: INDEX1
C       INTEGER, EXTERNAL :: FIND1

C...........   PARAMETERs and their descriptions:

        INTEGER, PARAMETER  :: ACAT = 1
        INTEGER, PARAMETER  :: MCAT = 2
        INTEGER, PARAMETER  :: PCAT = 3

        CHARACTER( 6), PARAMETER :: LOCCATS( 3 ) = ( / 'AREA  ', 'MOBILE', 'POINT ' / )

        CHARACTER(16), PARAMETER :: PROGNAME = 'WRTSUP' !  program name

C.........  Local varables

        INTEGER       I, S, V, D        ! indices and counters
        INTEGER       ICAT, IDY, IHR

        INTEGER       PMTH              ! monthly   profile from previous iteration
        INTEGER       PDOM              ! daily     profile from previous iteration
        INTEGER       PWEK              ! weekly    profile from previous iteration
        INTEGER       PMET              ! met based profile from previous iteration
        INTEGER       PMON              ! diurnal   profile from previous iteration
        INTEGER       PTUE              ! diurnal   profile from previous iteration
        INTEGER       PWED              ! diurnal   profile from previous iteration
        INTEGER       PTHU              ! diurnal   profile from previous iteration
        INTEGER       PFRI              ! diurnal   profile from previous iteration
        INTEGER       PSAT              ! diurnal   profile from previous iteration
        INTEGER       PSUN              ! diurnal   profile from previous iteration

        INTEGER       NMTH              ! active monthly-profile count for this source
        INTEGER       NDOM              ! active   daily-profile count for this source
        INTEGER       NWEK              ! active  weekly-profile count for this source
        INTEGER       NMET              ! active met based prof  count for this source
        INTEGER       NMON              ! active diurnal profile count for this source
        INTEGER       NTUE              ! active diurnal profile count for this source
        INTEGER       NWED              ! active diurnal profile count for this source
        INTEGER       NTHU              ! active diurnal profile count for this source
        INTEGER       NFRI              ! active diurnal profile count for this source
        INTEGER       NSAT              ! active diurnal profile count for this source
        INTEGER       NSUN              ! active diurnal profile count for this source

        INTEGER       MTHP( NVAR )      ! active   monthly profiles for this source
        INTEGER       DOMP( NVAR )      ! active     daily profiles for this source
        INTEGER       WEKP( NVAR )      ! active    weekly profiles for this source
        INTEGER       METP( NVAR )      ! active met based-profiles for this source
        INTEGER       HRLP( 7,NVAR )    ! active   diurnal profiles for this source

        LOGICAL       MTHFLAG           ! true: monthly       same for all pols
        LOGICAL       WEKFLAG           ! true: weekly        same for all pols
        LOGICAL       DOMFLAG           ! true: day-of-month  same for all pols
        LOGICAL       METFLAG           ! true: met based same for all pols
        LOGICAL       MONFLAG           ! true: Monday hourly same for all pols
        LOGICAL       TUEFLAG           ! true: Monday hourly same for all pols
        LOGICAL       WEDFLAG           ! true: Monday hourly same for all pols
        LOGICAL       THUFLAG           ! true: Monday hourly same for all pols
        LOGICAL       FRIFLAG           ! true: Monday hourly same for all pols
        LOGICAL       SATFLAG           ! true: Monday hourly same for all pols
        LOGICAL       SUNFLAG           ! true: Monday hourly same for all pols

        CHARACTER(TMPLEN3)      INVPRF        ! tmp inv tmp profile ID
        CHARACTER(TMPLEN3)      MTHPRF(NVAR)  ! monthly prof for inv source with day/hourly inv
        CHARACTER(TMPLEN3)      DOMPRF(NVAR)  ! dayofmon prof for inv source with day/hourly inv
        CHARACTER(TMPLEN3)      WEKPRF(NVAR)  ! weekly prof for inv source with day/hourly inv
        CHARACTER(TMPLEN3)      MONPRF(NVAR)  ! monday hourly prof for inv source with day/hourly inv
        CHARACTER(TMPLEN3)      TUEPRF(NVAR)  ! tuesday hourly prof for inv source with day/hourly inv
        CHARACTER(TMPLEN3)      WEDPRF(NVAR)  ! wednesday hourly prof for inv source with day/hourly inv
        CHARACTER(TMPLEN3)      THUPRF(NVAR)  ! thursday hourly prof for inv source with day/hourly inv
        CHARACTER(TMPLEN3)      FRIPRF(NVAR)  ! friday hourly prof for inv source with day/hourly inv
        CHARACTER(TMPLEN3)      SATPRF(NVAR)  ! saturday hourly prof for inv source with day/hourly inv
        CHARACTER(TMPLEN3)      SUNPRF(NVAR)  ! sunday hourly prof for inv source with day/hourly inv
        CHARACTER(TMPLEN3)      METPRF(NVAR)  ! met hourly prof for inv source with day/hourly inv
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

C.........  Create CSV output formats

        OUTFMT = '(A, ",", I8, ",", 2X,    2049( A,  :, "," ) )'   !  supports up through MXVARS3=2048 vbles

        DO S = 1, NSRC

C.............  Build source charaterisitics
            CALL MKSOURCE( ICAT, S, SOURCE )

C.............  Search index for inv source with day/hour inventory
            IDY = FIND1( S, NDYSRC, INDXD )
            IHR = FIND1( S, NHRSRC, INDXH )

C.............  Retrieve profile numbers for all pollutants
C.............  Implement profile-hierarchy and count active profiles
            MTHFLAG = ( MTHCOUNT .GT. 0 )
            DOMFLAG = ( DOMCOUNT .GT. 0 )
            WEKFLAG = ( WEKCOUNT .GT. 0 )
            METFLAG = ( METCOUNT .GT. 0 )
            MONFLAG = ( MONCOUNT .GT. 0 )
            TUEFLAG = ( TUECOUNT .GT. 0 )
            WEDFLAG = ( WEDCOUNT .GT. 0 )
            THUFLAG = ( THUCOUNT .GT. 0 )
            FRIFLAG = ( FRICOUNT .GT. 0 )
            SATFLAG = ( SATCOUNT .GT. 0 )
            SUNFLAG = ( SUNCOUNT .GT. 0 )

            NMTH = 0
            NDOM = 0
            NWEK = 0
            NMET = 0
            NMON = 0
            NTUE = 0
            NWED = 0
            NTHU = 0
            NFRI = 0
            NSAT = 0
            NSUN = 0

            DO V = 1, NVAR

                MTHP(V)   = MTHPROF( S,V )
                DOMP(V)   = DOMPROF( S,V )
                WEKP(V)   = WEKPROF( S,V )
                METP(V)   = METPROF( S,V )
                HRLP(:,V) = HRLPROF( S,:,V )

                IF ( METP(V) .GT. 0 ) THEN
                    MTHP(V)   = 0
                    DOMP(V)   = 0
                    WEKP(V)   = 0
                    HRLP(:,V) = 0
                END IF

                IF ( DOMP(V) .GT. 0 ) THEN
                    WEKP(V)   = 0
                END IF

                IF( MTHP( V )   .NE. MTHP( 1   ) ) MTHFLAG = .FALSE.
                IF( DOMP( V )   .NE. DOMP( 1   ) ) DOMFLAG = .FALSE.
                IF( WEKP( V )   .NE. WEKP( 1   ) ) WEKFLAG = .FALSE.
                IF( METP( V )   .NE. METP( 1   ) ) METFLAG = .FALSE.
                IF( HRLP( 1,V ) .NE. HRLP( 1,1 ) ) MONFLAG = .FALSE.
                IF( HRLP( 2,V ) .NE. HRLP( 2,1 ) ) TUEFLAG = .FALSE.
                IF( HRLP( 3,V ) .NE. HRLP( 3,1 ) ) WEDFLAG = .FALSE.
                IF( HRLP( 4,V ) .NE. HRLP( 4,1 ) ) THUFLAG = .FALSE.
                IF( HRLP( 5,V ) .NE. HRLP( 5,1 ) ) FRIFLAG = .FALSE.
                IF( HRLP( 6,V ) .NE. HRLP( 6,1 ) ) SATFLAG = .FALSE.
                IF( HRLP( 7,V ) .NE. HRLP( 7,1 ) ) SUNFLAG = .FALSE.

                IF( MTHP( V )   .GT. 0 )  NMTH = NMTH + 1
                IF( DOMP( V )   .GT. 0 )  NDOM = NDOM + 1
                IF( WEKP( V )   .GT. 0 )  NWEK = NWEK + 1
                IF( METP( V )   .GT. 0 )  NMET = NMET + 1
                IF( HRLP( 1,V ) .GT. 0 )  NMON = NMON + 1
                IF( HRLP( 2,V ) .GT. 0 )  NTUE = NTUE + 1
                IF( HRLP( 3,V ) .GT. 0 )  NWED = NWED + 1
                IF( HRLP( 4,V ) .GT. 0 )  NTHU = NTHU + 1
                IF( HRLP( 5,V ) .GT. 0 )  NFRI = NFRI + 1
                IF( HRLP( 6,V ) .GT. 0 )  NSAT = NSAT + 1
                IF( HRLP( 7,V ) .GT. 0 )  NSUN = NSUN + 1

                IF( IDY > 0 .OR. IHR > 0 ) THEN
                    IF( IDY > 0 ) WRITE( INVPRF, '(A,I8.8)' ) 'DY', S
                    IF( IHR > 0 ) WRITE( INVPRF, '(A,I8.8)' ) 'HR', S
                    MTHPRF( V ) = INVPRF
                    DOMPRF( V ) = INVPRF
                    WEKPRF( V ) = INVPRF
                    METPRF( V ) = INVPRF
                    MONPRF( V ) = INVPRF
                    TUEPRF( V ) = INVPRF
                    WEDPRF( V ) = INVPRF
                    THUPRF( V ) = INVPRF
                    FRIPRF( V ) = INVPRF
                    SATPRF( V ) = INVPRF
                    SUNPRF( V ) = INVPRF
                ELSE
                    MTHPRF( V ) = MTHIDP( MTHP(V) ) 
                    DOMPRF( V ) = DOMIDP( DOMP(V) ) 
                    WEKPRF( V ) = WEKIDP( WEKP(V) ) 
                    METPRF( V ) = SOURCE( 2:FIPLEN3+1 ) 
                    MONPRF( V ) = HRLIDP( HRLP(1,V) )
                    TUEPRF( V ) = HRLIDP( HRLP(2,V) )
                    WEDPRF( V ) = HRLIDP( HRLP(3,V) )
                    THUPRF( V ) = HRLIDP( HRLP(4,V) )
                    FRIPRF( V ) = HRLIDP( HRLP(5,V) )
                    SATPRF( V ) = HRLIDP( HRLP(6,V) )
                    SUNPRF( V ) = HRLIDP( HRLP(7,V) )
                END IF
            END DO

C.............  Write profile information by pollutant

            IF( MOD( TPFLAG(S), MTPRFAC ) .NE. 0 ) THEN
                CONTINUE
            ELSE IF( NMTH .EQ. 0 ) THEN
                CONTINUE
            ELSE IF( MTHFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"MTH"',    1, TRIM( SOURCE ), MTHPRF( 1 )
            ELSE
                WRITE( FDEV,OUTFMT ) '"MTH"', NVAR, TRIM( SOURCE ), ( TRIM( MTHPRF( V ) ),V=1,NVAR )
            END IF

            IF( NDOM .EQ. 0 ) THEN
                CONTINUE
            ELSE IF( DOMFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"DOM"',    1, TRIM( SOURCE ), DOMPRF( 1 )
            ELSE
                WRITE( FDEV,OUTFMT ) '"DOM"', NVAR, TRIM( SOURCE ), ( TRIM( DOMPRF( V ) ),V=1,NVAR )
            END IF

            IF( NWEK .EQ. 0  ) THEN
                CONTINUE
            ELSE IF( WEKFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"WEK"',    1, TRIM( SOURCE ), WEKPRF( 1 )
            ELSE
                WRITE( FDEV,OUTFMT ) '"WEK"', NVAR, TRIM( SOURCE ), ( TRIM( WEKPRF( V ) ),V=1,NVAR )
            END IF

            IF( NMET .EQ. 0 ) THEN
                CONTINUE
            ELSE IF( METFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"MET"',    1, TRIM( SOURCE ), METPRF( 1 )
            ELSE
                WRITE( FDEV,OUTFMT ) '"MET"', NVAR, TRIM( SOURCE ), ( TRIM( METPRF( V ) ), V=1,NVAR )
            END IF

            IF( NMON .EQ. 0 ) THEN
                CONTINUE
            ELSE IF( MONFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"MON"',    1, TRIM( SOURCE ), MONPRF( 1 )
            ELSE
                WRITE( FDEV,OUTFMT ) '"MON"', NVAR, TRIM( SOURCE ), ( TRIM( MONPRF( V ) ),V=1,NVAR )
            END IF

            IF( NTUE .EQ. 0 ) THEN
                CONTINUE
            ELSE IF( TUEFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"TUE"',    1, TRIM( SOURCE ), TUEPRF( 1 )
            ELSE
                WRITE( FDEV,OUTFMT ) '"TUE"', NVAR, TRIM( SOURCE ), ( TRIM( TUEPRF( V ) ),V=1,NVAR )
            END IF

            IF( NWED .EQ. 0 ) THEN
                CONTINUE
            ELSE IF( HRLPROF( S,1,1 ) .LT. 1 ) THEN
                CONTINUE
            ELSE IF( WEDFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"WED"',    1, TRIM( SOURCE ), WEDPRF( 1 )
            ELSE
                WRITE( FDEV,OUTFMT ) '"WED"', NVAR, TRIM( SOURCE ), ( TRIM( WEDPRF( V ) ),V=1,NVAR )
            END IF

            IF( NTHU .EQ. 0 ) THEN
                CONTINUE
            ELSE IF( THUFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"THU"',    1, TRIM( SOURCE ), THUPRF( 1 )
            ELSE
                WRITE( FDEV,OUTFMT ) '"THU"', NVAR, TRIM( SOURCE ), ( TRIM( THUPRF( V ) ),V=1,NVAR )
            END IF

            IF( NFRI .EQ. 0 ) THEN
                CONTINUE
            ELSE IF( FRIFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"FRI"',    1, TRIM( SOURCE ), FRIPRF( 1 )
            ELSE
                WRITE( FDEV,OUTFMT ) '"FRI"', NVAR, TRIM( SOURCE ), ( TRIM( FRIPRF( V ) ),V=1,NVAR )
            END IF

            IF( NSAT .EQ. 0 ) THEN
                CONTINUE
            ELSE IF( SATFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"SAT"',    1, TRIM( SOURCE ), SATPRF( 1 )
            ELSE
                WRITE( FDEV,OUTFMT ) '"SAT"', NVAR, TRIM( SOURCE ), ( TRIM( SATPRF( V ) ),V=1,NVAR )
            END IF

            IF( NSUN .EQ. 0 ) THEN
                CONTINUE
            ELSE IF( HRLPROF( S,1,1 ) .LT. 1 ) THEN
                CONTINUE
            ELSE IF( SUNFLAG ) THEN
                WRITE( FDEV,OUTFMT ) '"SUN"',    1, TRIM( SOURCE ), SUNPRF( 1 )
            ELSE
                WRITE( FDEV,OUTFMT ) '"SUN"', NVAR, TRIM( SOURCE ), ( TRIM( SUNPRF( V ) ),V=1,NVAR )
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
