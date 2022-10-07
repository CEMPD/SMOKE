
        SUBROUTINE ASGNTPRO

C***********************************************************************
C  subroutine body starts at line 145
Cc
C  DESCRIPTION:
C      For each source and pollutant or emission type, find the most specific
C      temporal profile that applies to that source. Do this using the
C      grouped tables of temporal cross references from RDTREF. The hierarchical
C      order is defined in this subroutine, and can be determined from the
C      in-source comments below. Once a profile code has been identified,
C      search for this code in the temporal profile tables (from RDTPROF) and
C      save the index to these tables for each source and pollutant.
C
C  PRECONDITIONS REQUIRED:
C     Call PROCTPRO() first
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 1/1999 by M. Houyoux
C
C       Revised ??/???? by ??
C
C       Version 07/2014 by C.Coats for  new GENTPRO CSV profiles and cross-references
C       Does source-to-profile mapping for all sources,species.
C       NOTES:
C       1) Met-based profiles are managed in PROCTPRO()
C       2) Uses profile-zero and xref-zero entries as sentinels for
C          the "profile/xref not supplied" cases.
C****************************************************************************
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
Ct
C***************************************************************************

C.........   MODULES for public variables
C.........   MODSOURC contains the source ararys
C.........   MODTMPRL contains the temporal profile tables
C.........   MODINFO contains the information about the source category

        USE MODSOURC, ONLY: CSOURC, CSCC, TPFLAG, IRCLAS, IVTYPE

        USE MODTMPRL, ONLY: NMON, NWEK, NHRL, NDOM, METPROF,
     &                      MTHPROF, WEKPROF, DOMPROF, HRLPROF,
     &                      MTHCOUNT, WEKCOUNT, DOMCOUNT,
     &                      MONCOUNT, TUECOUNT, WEDCOUNT, THUCOUNT,
     &                      FRICOUNT, SATCOUNT, SUNCOUNT, METCOUNT,
     &                      MTHPDEX, WEKPDEX, DOMPDEX,
     &                      MONPDEX, TUEPDEX, WEDPDEX, THUPDEX,
     &                      FRIPDEX, SATPDEX, SUNPDEX,
     &                      MTHKEYS, WEKKEYS, DOMKEYS,
     &                      MONKEYS, TUEKEYS, WEDKEYS, THUKEYS,
     &                      FRIKEYS, SATKEYS, SUNKEYS, METKEYS,
     &                      DAYFLAG, POLREFFLAG, METREFFLAG

        USE MODINFO, ONLY: NSRC, NCHARS, NIPPA, EANAM, LSCCEND, CATEGORY

        IMPLICIT NONE

C.........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api constant parameters

C.........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER(2), EXTERNAL :: CRLF
        INTEGER     , EXTERNAL :: ENVINT
        LOGICAL     , EXTERNAL :: ENVYN
        INTEGER     , EXTERNAL :: FIND1
        INTEGER     , EXTERNAL :: FINDC
        INTEGER     , EXTERNAL :: INDEX1
        LOGICAL     , EXTERNAL :: SETSCCTYPE
        LOGICAL     , EXTERNAL :: CHKEXPSCC
        LOGICAL     , EXTERNAL :: USEEXPGEO

C..........   Local parameters

        CHARACTER(9), PARAMETER :: DAYNAME( 7 ) =
     &      (/  'MONDAY   ', 'TUESDAY  ', 'WEDNESDAY', 'THURSDAY ',
     &          'FRIDAY   ', 'SATURDAY ', 'SUNDAY   '   /)

        CHARACTER( 1),      PARAMETER :: BLANK   = ' '
        CHARACTER(24),      PARAMETER :: ZEROS   = '000000000000000000000000'
        CHARACTER(24),      PARAMETER :: BADSTR  = '????????????????????????'
        CHARACTER(RWTLEN3), PARAMETER :: RWTZERO = ZEROS    !  zero roadway type
        CHARACTER(VIDLEN3), PARAMETER :: VIDZERO = ZEROS    !  zero vehicle type
        CHARACTER(FIPLEN3), PARAMETER :: CFIPZ   = ZEROS    !  zero Country/state code
        CHARACTER(SCCLEN3), PARAMETER :: TSCCZ   = ZEROS    !  zero TSCC
        CHARACTER(CHRLEN3), PARAMETER :: BLNK    = BLANK    !  tmp blank
        CHARACTER(PLTLEN3), PARAMETER :: BLNKPLT = BLANK    !  tmp blank

        CHARACTER(16),      PARAMETER :: PNAME   = 'ASGNTPRO' ! program name

C.........  Other local variables

        INTEGER         I, J, L2, S, V      !  counters and indices
        INTEGER         ISTAT
        INTEGER         IPROF
        INTEGER         HCOUNT

        INTEGER          ICOUNT, ERRCNT  !  count of errors
        INTEGER          F0, F1, F2, F3, F4, F5, F6  ! tmp find indices
        INTEGER          F0B     ! extra find index for mobile
        INTEGER          F2B     ! extra find index for mobile
        INTEGER          F4B     ! extra find index for mobile

        LOGICAL          EFLAG
        LOGICAL          MFLAG              !  true: use monthly profiles
        LOGICAL          WFLAG              !  true: use weekly  profiles
        LOGICAL          SCCFLAG            !  true: SCC type is different from previous

        INTEGER, SAVE :: MXERR              !  max error messages to output
        LOGICAL, SAVE :: FIRSTIME = .TRUE.  !  true: first time routine called
        LOGICAL, SAVE :: FULLSCC  = .FALSE. ! true: use only full SCC entries

        CHARACTER(5)        CPOS( 0:NIPPA ) !  string for sorted position of pol/act

        CHARACTER(300)      MESG            !  message buffer
        CHARACTER(300)      CBUF            !  CSRC line buffer
        CHARACTER(SRCLEN3)  CSRC            !  tmp source chars string
        CHARACTER(FIPLEN3)  CFIP            !  tmp (character) FIPS code
        CHARACTER(FIPLEN3)  CFIPL           !  tmp Country/state code
        CHARACTER(SCCLEN3)  TSCC            !  tmp 10-digit SCC
        CHARACTER(SCCLEN3)  TSCCL           !  tmp left digits of TSCC
        CHARACTER(SCCLEN3)  TSCC5           !  tmp left digits of TSCC
        CHARACTER(SCCLEN3)  CHKRWT          !  tmp roadway type only SCC
        CHARACTER(SCCLEN3)  CHKVID          !  tmp vehicle-type only SCC
        CHARACTER(RWTLEN3)  CRWT            !  tmp char roadway type
        CHARACTER(VIDLEN3)  CVID            !  tmp vehicle type
        CHARACTER(IOVLEN3)  CPOA            !  temporary pollutant/emission type
        CHARACTER(PLTLEN3)  CPLT            !  tmp plant ID
        CHARACTER(CHRLEN3)  CPNT            !  tmp point ID
        CHARACTER(CHRLEN3)  CSTK            !  tmp stack ID
        CHARACTER(CHRLEN3)  CSEG            !  tmp segment ID
        CHARACTER(CHRLEN3)  CPL5            !  tmp plt char 5 (==SCC)
        CHARACTER(CHRLEN3)  CPLL            !  tmp plt char 5 (==SCC)
        CHARACTER(CHRLEN3)  CPLZ            !  tmp plt char 5 (==SCC)


C.........  CSRCALL from BLDSRC(), listed in order of the search hierarchy
C.........  See SMOKE documentation, section 6.17:
C.........  https://www.cmascenter.org/smoke/documentation/3.5/html/ch06s17.html
C.........  Note that "no pollutant specified" initializes *PDEX arrays where possible;
C.........  then V=1,...,NIPPA override it
C.........  BADSTR initialization is retained for non-full-SCC entries
C.........  whenever FULLSCC is set.  (The logic is considerably simpler
C.........  to do all the searches, and use BADSTR to force the irrelevant
C.........  ones to fail.)

        CHARACTER(ALLLEN3) :: CSRC01 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC02 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC03 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC04 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC05 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC06 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC07 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC08 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC084= BADSTR
        CHARACTER(ALLLEN3) :: CSRC085= BADSTR
        CHARACTER(ALLLEN3) :: CSRC09 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC10 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC11 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC12 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC13 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC14 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC15 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC16 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC17 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC18 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC19 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC20 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC21 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC22 = BADSTR
        CHARACTER(ALLLEN3) :: CSRC23 = BADSTR

C***********************************************************************
C   begin body of subroutine ASGNTPRO

        EFLAG = .FALSE.

C............  For first time routine is called in all cases,

        IF( FIRSTIME ) THEN

C.............  Retrieve environment variables

            MESG = 'Use only full SCC matches'
            FULLSCC = ENVYN ( 'FULLSCC_ONLY', MESG, .FALSE., I )

C.............  Get error and warning limits from the environment

            MXERR  = ENVINT( ERRSET , ' ', 100, I )

            ALLOCATE( MTHPROF( NSRC,  NIPPA ),
     &                WEKPROF( NSRC,  NIPPA ),
     &                DOMPROF( NSRC,  NIPPA ),
     &                HRLPROF( NSRC,7,NIPPA ), STAT=ISTAT )

            IF ( ISTAT .NE. 0 ) THEN
                WRITE( MESG, '( A, I10 )' )
     &             'ERROR:  allocation failure.  STAT=', ISTAT
                CALL M3EXIT( PNAME, 0,0, MESG, 2 )
            END IF

            MTHPROF = IMISS3
            WEKPROF = IMISS3
            DOMPROF = IMISS3
            HRLPROF = IMISS3

            FIRSTIME = .FALSE.

        ENDIF

        CALL M3MSG2( 'Assigning temporal profiles to sources...' )

        ERRCNT = 0

        DO V = 0, NIPPA                 !  loop on pollutants: construct ASCII code
            WRITE( CPOS(V), '(I5.5)' ) V
        END DO

C.........  [Met-based profiles are managed in PROCTPRO()]

C.........  Set category-specific source characteristic combinations
C.........  Recall that profile-zero is the profile for the "none exist" case

        SELECT CASE ( CATEGORY )

            CASE ( 'AREA' )   ! Already set above

                DO S = 1, NSRC              !  loop on area sources

                    CSRC    = CSOURC( S )
                    TSCC    = CSCC( S )
                    TSCCL   = TSCC( 1:LSCCEND ) // ZEROS

                    CFIP    = CSRC( 1:FIPLEN3 )
                    IF ( USEEXPGEO() ) THEN
                        CFIPL = CFIP
                    ELSE
                        CFIPL = CFIP( 1:STALEN3 ) // ZEROS
                    END IF

C.....................  First, set all pollutants for pollutant independent part of hierarchy:

                    CALL BLDCSRC( CFIP,  TSCC,  BLANK, BLANK, BLANK, BLANK, BLANK, CPOS(0), CSRC07 )
                    CALL BLDCSRC( CFIPL, TSCC,  BLANK, BLANK, BLANK, BLANK, BLANK, CPOS(0), CSRC09 )
                    CALL BLDCSRC( CFIPZ, TSCC,  BLANK, BLANK, BLANK, BLANK, BLANK, CPOS(0), CSRC11 )
                    CALL BLDCSRC( CFIP,  TSCCZ, BLANK, BLANK, BLANK, BLANK, BLANK, CPOS(0), CSRC13 )
                    CALL BLDCSRC( CFIPL, TSCCZ, BLANK, BLANK, BLANK, BLANK, BLANK, CPOS(0), CSRC14 )
                    CALL BLDCSRC( CFIPZ, TSCCZ, BLANK, BLANK, BLANK, BLANK, BLANK, CPOS(0), CSRC15 )    !  ultimate fallback

                    IF ( .NOT.FULLSCC .AND. .NOT.CHKEXPSCC( TSCC ) ) THEN
                        CALL BLDCSRC( CFIP,  TSCCL, BLANK, BLANK, BLANK, BLANK, BLANK, CPOS(0), CSRC08 )
                        CALL BLDCSRC( CFIPL, TSCCL, BLANK, BLANK, BLANK, BLANK, BLANK, CPOS(0), CSRC10 )
                        CALL BLDCSRC( CFIPZ, TSCCL, BLANK, BLANK, BLANK, BLANK, BLANK, CPOS(0), CSRC12 )
                    END IF

C.....................  Find month-of-year profile:

                    IPROF =                     FINDC( CSRC07, MTHCOUNT,  MTHKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, MTHCOUNT,  MTHKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC09, MTHCOUNT,  MTHKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC10, MTHCOUNT,  MTHKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, MTHCOUNT,  MTHKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC12, MTHCOUNT,  MTHKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, MTHCOUNT,  MTHKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, MTHCOUNT,  MTHKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, MTHCOUNT,  MTHKEYS )
                    MTHPROF( S,: ) = MTHPDEX( MAX( IPROF,0 ) )

C.....................  Find day-of-month profile:

                    IPROF =                     FINDC( CSRC07, DOMCOUNT,  DOMKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, DOMCOUNT,  DOMKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC09, DOMCOUNT,  DOMKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC10, DOMCOUNT,  DOMKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, DOMCOUNT,  DOMKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC12, DOMCOUNT,  DOMKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, DOMCOUNT,  DOMKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, DOMCOUNT,  DOMKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, DOMCOUNT,  DOMKEYS )
                    DOMPROF( S,: ) = DOMPDEX( MAX( IPROF,0 ) )

C.....................  Find day-of-week profile:

                    IPROF =                     FINDC( CSRC07, WEKCOUNT,  WEKKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, WEKCOUNT,  WEKKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC09, WEKCOUNT,  WEKKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC10, WEKCOUNT,  WEKKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, WEKCOUNT,  WEKKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC12, WEKCOUNT,  WEKKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, WEKCOUNT,  WEKKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, WEKCOUNT,  WEKKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, WEKCOUNT,  WEKKEYS )
                    WEKPROF( S,: ) = WEKPDEX( MAX( IPROF,0 ) )

C.....................  Find hour-of-day profile for each day of the week:

                    IPROF =                     FINDC( CSRC07, MONCOUNT,  MONKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, MONCOUNT,  MONKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC09, MONCOUNT,  MONKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC10, MONCOUNT,  MONKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, MONCOUNT,  MONKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC12, MONCOUNT,  MONKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, MONCOUNT,  MONKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, MONCOUNT,  MONKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, MONCOUNT,  MONKEYS )
                    HRLPROF( S,1,: ) = MONPDEX( MAX( IPROF,0 ) )

                    IPROF =                     FINDC( CSRC07, TUECOUNT,  TUEKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, TUECOUNT,  TUEKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC09, TUECOUNT,  TUEKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC10, TUECOUNT,  TUEKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, TUECOUNT,  TUEKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC12, TUECOUNT,  TUEKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, TUECOUNT,  TUEKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, TUECOUNT,  TUEKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, TUECOUNT,  TUEKEYS )
                    HRLPROF( S,2,: ) = TUEPDEX( MAX( IPROF,0 ) )

                    IPROF =                     FINDC( CSRC07, WEDCOUNT,  WEDKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, WEDCOUNT,  WEDKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC09, WEDCOUNT,  WEDKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC10, WEDCOUNT,  WEDKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, WEDCOUNT,  WEDKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC12, WEDCOUNT,  WEDKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, WEDCOUNT,  WEDKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, WEDCOUNT,  WEDKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, WEDCOUNT,  WEDKEYS )
                    HRLPROF( S,3,: ) = WEDPDEX( MAX( IPROF,0 ) )

                    IPROF =                     FINDC( CSRC07, THUCOUNT,  THUKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, THUCOUNT,  THUKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC09, THUCOUNT,  THUKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC10, THUCOUNT,  THUKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, THUCOUNT,  THUKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC12, THUCOUNT,  THUKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, THUCOUNT,  THUKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, THUCOUNT,  THUKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, THUCOUNT,  THUKEYS )
                    HRLPROF( S,4,: ) = THUPDEX( MAX( IPROF,0 ) )

                    IPROF =                     FINDC( CSRC07, FRICOUNT,  FRIKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, FRICOUNT,  FRIKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC09, FRICOUNT,  FRIKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC10, FRICOUNT,  FRIKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, FRICOUNT,  FRIKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC12, FRICOUNT,  FRIKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, FRICOUNT,  FRIKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, FRICOUNT,  FRIKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, FRICOUNT,  FRIKEYS )
                    HRLPROF( S,5,: ) = FRIPDEX( MAX( IPROF,0 ) )

                    IPROF =                     FINDC( CSRC07, SATCOUNT,  SATKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, SATCOUNT,  SATKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC09, SATCOUNT,  SATKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC10, SATCOUNT,  SATKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, SATCOUNT,  SATKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC12, SATCOUNT,  SATKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, SATCOUNT,  SATKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, SATCOUNT,  SATKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, SATCOUNT,  SATKEYS )
                    HRLPROF( S,6,: ) = SATPDEX( MAX( IPROF,0 ) )

                    IPROF =                     FINDC( CSRC07, SUNCOUNT,  SUNKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, SUNCOUNT,  SUNKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC09, SUNCOUNT,  SUNKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC10, SUNCOUNT,  SUNKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, SUNCOUNT,  SUNKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC12, SUNCOUNT,  SUNKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, SUNCOUNT,  SUNKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, SUNCOUNT,  SUNKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, SUNCOUNT,  SUNKEYS )
                    HRLPROF( S,7,: ) = SUNPDEX( MAX( IPROF,0 ) )

C.....................  Now, overrides for pollutant dependent part of hierarchy:

                    DO V = 1, NIPPA         !  loop on pollutants

                        IF ( .NOT.POLREFFLAG( V ) )  CYCLE

                        CALL BLDCSRC( CFIP,  TSCC, BLANK, BLANK, BLANK, BLANK, BLANK, CPOS(V), CSRC01 )
                        CALL BLDCSRC( CFIPL, TSCC, BLANK, BLANK, BLANK, BLANK, BLANK, CPOS(V), CSRC03 )
                        CALL BLDCSRC( CFIPZ, TSCC, BLANK, BLANK, BLANK, BLANK, BLANK, CPOS(V), CSRC05 )

                        IF ( .NOT.FULLSCC .AND. .NOT.CHKEXPSCC( TSCC ) ) THEN
                            CALL BLDCSRC( CFIP,  TSCCL, BLANK, BLANK, BLANK, BLANK, BLANK, CPOS(V), CSRC02 )
                            CALL BLDCSRC( CFIPL, TSCCL, BLANK, BLANK, BLANK, BLANK, BLANK, CPOS(V), CSRC04 )
                            CALL BLDCSRC( CFIPZ, TSCCL, BLANK, BLANK, BLANK, BLANK, BLANK, CPOS(V), CSRC06 )
                        END IF

C.........................  Find month-of-year profile:

                        IPROF =                     FINDC( CSRC01, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC03, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC04, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .GT. 0 ) MTHPROF( S,V ) = MTHPDEX( IPROF )

C.........................  Find day-of-month profile:

                        IPROF =                     FINDC( CSRC01, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC03, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC04, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .GT. 0 ) DOMPROF( S,V ) = DOMPDEX( IPROF )

C.........................  Find day-of-week profile:

                        IPROF =                     FINDC( CSRC01, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC03, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC04, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .GT. 0 ) WEKPROF( S,V ) = WEKPDEX( IPROF )

C.........................  Find hour-of-day profile for each day of the week:

                        IPROF =                     FINDC( CSRC01, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC03, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC04, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, MONCOUNT,  MONKEYS )
                        IF ( IPROF .GT. 0 ) HRLPROF( S,1,V ) = MONPDEX( IPROF )

                        IPROF =                     FINDC( CSRC01, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC03, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC04, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .GT. 0 ) HRLPROF( S,2,V ) = TUEPDEX( IPROF )

                        IPROF =                     FINDC( CSRC01, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC03, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC04, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .GT. 0 ) HRLPROF( S,3,V ) = WEDPDEX( IPROF )

                        IPROF =                     FINDC( CSRC01, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC03, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC04, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, THUCOUNT,  THUKEYS )
                        IF ( IPROF .GT. 0 ) HRLPROF( S,4,V ) = THUPDEX( IPROF )

                        IPROF =                     FINDC( CSRC01, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC03, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC04, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .GT. 0 ) HRLPROF( S,5,V ) = FRIPDEX( IPROF )

                        IPROF =                     FINDC( CSRC01, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC03, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC04, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, SATCOUNT,  SATKEYS )
                        IF ( IPROF .GT. 0 ) HRLPROF( S,6,V ) = SATPDEX( IPROF )

                        IPROF =                     FINDC( CSRC01, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC03, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC04, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .GT. 0 ) HRLPROF( S,7,V ) = SUNPDEX( IPROF )

                    END DO                  !  end loop on pollutants

                END DO                      !  end loop- on area sources

            CASE ( 'MOBILE' )

                DO S = 1, NSRC              !  loop on mobile sources

                    CSRC    = CSOURC( S )
                    TSCC    = CSCC( S )
                    TSCCL   = TSCC( 1:LSCCEND ) // ZEROS

                    CFIP    = CSRC( 1:FIPLEN3 )
                    IF ( USEEXPGEO() ) THEN
                        CFIPL = CFIP
                    ELSE
                        CFIPL = CFIP( 1:STALEN3 ) // ZEROS
                    END IF

C.....................  First, set all pollutants for pollutant independent part of hierarchy:

                    CALL BLDCSRC( CFIP,  TSCC,  BLNK, BLNK, BLNK, BLNK, BLNK, CPOS(0), CSRC10 )
                    CALL BLDCSRC( CFIPL, TSCC,  BLNK, BLNK, BLNK, BLNK, BLNK, CPOS(0), CSRC13 )
                    CALL BLDCSRC( CFIPZ, TSCC,  BLNK, BLNK, BLNK, BLNK, BLNK, CPOS(0), CSRC15 )
                    CALL BLDCSRC( CFIP,  TSCCZ, BLNK, BLNK, BLNK, BLNK, BLNK, CPOS(0), CSRC19 )
                    CALL BLDCSRC( CFIPL, TSCCZ, BLNK, BLNK, BLNK, BLNK, BLNK, CPOS(0), CSRC20 )
                    CALL BLDCSRC( CFIPZ, TSCCZ, BLNK, BLNK, BLNK, BLNK, BLNK, CPOS(0), CSRC21 )    !  ultimate fallback

                    IF ( .NOT. FULLSCC .AND. .NOT.CHKEXPSCC( TSCC ) ) THEN
                        CALL BLDCSRC( CFIP,  TSCCL, BLNK, BLNK, BLNK, BLNK, BLNK, CPOS(0), CSRC11 )
                        CALL BLDCSRC( CFIPL, TSCCL, BLNK, BLNK, BLNK, BLNK, BLNK, CPOS(0), CSRC14 )
                        CALL BLDCSRC( CFIPZ, TSCCL, BLNK, BLNK, BLNK, BLNK, BLNK, CPOS(0), CSRC16 )
                    END IF

C.....................  Find month-of-year profile:

                    IPROF =                     FINDC( CSRC10, MTHCOUNT,  MTHKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, MTHCOUNT,  MTHKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, MTHCOUNT,  MTHKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, MTHCOUNT,  MTHKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, MTHCOUNT,  MTHKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC16, MTHCOUNT,  MTHKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC19, MTHCOUNT,  MTHKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC20, MTHCOUNT,  MTHKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC21, MTHCOUNT,  MTHKEYS )
                    MTHPROF( S,: ) = MTHPDEX( MAX( IPROF,0 ) )

C.....................  Find day-of-month profile:

                    IPROF =                     FINDC( CSRC10, DOMCOUNT,  DOMKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, DOMCOUNT,  DOMKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, DOMCOUNT,  DOMKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, DOMCOUNT,  DOMKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, DOMCOUNT,  DOMKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC16, DOMCOUNT,  DOMKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC19, DOMCOUNT,  DOMKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC20, DOMCOUNT,  DOMKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC21, DOMCOUNT,  DOMKEYS )
                    DOMPROF( S,: ) = DOMPDEX( MAX( IPROF,0 ) )

C.....................  Find day-of-week profile:

                    IPROF =                     FINDC( CSRC10, WEKCOUNT,  WEKKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, WEKCOUNT,  WEKKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, WEKCOUNT,  WEKKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, WEKCOUNT,  WEKKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, WEKCOUNT,  WEKKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC16, WEKCOUNT,  WEKKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC19, WEKCOUNT,  WEKKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC20, WEKCOUNT,  WEKKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC21, WEKCOUNT,  WEKKEYS )
                    WEKPROF( S,: ) = WEKPDEX( MAX( IPROF,0 ) )

C.....................  Find hour-of-day profile for each day of the week:

                    IPROF =                     FINDC( CSRC10, MONCOUNT,  MONKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, MONCOUNT,  MONKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, MONCOUNT,  MONKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, MONCOUNT,  MONKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, MONCOUNT,  MONKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC16, MONCOUNT,  MONKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC19, MONCOUNT,  MONKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC20, MONCOUNT,  MONKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC21, MONCOUNT,  MONKEYS )
                    HRLPROF( S,1,: ) = MONPDEX( MAX( IPROF,0 ) )

                    IPROF =                     FINDC( CSRC10, TUECOUNT,  TUEKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, TUECOUNT,  TUEKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, TUECOUNT,  TUEKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, TUECOUNT,  TUEKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, TUECOUNT,  TUEKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC16, TUECOUNT,  TUEKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC19, TUECOUNT,  TUEKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC20, TUECOUNT,  TUEKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC21, TUECOUNT,  TUEKEYS )
                    HRLPROF( S,2,: ) = TUEPDEX( MAX( IPROF,0 ) )

                    IPROF =                     FINDC( CSRC10, WEDCOUNT,  WEDKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, WEDCOUNT,  WEDKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, WEDCOUNT,  WEDKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, WEDCOUNT,  WEDKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, WEDCOUNT,  WEDKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC16, WEDCOUNT,  WEDKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC19, WEDCOUNT,  WEDKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC20, WEDCOUNT,  WEDKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC21, WEDCOUNT,  WEDKEYS )
                    HRLPROF( S,3,: ) = WEDPDEX( MAX( IPROF,0 ) )

                    IPROF =                     FINDC( CSRC10, THUCOUNT,  THUKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, THUCOUNT,  THUKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, THUCOUNT,  THUKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, THUCOUNT,  THUKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, THUCOUNT,  THUKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC16, THUCOUNT,  THUKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC19, THUCOUNT,  THUKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC20, THUCOUNT,  THUKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC21, THUCOUNT,  THUKEYS )
                    HRLPROF( S,4,: ) = THUPDEX( MAX( IPROF,0 ) )

                    IPROF =                     FINDC( CSRC10, FRICOUNT,  FRIKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, FRICOUNT,  FRIKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, FRICOUNT,  FRIKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, FRICOUNT,  FRIKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, FRICOUNT,  FRIKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC16, FRICOUNT,  FRIKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC19, FRICOUNT,  FRIKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC20, FRICOUNT,  FRIKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC21, FRICOUNT,  FRIKEYS )
                    HRLPROF( S,5,: ) = FRIPDEX( MAX( IPROF,0 ) )

                    IPROF =                     FINDC( CSRC10, SATCOUNT,  SATKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, SATCOUNT,  SATKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, SATCOUNT,  SATKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, SATCOUNT,  SATKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, SATCOUNT,  SATKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC16, SATCOUNT,  SATKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC19, SATCOUNT,  SATKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC20, SATCOUNT,  SATKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC21, SATCOUNT,  SATKEYS )
                    HRLPROF( S,6,: ) = SATPDEX( MAX( IPROF,0 ) )

                    IPROF =                     FINDC( CSRC10, SUNCOUNT,  SUNKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, SUNCOUNT,  SUNKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, SUNCOUNT,  SUNKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, SUNCOUNT,  SUNKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, SUNCOUNT,  SUNKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC16, SUNCOUNT,  SUNKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC19, SUNCOUNT,  SUNKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC20, SUNCOUNT,  SUNKEYS )
                    IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC21, SUNCOUNT,  SUNKEYS )
                    HRLPROF( S,7,: ) = SUNPDEX( MAX( IPROF,0 ) )

C.....................  Now, overrides for pollutant dependent part of hierarchy:

                    DO V = 1, NIPPA         !  loop on pollutants

                        IF ( .NOT.POLREFFLAG( V ) )  CYCLE

                        CALL BLDCSRC( CFIP,  TSCC,  BLNK, BLNK, BLNK, BLNK, BLNK, CPOS(V), CSRC01 )
                        CALL BLDCSRC( CFIPL, TSCC,  BLNK, BLNK, BLNK, BLNK, BLNK, CPOS(V), CSRC05 )
                        CALL BLDCSRC( CFIPZ, TSCC,  BLNK, BLNK, BLNK, BLNK, BLNK, CPOS(V), CSRC07 )

                        IF ( .NOT. FULLSCC .AND. .NOT.CHKEXPSCC( TSCC ) ) THEN
                            CALL BLDCSRC( CFIP,  TSCCL, BLNK, BLNK, BLNK, BLNK, BLNK, CPOS(V), CSRC02 )
                            CALL BLDCSRC( CFIPL, TSCCL, BLNK, BLNK, BLNK, BLNK, BLNK, CPOS(V), CSRC06 )
                            CALL BLDCSRC( CFIPZ, TSCCL, BLNK, BLNK, BLNK, BLNK, BLNK, CPOS(V), CSRC08 )
                        END IF

C.........................  Find month-of-year profile:

                        IPROF =                     FINDC( CSRC01, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC07, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .GT. 0 ) MTHPROF( S,V ) = MTHPDEX( IPROF )

C.........................  Find day-of-month profile:

                        IPROF =                     FINDC( CSRC01, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC07, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .GT. 0 ) DOMPROF( S,V ) = DOMPDEX( IPROF )

C.........................  Find day-of-week profile:

                        IPROF =                     FINDC( CSRC01, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC07, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .GT. 0 ) WEKPROF( S,V ) = WEKPDEX( IPROF )

C.........................  Find hour-of-day profile for each day of the week:

                        IPROF =                     FINDC( CSRC01, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC07, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, MONCOUNT,  MONKEYS )
                        IF ( IPROF .GT. 0 ) HRLPROF( S,1,V ) = MONPDEX( IPROF )

                        IPROF =                     FINDC( CSRC01, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC07, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .GT. 0 ) HRLPROF( S,2,V ) = TUEPDEX( IPROF )

                        IPROF =                     FINDC( CSRC01, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC07, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .GT. 0 ) HRLPROF( S,3,V ) = WEDPDEX( IPROF )

                        IPROF =                     FINDC( CSRC01, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC07, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, THUCOUNT,  THUKEYS )
                        IF ( IPROF .GT. 0 ) HRLPROF( S,4,V ) = THUPDEX( IPROF )

                        IPROF =                     FINDC( CSRC01, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC07, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .GT. 0 ) HRLPROF( S,5,V ) = FRIPDEX( IPROF )

                        IPROF =                     FINDC( CSRC01, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC07, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, SATCOUNT,  SATKEYS )
                        IF ( IPROF .GT. 0 ) HRLPROF( S,6,V ) = SATPDEX( IPROF )

                        IPROF =                     FINDC( CSRC01, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC07, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .GT. 0 ) HRLPROF( S,7,V ) = SUNPDEX( IPROF )

                    END DO                  !  end loop on pollutants

                END DO                      !  end loop on mobile sources

            CASE ( 'POINT' )

                DO S = 1, NSRC              !  loop on point sources

                    CSRC    = CSOURC( S )
                    TSCC    = CSCC( S )
                    TSCC5   = TSCC( 1:SCCEXPLEN3+5 ) // ZEROS

                    CFIP    = CSRC( 1:FIPLEN3 )
                    IF ( USEEXPGEO() ) THEN
                        CFIPL = CFIP
                    ELSE
                        CFIPL = CFIP( 1:STALEN3 ) // ZEROS
                    END IF

                    CPLT = ADJUSTL( CSRC( PTBEGL3(2):PTENDL3(2) ) )
                    CPNT = ADJUSTL( CSRC( PTBEGL3(3):PTENDL3(3) ) )
                    CSTK = ADJUSTL( CSRC( PTBEGL3(4):PTENDL3(4) ) )
                    CSEG = ADJUSTL( CSRC( PTBEGL3(5):PTENDL3(5) ) )
                    CPL5 = ADJUSTL( CSRC( PTBEGL3(6):PTENDL3(6) ) )
                    CPLL = TSCC5                          ! zero scc as a plant-characteristic
                    CPLZ = TSCCZ                          ! zero scc as a plant-characteristic

C.....................  pollutant-independent search targets:

                    CALL BLDCSRC( CFIP,  CPLT,    CPNT, CSTK, CSEG, CPL5, TSCC,  CPOS(0), CSRC05 )
                    CALL BLDCSRC( CFIP,  CPLT,    CPNT, CSTK, BLNK, CPL5, TSCC,  CPOS(0), CSRC06 )
                    CALL BLDCSRC( CFIP,  CPLT,    CPNT, BLNK, BLNK, CPL5, TSCC,  CPOS(0), CSRC07 )
                    CALL BLDCSRC( CFIP,  CPLT,    BLNK, BLNK, BLNK, CPL5, TSCC,  CPOS(0), CSRC08 )
                    CALL BLDCSRC( CFIP,  CPLT,    CPNT, BLNK, BLNK, CPLZ, TSCCZ, CPOS(0), CSRC084 )
                    CALL BLDCSRC( CFIP,  CPLT,    BLNK, BLNK, BLNK, CPLZ, TSCCZ, CPOS(0), CSRC085 )
                    CALL BLDCSRC( CFIP,  BLNKPLT, BLNK, BLNK, BLNK, CPL5, TSCC,  CPOS(0), CSRC15 )
                    CALL BLDCSRC( CFIPL, BLNKPLT, BLNK, BLNK, BLNK, CPL5, TSCC,  CPOS(0), CSRC17 )
                    CALL BLDCSRC( CFIPZ, BLNKPLT, BLNK, BLNK, BLNK, CPL5, TSCC,  CPOS(0), CSRC19 )
                    CALL BLDCSRC( CFIP,  BLNKPLT, BLNK, BLNK, BLNK, CPLZ, TSCCZ, CPOS(0), CSRC21 )
                    CALL BLDCSRC( CFIPL, BLNKPLT, BLNK, BLNK, BLNK, CPLZ, TSCCZ, CPOS(0), CSRC22 )
                    CALL BLDCSRC( CFIPZ, BLNKPLT, BLNK, BLNK, BLNK, CPLZ, TSCCZ, CPOS(0), CSRC23 )    !  ultimate fallback

                    IF ( .NOT. FULLSCC .AND. .NOT.CHKEXPSCC( TSCC ) ) THEN
                        CALL BLDCSRC( CFIP,  BLNKPLT, BLNK, BLNK, BLNK, CPLL, TSCC5, CPOS(0), CSRC16 )
                        CALL BLDCSRC( CFIPL, BLNKPLT, BLNK, BLNK, BLNK, CPLL, TSCC5, CPOS(0), CSRC18 )
                        CALL BLDCSRC( CFIPZ, BLNKPLT, BLNK, BLNK, BLNK, CPLL, TSCC5, CPOS(0), CSRC20 )
                    END IF

C.....................  Note that pollutant dependent and pollutant independent
C.....................  parts of the point source heirarchy search are tangled together:

                    DO V = 1, NIPPA         !  loop on pollutants

                        CALL BLDCSRC( CFIP,  CPLT,    CPNT, CSTK, CSEG, CPL5, TSCC, CPOS(V), CSRC01 )
                        CALL BLDCSRC( CFIP,  CPLT,    CPNT, CSTK, BLNK, CPL5, TSCC, CPOS(V), CSRC02 )
                        CALL BLDCSRC( CFIP,  CPLT,    CPNT, BLNK, BLNK, CPL5, TSCC, CPOS(V), CSRC03 )
                        CALL BLDCSRC( CFIP,  CPLT,    BLNK, BLNK, BLNK, CPL5, TSCC, CPOS(V), CSRC04 )
                        CALL BLDCSRC( CFIP,  BLNKPLT, BLNK, BLNK, BLNK, CPL5, TSCC, CPOS(V), CSRC09 )
                        CALL BLDCSRC( CFIPL, BLNKPLT, BLNK, BLNK, BLNK, CPL5, TSCC, CPOS(V), CSRC11 )
                        CALL BLDCSRC( CFIPZ, BLNKPLT, BLNK, BLNK, BLNK, CPL5, TSCC, CPOS(V), CSRC13 )

                        IF ( .NOT. FULLSCC .AND. .NOT.CHKEXPSCC( TSCC ) ) THEN
                            CALL BLDCSRC( CFIP,  BLNKPLT, BLNK, BLNK, BLNK, CPLL, TSCC5, CPOS(V), CSRC10 )
                            CALL BLDCSRC( CFIPL, BLNKPLT, BLNK, BLNK, BLNK, CPLL, TSCC5, CPOS(V), CSRC12 )
                            CALL BLDCSRC( CFIPZ, BLNKPLT, BLNK, BLNK, BLNK, CPLL, TSCC5, CPOS(V), CSRC14 )
                        END IF

C.........................  Find month-of-year profile:

                        IPROF =                     FINDC( CSRC01, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC03, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC04, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC07, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC084,MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC085,MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC09, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC10, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC12, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC16, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC17, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC18, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC19, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC20, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC21, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC22, MTHCOUNT,  MTHKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC23, MTHCOUNT,  MTHKEYS )
                        MTHPROF( S,V ) = MTHPDEX( MAX( IPROF,0 ) )

C.........................  Find day-of-month profile:

                        IPROF =                     FINDC( CSRC01, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC03, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC04, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC07, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC084,DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC085,DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC09, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC10, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC12, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC16, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC17, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC18, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC19, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC20, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC21, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC22, DOMCOUNT,  DOMKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC23, DOMCOUNT,  DOMKEYS )
                        DOMPROF( S,V ) = DOMPDEX( MAX( IPROF,0 ) )

C.........................  Find day-of-week profile:

                        IPROF =                     FINDC( CSRC01, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC03, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC04, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC07, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC084,WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC085,WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC09, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC10, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC12, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC16, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC17, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC18, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC19, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC20, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC21, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC22, WEKCOUNT,  WEKKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC23, WEKCOUNT,  WEKKEYS )
                        WEKPROF( S,V ) = WEKPDEX( MAX( IPROF,0 ) )

C.........................  Find hour-of-day profiles for each day of the week:

                        IPROF =                     FINDC( CSRC01, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC03, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC04, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC07, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC084,MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC085,MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC09, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC10, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC12, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC16, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC17, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC18, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC19, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC20, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC21, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC22, MONCOUNT,  MONKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC23, MONCOUNT,  MONKEYS )
                        HRLPROF( S,1,V ) = MONPDEX( MAX( IPROF,0 ) )

                        IPROF =                     FINDC( CSRC01, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC03, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC04, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC07, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC084,TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC085,TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC09, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC10, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC12, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC16, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC17, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC18, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC19, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC20, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC21, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC22, TUECOUNT,  TUEKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC23, TUECOUNT,  TUEKEYS )
                        HRLPROF( S,2,V ) = TUEPDEX( MAX( IPROF,0 ) )

                        IPROF =                     FINDC( CSRC01, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC03, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC04, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC07, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC084,WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC085,WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC09, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC10, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC12, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC16, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC17, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC18, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC19, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC20, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC21, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC22, WEDCOUNT,  WEDKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC23, WEDCOUNT,  WEDKEYS )
                        HRLPROF( S,3,V ) = WEDPDEX( MAX( IPROF,0 ) )

                        IPROF =                     FINDC( CSRC01, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC03, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC04, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC07, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC084,THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC085,THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC09, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC10, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC12, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC16, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC17, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC18, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC19, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC20, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC21, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC22, THUCOUNT,  THUKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC23, THUCOUNT,  THUKEYS )
                        HRLPROF( S,4,V ) = THUPDEX( MAX( IPROF,0 ) )

                        IPROF =                     FINDC( CSRC01, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC03, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC04, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC07, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC084,FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC085,FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC09, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC10, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC12, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC16, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC17, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC18, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC19, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC20, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC21, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC22, FRICOUNT,  FRIKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC23, FRICOUNT,  FRIKEYS )
                        HRLPROF( S,5,V ) = FRIPDEX( MAX( IPROF,0 ) )

                        IPROF =                     FINDC( CSRC01, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC03, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC04, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC07, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC084,SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC085,SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC09, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC10, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC12, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC16, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC17, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC18, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC19, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC20, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC21, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC22, SATCOUNT,  SATKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC23, SATCOUNT,  SATKEYS )
                        HRLPROF( S,6,V ) = SATPDEX( MAX( IPROF,0 ) )

                        IPROF =                     FINDC( CSRC01, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC02, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC03, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC04, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC05, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC06, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC07, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC08, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC084,SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC085,SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC09, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC10, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC11, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC12, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC13, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC14, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC15, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC16, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC17, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC18, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC19, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC20, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC21, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC22, SUNCOUNT,  SUNKEYS )
                        IF ( IPROF .LE. 0 ) IPROF = FINDC( CSRC23, SUNCOUNT,  SUNKEYS )
                        HRLPROF( S,7,V ) = SUNPDEX( MAX( IPROF,0 ) )

                    END DO                  !  end loop on pollutants

                END DO                      !  end loop on point sources

            CASE DEFAULT

                MESG = 'ERROR:  unrecognized category "' // TRIM( CATEGORY ) // '"'
                CALL M3EXIT( PNAME, 0,0, MESG, 2 )

        END SELECT

C.........  Check that all sources, pollutants have assigned profiles:

        IF ( MTHCOUNT .EQ. 0 ) THEN
            CALL M3MESG( 'ASGNTPRO:  Not using Month-of-Year profiles' )
        END IF

        IF ( WEKCOUNT + DOMCOUNT .EQ. 0 ) THEN
            CALL M3MESG( 'ASGNTPRO:  Not using Day-of-Month nor Day-of-Week profiles' )
        END IF

        HCOUNT = METCOUNT + MONCOUNT + TUECOUNT + WEDCOUNT + THUCOUNT + FRICOUNT + SATCOUNT + SUNCOUNT
        IF ( HCOUNT .EQ. 0 ) THEN
            CALL M3MESG( 'ASGNTPRO:  Not using Hour-of-Day nor Met Based profiles' )
        END IF

        DO S = 1, NSRC

            IF ( MOD( TPFLAG(S), MTPRFAC ) .EQ. 0 ) THEN

                ICOUNT = 0
                DO V = 1, NIPPA
                    IF ( MTHPROF( S,V ) .LE. 0 ) THEN
                        ICOUNT = ICOUNT + 1
                    END IF
                END DO          !  end loop on pollutants, V

                IF ( ICOUNT .EQ. NIPPA ) THEN
                    ERRCNT = ERRCNT + 1
                    CALL FMTCSRC( CSOURC(S), NCHARS, CBUF, L2 )
                    WRITE( MESG, '( A, I8, 3( 1X, A ) )' )
     &                    'ERROR:  No month-of-year profile found for source', S,
     &                    ':  ', TRIM( CBUF )
                        CALL M3MESG( MESG )
                ELSE IF ( ICOUNT .GT. 0 ) THEN
                    ERRCNT = ERRCNT + 1
                    CALL FMTCSRC( CSOURC(S), NCHARS, CBUF, L2 )
                    DO V = 1, NIPPA
                        IF ( MTHPROF( S,V ) .LE. 0 ) THEN
                            WRITE( MESG, '( A, I8, 4( 1X, A ) )' )
     &                        'ERROR:  No month-of-year profile found for source', S,
     &                        'pollutant', TRIM( EANAM(V) ), 'source', TRIM( CBUF )
                            CALL M3MESG( MESG )
                        END IF
                    END DO
                END IF

            END IF      !  if monthly profiles needed


            IF ( MOD( TPFLAG(S), WDTPFAC ) .EQ. 0 .OR.
     &           MOD( TPFLAG(S), WTPRFAC ) .EQ. 0 ) THEN

                ICOUNT = 0
                DO V = 1, NIPPA
                    IF ( ( WEKPROF( S,V ) .LE. 0 ) .AND.
     &                   ( DOMPROF( S,V ) .LE. 0 ) ) THEN
                        ICOUNT = ICOUNT + 1
                    END IF
                END DO          !  end loop on pollutants, V

                IF ( ICOUNT .EQ. NIPPA ) THEN
                    ERRCNT = ERRCNT + 1
                    CALL FMTCSRC( CSOURC(S), NCHARS, CBUF, L2 )
                    WRITE( MESG, '( A, I8, 3( 1X, A ) )' )
     &                  'ERROR:  No day-of-month nor day-of-week for source', S,
     &                  ':  ', TRIM( CBUF )
                        CALL M3MESG( MESG )
                ELSE IF ( ICOUNT .GT. 0 ) THEN
                    ERRCNT = ERRCNT + 1
                    CALL FMTCSRC( CSOURC(S), NCHARS, CBUF, L2 )
                    DO V = 1, NIPPA
                        IF ( WEKPROF( S,V ) .LE. 0.AND.
     &                       DOMPROF( S,V ) .LE. 0 ) THEN
                            WRITE( MESG, '( A, I8, 4( 1X, A ) )' )
     &                        'ERROR:  No day-of-month nor day-of-week profile found for source', S,
     &                        'pollutant', TRIM( EANAM(V) ), 'source', TRIM( CBUF )
                            CALL M3MESG( MESG )
                        END IF
                    END DO
                END IF

            END IF      !  if weekly profiles needed

            DO I = 1, 7

                IF ( DAYFLAG(I)  ) THEN

                    ICOUNT = 0
                    DO V = 1, NIPPA
                        IF ( HRLPROF( S,I,V ) .LE. 0  .AND.
     &                       METPROF( S,V )   .LE. 0 ) THEN
                            ICOUNT = ICOUNT + 1
                        END IF
                    END DO          !  end loop on pollutants, V

                    IF ( HCOUNT .EQ. 0 ) THEN
                        CONTINUE
                    ELSE IF ( ICOUNT .EQ. NIPPA ) THEN
                        ERRCNT = ERRCNT + 1
                        CALL FMTCSRC( CSOURC(S), NCHARS, CBUF, L2 )
                        WRITE( MESG, '( A, I8, 4( 1X, A ) )' )
     &                        'ERROR:  No hour-of-day nor Met based profile profile found for source', S,
     &                        ':  ', TRIM( CBUF ), 'day', DAYNAME(I)
                            CALL M3MESG( MESG )
                    ELSE IF ( ICOUNT .GT. 0 ) THEN
                        ERRCNT = ERRCNT + 1
                        CALL FMTCSRC( CSOURC(S), NCHARS, CBUF, L2 )
                        DO V = 1, NIPPA
                            IF ( HRLPROF( S,I,V ) .LE. 0  .AND.
     &                           METPROF( S,V )   .LE. 0 ) THEN
                                WRITE( MESG, '( A, I8, 6( 1X, A ) )' )
     &                            'ERROR:  No hour-of-day nor Met based profile profile found for source', S,
     &                            'pollutant', TRIM( EANAM(V) ), 'source', TRIM( CBUF ),
     &                            'day', DAYNAME(I)
                                CALL M3MESG( MESG )
                            END IF
                        END DO
                    END IF

                END IF      !  if dayflag

            END DO      !  end loop on days

            IF ( ERRCNT .GT. MXERR ) THEN
                CALL M3MESG( ' ' )
                CALL M3MESG( 'ASGNTPRO:  Maximum number of errors exceeded.' )
                EXIT
            END IF

        END DO      !  end loop on sources, S

        IF( ERRCNT .GT. 0 ) THEN
            MESG = 'Problem assigning temporal profiles to sources'
            CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
        END IF

        DEALLOCATE( MTHPDEX,
     &              WEKPDEX,
     &              DOMPDEX,
     &              MONPDEX,
     &              TUEPDEX,
     &              WEDPDEX,
     &              THUPDEX,
     &              FRIPDEX,
     &              SATPDEX,
     &              SUNPDEX,
     &              MTHKEYS,
     &              WEKKEYS,
     &              DOMKEYS,
     &              MONKEYS,
     &              TUEKEYS,
     &              WEDKEYS,
     &              THUKEYS,
     &              FRIKEYS,
     &              SATKEYS,
     &              SUNKEYS )

        RETURN

        END SUBROUTINE ASGNTPRO
