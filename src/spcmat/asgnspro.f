
        SUBROUTINE ASGNSPRO( MASSOUT, MOLEOUT, REPORT, NSRCIN, UDEV,
     &                       ENAM, LVALCHK, MASSMATX, MOLEMATX )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      For each source and current pollutant, find the most specific speciation
C      profile that applies to that source. Do this using the grouped tables of
C      speciation cross references from RDSREF.  The hierarchical order is
C      defined in this subroutine, and can be determined from the in-source
C      comments below. Once a profile code has been identified, search for this
C      code in the speciation profile tables (from RDSPROF) and use this profile
C      to update the previously initialized speciation matrices.
C
C  PRECONDITIONS REQUIRED:
C      Expects cross-reference tables to be set to EMCMISS3 if not defined
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 2/99 by M. Houyoux
C
C***************************************************************************
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
C...........   This module contains the source arrays
        USE MODSOURC, ONLY: CSOURC, CSCC, IRCLAS, IVTYPE, CISIC, CMACT

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVIFIP, INVCFIP

C...........   This module contains the cross-reference tables
        USE MODXREF, ONLY: TXCNT, CHRT02, CHRT03, CHRT04,
     &          CHRT05, CHRT06, CHRT07, CHRT08, CHRT09, CHRT10,
     &          CHRT11, CHRT12, CHRT13, CHRT14, CHRT15, CHRT16,
     &          CHRT26, CHRT27, CHRT28, CHRT29, CHRT30, CHRT31,
     &          CHRT32, CHRT33, CHRT34, CHRT35, CHRT36, CHRT37,
     &          CSPT01, CSPT02, CSPT03, CSPT04,
     &          CSPT05, CSPT06, CSPT07, CSPT08, CSPT09, CSPT10,
     &          CSPT11, CSPT12, CSPT13, CSPT14, CSPT15, CSPT16,
     &          CSPT26, CSPT27, CSPT28, CSPT29, CSPT30, CSPT31,
     &          CSPT32, CSPT33, CSPT34, CSPT35, CSPT36, CSPT37

C...........   This module contains the speciation profile tables
        USE MODSPRO, ONLY: MXSPEC, CNVRT01, CNVRT02, CNVRT03, CNVRT04,
     &                     NCNV1, NCNV2, NCNV3, NCNV4, NSPROF, SPROFN,
     &                     CNVFC00, CNVFC01, CNVFC02, CNVFC03, CNVFC04,
     &                     IDXSPRO, IDXSSPEC, NSPECIES, CNVFLAG,
     &                     MASSFACT, MOLEFACT, CMBCNT, CMBNP, CMBSPCD,
     &                     CMBWGHT, CMBMAX, CMBDEX

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, NCHARS, JSCC, NIPPA, EANAM,
     &                     LSCCEND

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api constant parameters

C...........   PARAMETERs

        CHARACTER(1),  PARAMETER :: QUOTE    = "'"
        CHARACTER(16), PARAMETER :: PROGNAME = 'ASGNSPRO' ! program name

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)    CRLF
        LOGICAL         ENVYN
        INTEGER         FIND1
        INTEGER         FINDC
        INTEGER         INDEX1
        INTEGER         ENVINT
        INTEGER         PROMPTFFILE
        LOGICAL         SETSCCTYPE

        EXTERNAL CRLF, ENVYN, FINDC, FIND1, INDEX1, PROMPTFFILE,
     &           SETSCCTYPE, ENVINT

C.........  SUBROUTINE ARGUMENTS
        LOGICAL     , INTENT    (IN) :: MASSOUT        ! true: create mass-based
        LOGICAL     , INTENT    (IN) :: MOLEOUT        ! true: create mole-based
        LOGICAL     , INTENT    (IN) :: REPORT         ! true: rep defaults
        INTEGER     , INTENT    (IN) :: NSRCIN         ! number of sources
        INTEGER     , INTENT    (IN) :: UDEV           ! suplmt file unit no.
        CHARACTER(*), INTENT    (IN) :: ENAM    ! pol/emis type name of interest
        LOGICAL     , INTENT    (IN) :: LVALCHK ! true: check value of pol/emis from EI
        REAL        , INTENT(IN OUT) :: MASSMATX( NSRCIN,* )! mass spec matx
        REAL        , INTENT(IN OUT) :: MOLEMATX( NSRCIN,* )! mole spec matx

C.........  Local arrays
        REAL                 :: EMISTMP( NSRCIN )  ! tmp emis for 1 pollutant by source
        CHARACTER( SPNLEN3 ) :: CCODE  ( CMBMAX )
        REAL                 :: CWEIGHT( CMBMAX )
        INTEGER, SAVE        :: ERRCNT( 8 )
        INTEGER, SAVE        :: WARNCNT( 2 )

C.........  Other local variables
        INTEGER         K, L, L2, NP, S, V    !  counters and indices

        INTEGER          F0, F1, F2, F3, F4, F5, F6  ! tmp find indices
        INTEGER       :: F0B = 0      ! extra find index for mobile
        INTEGER       :: F2B = 0      ! extra find index for mobile
        INTEGER       :: F4B = 0      ! extra find index for mobile
        INTEGER          IOS          ! i/o status
        INTEGER          NCHKCHR      ! position of last non-SCC src char
        INTEGER          NCOUT        ! no. output source chars for mesgs
        INTEGER, SAVE :: MXWARN       ! maximum number of warnings of each type to write
        INTEGER, SAVE :: MXERR        ! maximum number of errors of each type to write

        LOGICAL       :: EFLAG    = .FALSE. ! true: error detected
        LOGICAL       :: LRDCOMBO           ! true: for first COMBO for a pollutant/emis type
        LOGICAL, SAVE :: CNVERROR = .TRUE.  ! true: when default (=1.0) conv fac is applied
        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time subrtn called
        LOGICAL, SAVE :: MACTFLAG = .FALSE. ! true: MACT codes available in inventory
        LOGICAL, SAVE :: REPDEFLT = .TRUE.  ! true: report when defaults used
        LOGICAL, SAVE :: LNOZERO  = .FALSE. ! true: skip speication for zero emission
        LOGICAL, SAVE :: SICFLAG  = .FALSE. ! true: SIC available in inventory
        LOGICAL          SCCFLAG            ! true: SCC type is different from previous

        CHARACTER(300)       BUFFER  ! source fields buffer
        CHARACTER(300)       MESG    ! message buffer
        CHARACTER(FIPLEN3)   CFIP    ! tmp (character) FIPS code
        CHARACTER(STALEN3)   CSTA    ! tmp Country/state code
        CHARACTER(SRCLEN3)   CSRC    ! tmp source chars string
        CHARACTER(SCCLEN3)   TSCC    ! tmp 10-digit SCC
        CHARACTER(SCCLEN3)   SCCORIG ! tmp original SCC for mobile processing
        CHARACTER(SCCLEN3)   TSCCL   ! tmp left digits of TSCC
        CHARACTER(SCCLEN3)   TSCCINIT! tmp initial 10-digit SCC
        CHARACTER(SPNLEN3)   SPCODE  ! tmp speciation profile code
        CHARACTER(SCCLEN3)   CHKRWT  ! tmp roadway type only SCC -- BUG:  UNINITIALIZED -- CJC
        CHARACTER(SS5LEN3):: CSRC5=' '! tmp source chars through char5
        CHARACTER(SS4LEN3):: CSRC4=' '! tmp source chars through char4
        CHARACTER(SS3LEN3):: CSRC3=' '! tmp source chars through char3
        CHARACTER(SS2LEN3):: CSRC2=' '! tmp source chars through char2
        CHARACTER(SS1LEN3):: CSRC1=' '! tmp source chars through char1
        CHARACTER(SS5LEN3):: CHK16=' '! tmp source chars through char5// SCC
        CHARACTER(SS4LEN3):: CHK15=' '! tmp source chars through char4// SCC
        CHARACTER(SS3LEN3):: CHK14=' '! tmp source chars through char3// SCC
        CHARACTER(SS2LEN3):: CHK13=' '! tmp source chars through char2// SCC
        CHARACTER(SS1LEN3):: CHK12=' '! tmp source chars through char1// SCC
        CHARACTER(SS0LEN3):: CHK11=' '! tmp FIPS // Plant // SCC
        CHARACTER(FPLLEN3):: CHK10=' '! tmp FIPS code // plant id
        CHARACTER(FPSLEN3):: CHK09=' '! tmp FIPS code // SCC
        CHARACTER(FPSLEN3):: CHK08=' '! tmp FIPS code // left SCC
        CHARACTER(FPSLEN3):: CHK08B=' '! tmp FIPS code // veh ID SCC
        CHARACTER(STSLEN3):: CHK06=' '! tmp Country/state code // SCC
        CHARACTER(STSLEN3):: CHK05=' '! tmp Country/state code // left SCC
        CHARACTER(STSLEN3):: CHK05B=' '! tmp Country/state code// veh ID SCC
        CHARACTER(SCCLEN3):: CHK02B=' '! tmp veh ID SCC
        CHARACTER(STILEN3):: CHK28=' '! tmp Country/state code // left SIC
        CHARACTER(STILEN3):: CHK29=' '! tmp Country/state code // SIC
        CHARACTER(FPILEN3):: CHK30=' '! tmp FIPS code // left SIC
        CHARACTER(FPILEN3):: CHK31=' '! tmp FIPS code // SIC
        CHARACTER(MSCLEN3):: CHK33=' '! tmp SCC // MACT
        CHARACTER(MSTLEN3):: CHK34=' '! tmp Country/state code // MACT
        CHARACTER(MSSLEN3):: CHK35=' '! tmp Country/state code // SCC // MACT
        CHARACTER(MFPLEN3):: CHK36=' '! tmp FIPS code // MACT
        CHARACTER(MFSLEN3):: CHK37=' '! tmp FIPS code // SCC // MACT
        CHARACTER(MACLEN3)   CMCT    ! tmp MACT code
        CHARACTER(SICLEN3)   CSIC    ! tmp SIC code
        CHARACTER(SICLEN3)   CSICL   ! tmp left SIC code

C***********************************************************************
C   begin body of subroutine ASGNSPRO

C.........  For first time routine is called in all cases,
        IF( FIRSTIME ) THEN

C.............  Get maximum number of warnings/errors

            MXWARN = ENVINT( WARNSET , ' ', 100, IOS )
            IF ( IOS .GT. 0 ) THEN
                MESG = 'Bad environment variable "' // WARNSET // '"'
                CALL M3EXIT( PROGNAME, 0,0, MESG, 2 )
            END IF

            MXERR  = ENVINT( ERRSET  , ' ', 100, IOS )
            IF ( IOS .GT. 0 ) THEN
                MESG = 'Bad environment variable "' // ERRSET // '"'
                CALL M3EXIT( PROGNAME, 0,0, MESG, 2 )
            END IF
            ERRCNT  = 0   ! array
            WARNCNT = 0   ! array

C.............  Retrieve environment variables
            MESG = 'Switch for reporting default speciation profiles'
            REPDEFLT = ENVYN ( 'REPORT_DEFAULTS', MESG, .TRUE., IOS )
            IF ( IOS .GT. 0 ) THEN
                MESG = 'Bad environment variable "REPORT_DEFAULTS"'
                CALL M3EXIT( PROGNAME, 0,0, MESG, 2 )
            END IF

            MESG = 'Treat applying default conversion factor (=1.0) as an error'
            CNVERROR = ENVYN ( 'DEFAULT_CONV_FAC_ERROR', MESG, .FALSE., IOS )
            IF ( IOS .GT. 0 ) THEN
                MESG = 'Bad environment variable "DEFAULT_CONV_FAC_ERROR"'
                CALL M3EXIT( PROGNAME, 0,0, MESG, 2 )
            END IF

C.............  Retrieve skipping speciation profile flags for zero emissions
            MESG = 'Setting for assigning profile for zero emissions'
            LNOZERO = ENVYN ( 'NO_SPC_ZERO_EMIS', MESG, .FALSE., IOS )
            IF ( IOS .GT. 0 ) THEN
                MESG = 'Bad environment variable "NO_SPC_ZERO_EMIS"'
                CALL M3EXIT( PROGNAME, 0,0, MESG, 2 )
            END IF

C.............  Figure out if SIC and/or MACT codes are available
            IF ( ASSOCIATED ( CISIC ) ) SICFLAG  = .TRUE.
            IF ( ASSOCIATED ( CMACT ) ) MACTFLAG = .TRUE.

            FIRSTIME = .FALSE.

        ENDIF

C.........  Initialize that combo file would need to be read
C           for this pollutant at first instance of COMBO in GSREF
        LRDCOMBO = .TRUE.

C.........  Set number of output fields for FMTCSRC to use
        SELECT CASE ( CATEGORY )
        CASE ( 'AREA' )
            NCOUT = 1
        CASE ( 'MOBILE' )
            NCOUT = NCHARS
        CASE ( 'POINT' )
            NCOUT = NCHARS
        END SELECT

C.........  Initialize matrices to 0.
        IF( MASSOUT ) THEN
            MASSMATX( :,1:MXSPEC ) = 0.    ! array
        END IF

        IF( MOLEOUT ) THEN
            MOLEMATX( :,1:MXSPEC ) = 0.    ! array
        END IF

C.........  Read emissions for current pollutant. If the emissions for a
C           source are zero, then the profile is not needed and the
C           massmatx and molematx calculations can be skipped.
C.........  Only do this when emissions (not VMT) in inventory (LVALCHK=TRUE)
        IF( LVALCHK ) CALL RDMAPPOL( NSRCIN, 1, 1, ENAM, EMISTMP )

C.........  Write pollutant of interest to the supplemental file
        WRITE( UDEV, '(A)' ) '"' // ENAM // '"'

C.........  Initialize index check
        NCHKCHR = NCHARS
        IF( JSCC .GT. 0 ) NCHKCHR = NCHARS - 1

C.........  Find index in complete list of pollutants and set length of name
        V  = INDEX1( ENAM, NIPPA, EANAM )

        DO S = 1, NSRCIN

C.............  If emissions are zero for this source, then skip it
          IF( LVALCHK .AND. LNOZERO .AND. EMISTMP( S ) <= 0 ) THEN

            CSRC  = CSOURC( S )
            WRITE( UDEV, '( 5X, 3A, 2X, A, I9 )' )
     &             QUOTE, TRIM( CSRC ), QUOTE, 'NFRAC=', 1
            WRITE( UDEV, '( 5X, 3A, F10.7 )' ) QUOTE, '      Drop', QUOTE, 1.0  ! assign "fake" speciation profile

          ELSE      !  if emissions are not zero for this source

            CSRC  = CSOURC( S )
            CFIP  = CSRC( 1:FIPLEN3 )
            CSTA  = CFIP( 1:STALEN3 )
            TSCC  = CSCC( S )
            SCCORIG = TSCC

C.............  Set type of SCC
            SCCFLAG = SETSCCTYPE ( TSCC )
            TSCCL = TSCC( 1:LSCCEND )

            CHK09 = CFIP // TSCC
            CHK08 = CFIP // TSCCL
            CHK06 = CSTA // TSCC
            CHK05 = CSTA // TSCCL

            IF( SICFLAG ) THEN
                CSIC = CISIC( S )
                CSICL = CSIC( 1:SICEXPLEN3+2 )
                CHK28 = CSTA // CSICL
                CHK29 = CSTA // CSIC
                CHK30 = CFIP // CSICL
                CHK31 = CFIP // CSIC
            END IF

            IF( MACTFLAG ) THEN
                CMCT  = CMACT( S )
                CHK33 = TSCC // CMCT
                CHK34 = CSTA // CMCT
                CHK35 = CSTA // TSCC // CMCT
                CHK36 = CFIP // CMCT
                CHK37 = CFIP // TSCC // CMCT
            END IF

            TSCCINIT = TSCC

C.............  Create selection
            SELECT CASE ( CATEGORY )

            CASE ( 'AREA' )   ! Already set above

            CASE ( 'MOBILE' )

            CASE ( 'POINT' )

                CHK16   = CSRC( 1:PTENDL3( 7 ) ) // TSCC
                CHK15   = CSRC( 1:PTENDL3( 6 ) ) // TSCC
                CHK14   = CSRC( 1:PTENDL3( 5 ) ) // TSCC
                CHK13   = CSRC( 1:PTENDL3( 4 ) ) // TSCC
                CHK12   = CSRC( 1:PTENDL3( 3 ) ) // TSCC
                CHK11   = CSRC( 1:PTENDL3( 2 ) ) // TSCC
                CHK10   = CSRC( 1:PTENDL3( 2 ) )

                CSRC5   = CSRC( 1:PTENDL3( 7 ) )
                CSRC4   = CSRC( 1:PTENDL3( 6 ) )
                CSRC3   = CSRC( 1:PTENDL3( 5 ) )
                CSRC2   = CSRC( 1:PTENDL3( 4 ) )
                CSRC1   = CSRC( 1:PTENDL3( 3 ) )

            CASE DEFAULT

            END SELECT

C.........................................................................
C.............  Now find and apply speciation profiles data
C.........................................................................

C.............  In the tables used in the following heirarchy, all cross-
C               reference entries are by definition, pollutant- specific.
C               The cross-reference tables (e.g,, CHRT02 come from MODXREF)

C.............  Try for pollutant-specific CHAR5 non-blank// SCC match; then
C                       pollutant-specific CHAR4 non-blank// SCC or blank match; then
C                       pollutant-specific CHAR3 non-blank// SCC or blank match; then
C                       pollutant-specific CHAR2 non-blank// SCC or blank match; then
C                       pollutant-specific CHAR1 non-blank// SCC or blank match; then
C                       pollutant-specific PLANT non-blank// SCC match; then
C                       pollutant-specific PLANT non-blank       match

            F6 = 0
            F5 = 0
            F4 = 0
            F3 = 0
            F2 = 0
            SELECT CASE( NCHKCHR )
            CASE( 7 )
                F6 = FINDC( CHK16, TXCNT( 16 ), CHRT16 )
            CASE( 6 )
                F5 = FINDC( CHK15, TXCNT( 15 ), CHRT15 )
            CASE( 5 )
                F4 = FINDC( CHK14, TXCNT( 14 ), CHRT14 )
            CASE( 4 )
                F3 = FINDC( CHK13, TXCNT( 13 ), CHRT13 )
            CASE( 3 )
                F2 = FINDC( CHK12, TXCNT( 12 ), CHRT12 )
            END SELECT

            IF( F6 .LE. 0 ) F6 = FINDC( CSRC5, TXCNT( 16 ), CHRT16 )
            IF( F5 .LE. 0 ) F5 = FINDC( CSRC4, TXCNT( 15 ), CHRT15 )
            IF( F4 .LE. 0 ) F4 = FINDC( CSRC3, TXCNT( 14 ), CHRT14 )
            IF( F3 .LE. 0 ) F3 = FINDC( CSRC2, TXCNT( 13 ), CHRT13 )
            IF( F2 .LE. 0 ) F2 = FINDC( CSRC1, TXCNT( 12 ), CHRT12 )
            F1 = FINDC( CHK11, TXCNT( 11 ), CHRT11 )
            F0 = FINDC( CHK10, TXCNT( 10 ), CHRT10 )

            IF( F6 .GT. 0 .AND. CSPT16(F6,V) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT16( F6,V )
                CALL SETSOURCE_SMATS
                CYCLE                       !  to end of sources-loop

            ELSEIF( F5 .GT. 0 .AND. CSPT15(F5,V) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT15( F5,V )
                CALL SETSOURCE_SMATS
                CYCLE                       !  to end of sources-loop

            ELSEIF( F4 .GT. 0 .AND. CSPT14(F4,V) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT14( F4,V )
                CALL SETSOURCE_SMATS
                CYCLE                       !  to end of sources-loop

            ELSEIF( F3 .GT. 0 .AND. CSPT13(F3,V) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT13( F3,V )
                CALL SETSOURCE_SMATS
                CYCLE                       !  to end of sources-loop

            ELSEIF( F2 .GT. 0 .AND. CSPT12(F2,V) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT12( F2,V )
                CALL SETSOURCE_SMATS
                CYCLE                       !  to end of sources-loop

            ELSEIF( F1 .GT. 0 .AND. CSPT11(F1,V) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT11( F1,V )
                CALL SETSOURCE_SMATS
                CYCLE                       !  to end of sources-loop

            ELSEIF( F0 .GT. 0 .AND. CSPT10(F0,V) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT10( F0,V )
                CALL SETSOURCE_SMATS
                CYCLE                       !  to end of sources-loop

            END IF

C.............  If MACT available in inventory...
            IF ( MACTFLAG ) THEN

C.............  Try for pollutant-specific FIPS code, SCC match, & MACT code; then
C                       pollutant-specific FIPS code & MACT code; then
C                       pollutant-specific Cy/st code, SCC match, & MACT code; then
C                       pollutant-specific Cy/st code & MACT code; then
C                       pollutant-specific SCC match & MACT code; then
C                       pollutant-specific MACT code
                F5 = FINDC( CHK37, TXCNT( 37 ), CHRT37 )
                F4 = FINDC( CHK36, TXCNT( 36 ), CHRT36 )
                F3 = FINDC( CHK35, TXCNT( 35 ), CHRT35 )
                F2 = FINDC( CHK34, TXCNT( 34 ), CHRT34 )
                F1 = FINDC( CHK33, TXCNT( 33 ), CHRT33 )
                F0 = FINDC( CMCT , TXCNT( 32 ), CHRT32 )

                IF( F5 .GT. 0 .AND. CSPT37( F5,V ) .NE. EMCMISS3 ) THEN
                    SPCODE = CSPT37( F5,V )
                    CALL SETSOURCE_SMATS
                    CYCLE                       !  to end of sources-loop

                ELSEIF(F4 .GT. 0 .AND. CSPT36(F4,V) .NE. EMCMISS3 ) THEN
                    SPCODE = CSPT36( F4,V )
                    CALL SETSOURCE_SMATS
                    CYCLE                       !  to end of sources-loop

                ELSEIF(F3 .GT. 0 .AND. CSPT35(F3,V) .NE. EMCMISS3 ) THEN
                    SPCODE = CSPT35( F3,V )
                    CALL SETSOURCE_SMATS
                    CYCLE                       !  to end of sources-loop

                ELSEIF(F2 .GT. 0 .AND. CSPT34(F2,V) .NE. EMCMISS3 ) THEN
                    SPCODE = CSPT34( F2,V )
                    CALL SETSOURCE_SMATS
                    CYCLE                       !  to end of sources-loop

                ELSEIF(F1 .GT. 0 .AND. CSPT33(F1,V) .NE. EMCMISS3 ) THEN
                    SPCODE = CSPT33( F1,V )
                    CALL SETSOURCE_SMATS
                    CYCLE                       !  to end of sources-loop

                ELSEIF(F0 .GT. 0 .AND. CSPT32(F0,V) .NE. EMCMISS3 ) THEN
                    SPCODE = CSPT32( F0,V )
                    CALL SETSOURCE_SMATS
                    CYCLE                       !  to end of sources-loop

                END IF
            END IF

C.............  If SIC available in inventory...
            IF ( SICFLAG ) THEN

C.............  Try for pollutant-specific FIPS code & SIC match; then
C                       pollutant-specific FIPS code & left SIC match; then
C                       pollutant-specific Cy/st code & SIC match; then
C                       pollutant-specific Cy/st code & left SIC match; then
C                       pollutant-specific SIC match; then
C                       pollutant-specific left SIC match

                F5 = FINDC( CHK31, TXCNT( 31 ), CHRT31 )
                F4 = FINDC( CHK30, TXCNT( 30 ), CHRT30 )
                F3 = FINDC( CHK29, TXCNT( 29 ), CHRT29 )
                F2 = FINDC( CHK28, TXCNT( 28 ), CHRT28 )
                F1 = FINDC( CSIC , TXCNT( 27 ), CHRT27 )
                F0 = FINDC( CSICL, TXCNT( 26 ), CHRT26 )

                IF( F5 .GT. 0 .AND. CSPT31( F5,V ) .NE. EMCMISS3 ) THEN
                    SPCODE = CSPT31( F5,V )
                    CALL SETSOURCE_SMATS
                    CYCLE                       !  to end of sources-loop

                ELSEIF(F4 .GT. 0 .AND. CSPT30(F4,V) .NE. EMCMISS3 ) THEN
                    SPCODE = CSPT30( F4,V )
                    CALL SETSOURCE_SMATS
                    CYCLE                       !  to end of sources-loop

                ELSEIF(F3 .GT. 0 .AND. CSPT29(F3,V) .NE. EMCMISS3 ) THEN
                    SPCODE = CSPT29( F3,V )
                    CALL SETSOURCE_SMATS
                    CYCLE                       !  to end of sources-loop

                ELSEIF(F2 .GT. 0 .AND. CSPT28(F2,V) .NE. EMCMISS3 ) THEN
                    SPCODE = CSPT28( F2,V )
                    CALL SETSOURCE_SMATS
                    CYCLE                       !  to end of sources-loop

                ELSEIF(F1 .GT. 0 .AND. CSPT27(F1,V) .NE. EMCMISS3 ) THEN
                    SPCODE = CSPT27( F1,V )
                    CALL SETSOURCE_SMATS
                    CYCLE                       !  to end of sources-loop

                ELSEIF(F0 .GT. 0 .AND. CSPT26(F0,V) .NE. EMCMISS3 ) THEN
                    SPCODE = CSPT26( F0,V )
                    CALL SETSOURCE_SMATS
                    CYCLE                       !  to end of sources-loop

                END IF
            END IF

C.............  Try for pollutant-specific FIPS code & SCC match; then
C                       pollutant-specific FIPS code & left SCC match; then
C                       pollutant-specific Cy/st code & SCC match; then
C                       pollutant-specific Cy/st code & left SCC match; then
C                       pollutant-specific SCC match; then
C                       pollutant-specific left SCC match

            F5 = FINDC( CHK09, TXCNT( 9 ), CHRT09 )
            F4 = FINDC( CHK08, TXCNT( 8 ), CHRT08 )
            F3 = FINDC( CHK06, TXCNT( 6 ), CHRT06 )
            F2 = FINDC( CHK05, TXCNT( 5 ), CHRT05 )
            F1 = FINDC( TSCC , TXCNT( 3 ), CHRT03 )
            F0 = FINDC( TSCCL, TXCNT( 2 ), CHRT02 )

C............. Check for mobile-specific matches that use a TSCC with
C              road class of zero and vehicle type. The assignment of
C              temporal profile based on  a vehicle type and no road class
C              comes after the road class only match (or TSCCL in CHRT08,
C              for example) but the match uses the full TSCC (or CHRT09, for
C              example).

            IF( F5 .GT. 0 .AND. CSPT09(F5,V) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT09( F5,V )
                CALL SETSOURCE_SMATS
                CYCLE                       !  to end of sources-loop

            ELSEIF( F4 .GT. 0 .AND. CSPT08(F4,V) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT08( F4,V )
                CALL SETSOURCE_SMATS
                CYCLE                       !  to end of sources-loop

            ELSEIF( F4B .GT. 0 .AND. CSPT09(F4B,V) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT09( F4B,V )
                CALL SETSOURCE_SMATS
                CYCLE                       !  to end of sources-loop

            ELSEIF( F3 .GT. 0 .AND. CSPT06(F3,V) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT06( F3,V )
                CALL SETSOURCE_SMATS
                CYCLE                       !  to end of sources-loop

            ELSEIF( F2 .GT. 0 .AND. CSPT05(F2,V) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT05( F2,V )
                CALL SETSOURCE_SMATS
                CYCLE                       !  to end of sources-loop

            ELSEIF( F2B .GT. 0 .AND. CSPT06(F2B,V) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT06( F2B,V )
                CALL SETSOURCE_SMATS
                CYCLE                       !  to end of sources-loop

            ELSEIF( F1 .GT. 0 .AND. CSPT03(F1,V) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT03( F1,V )
                CALL SETSOURCE_SMATS
                CYCLE                       !  to end of sources-loop

            ELSEIF( F0 .GT. 0 .AND. CSPT02(F0,V) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT02( F0,V )
                CALL SETSOURCE_SMATS
                CYCLE                       !  to end of sources-loop

            ELSEIF( F0B .GT. 0 .AND. CSPT03(F0B,V) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT03( F0B,V )
                CALL SETSOURCE_SMATS
                CYCLE                       !  to end of sources-loop

            END IF

C.............  Try for any FIPS code match
            F0 = FINDC( CFIP, TXCNT( 7 ), CHRT07 )

            IF( F0 .GT. 0 ) THEN
                SPCODE = CSPT07( F0,V )
                CALL SETSOURCE_SMATS
                CYCLE                       !  to end of sources-loop
            END IF

C.............  Try for any country/state code match (not, pol-specific)
            F0 = FINDC( CSTA, TXCNT( 4 ), CHRT04 )

            IF( F0 .GT. 0 ) THEN
                SPCODE = CSPT04( F0,V )
                CALL SETSOURCE_SMATS
                CYCLE                       !  to end of sources-loop
            END IF

C.............  For default speciation profile, make sure that it has been
C               defined for the current pollutant, and that we want to report
C               the use of defaults.
            IF( CSPT01( V ) .NE. EMCMISS3 .AND.
     &          REPDEFLT .AND. REPORT           ) THEN
                SPCODE = CSPT01( V )

                CALL FMTCSRC( CSRC, NCOUT, BUFFER, L2 )

                MESG = 'NOTE: Using default speciation profile "' //
     &                 SPCODE // '" for:'//
     &                 CRLF() // BLANK10 // BUFFER( 1:L2 ) //
     &                 CRLF() // BLANK10 //
     &                 'SCC: ' // TSCCINIT // ' POL: ' // EANAM( V )
                CALL M3MESG( MESG )

                CALL SETSOURCE_SMATS

            ELSEIF( CSPT01( V ) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT01( V )
                CALL SETSOURCE_SMATS

            ELSE

                EFLAG = .TRUE.
                IF( ERRCNT(1) < MXERR ) CALL REPORT_MISSING_DEFAULT
                ERRCNT(1) = ERRCNT(1) + 1

            END IF    !  if default profile code is available or not

          END IF

        END DO        !  end loop on source, S

        IF( EFLAG ) THEN
            MESG = 'Problem assigning speciation profiles to sources'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94300   FORMAT( A, I2.2, A, I2.2, A )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subroutine writes the message when a default
C               speciation profile is unavailable for a given pollutant
            SUBROUTINE REPORT_MISSING_DEFAULT

                CALL FMTCSRC( CSRC, NCOUT, BUFFER, L2 )

                MESG = 'ERROR: No speciation cross-reference ' //
     &                 'available (and no default) for:' //
     &                 CRLF() // BLANK10 // BUFFER( 1:L2 ) //
     &                 CRLF() // BLANK10 //
     &                 'SCC: ' // TSCCINIT // ' POL: ' // EANAM( V )

                CALL M3MESG( MESG )

            END SUBROUTINE REPORT_MISSING_DEFAULT

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C.............  This internal subprogram searches for the speciation profile
C               code in the abriged list (from MODSPRO) and if found, applies
C               the speciation factors for all species in that profile to
C               the speciation matrices.
C.............  Most variables are defined through host association.
            SUBROUTINE SETSOURCE_SMATS

C.............  Local variables
            INTEGER   F, I, J, K, N  ! indices and counters

            INTEGER   ITBL        ! position in full table of current profile
            INTEGER   NTBL        ! number of species of current profile
            INTEGER   NPCOMBO     ! tmp no. of profiles in combination profile
            REAL      CNVFAC      ! tmp pol-to-pol conversion factor

            LOGICAL   VALID, USECOMBO, USEFRACS

            INTEGER, SAVE :: CDEV = 0       ! unit number for GSPRO_COMBO
            LOGICAL, SAVE :: FIRSTCOMBO = .TRUE.  ! true: until first COMBO entry found in GSREF

            CHARACTER(24) :: PROGNAME = 'ASGNSPRO/SETSOURCE_SMATS'

C----------------------------------------------------------------------

C.............  If speciation profile code assigned is a combination profile...
C.............  This internal subroutine checks to see if the assignment has
C               been made to a COMBO profile and takes action accordingly.
C               Its main jobs are to read in the GSREF_COMBO file, find the
C               combination profile in the dataset, and to setup the number
C               of profiles in a combination profile.
C.............  NOTE:SREFDEX( SPCODE ) =
C                   0 -- not an SREF* profile
C                   -9999:  out-of-range/invalid SREF* profile
C                   1-SREFCNT:  SREF* subscript
            F = CMBDEX( SPCODE )

            IF ( F .LT. 0 ) THEN

                MESG = 'INTERNAL ERROR:  invalid CMB-fraction profile '
     &                 // SPCODE
                CALL M3EXIT( PROGNAME, 0,0, MESG, 2 )

            ELSE IF ( F .GT. 0 ) THEN

                NPCOMBO = CMBNP( F )
                DO I = 1, NPCOMBO
                    CCODE  ( I ) = ADJUSTR( CMBSPCD( F,I ) )
                    CWEIGHT( I ) = CMBWGHT( F,I )
                END DO
                USECOMBO = .FALSE.
                USEFRACS = .TRUE.
            ELSE IF ( ADJUSTL( SPCODE ) .EQ. 'COMBO' ) THEN
                    USECOMBO = .TRUE.
                    USEFRACS = .FALSE.

C................. For the first combo profile encountered, open the GSPRO_COMBO
C                  file. Arrays are passed back via the MODSPRO module
                IF ( FIRSTCOMBO ) THEN

                    CDEV = PROMPTFFILE(
     &                'Enter logical name for COMBINATION PROFILE file',
     &           .     TRUE., .TRUE., 'GSPRO_COMBO', PROGNAME )
                    IF ( CDEV .LT. 0 ) THEN
                        MESG = 'ERROR:  could not open "GSPRO_COMBO"'
                        CALL M3EXIT( PROGNAME, 0,0, MESG, 2 )
                    END IF
                    FIRSTCOMBO = .FALSE.

                END IF

C.................  If first time for this pollutant, read the combo file
C.................  It will store only for pollutant/type ENAM and SPCMAT_PERIOD
C                        from the environment
                IF ( LRDCOMBO ) THEN
                    CALL RDCOMBO( CDEV, ENAM )
                    LRDCOMBO = .FALSE.

                END IF

C.................  Match the FIPS code to the combination profile
                F = FINDC( CFIP, NINVIFIP, INVCFIP )

C.................  If FIPS code is not found, then internal error since these
C                   arrays are supposed to be built from inventory list of FIPs codes
                IF ( F <= 0 ) THEN
                    MESG = 'INTERNAL ERROR: FIPS code not found '//
     &                     'in INVCFIP array from MODLISTS for:'//
     &                     CRLF() // BLANK10 // 'FIP: '// TRIM( CFIP )
                    IF( ERRCNT(2) < MXERR ) CALL M3MESG( MESG )
                    ERRCNT(2) = ERRCNT(2) + 1
                    EFLAG = .TRUE.

C.................  Otherwise, store information from combination data, if found
                ELSE

                    NPCOMBO = CMBNP( F )

C.....................  If FIPs code found, set the arrays for SETSOURCE_SMATS
                    IF ( NPCOMBO > 0 ) THEN

                        DO I = 1, NPCOMBO
                            CCODE  ( I ) = ADJUSTR( CMBSPCD( F,I ) )
                            CWEIGHT( I ) = CMBWGHT( F,I )
                        END DO

C.....................  Otherwise, a COMBO was specified, but none found
C.....................  Note that NPCOMBO = 0 should cause loop in SETSOURCE_SMAT
C                          to be skipped
                    ELSE
                        MESG = 'ERROR: No FIPS-specific entry '//
     &                         'in GSPRO_COMBO for:'//
     &                         CRLF() // BLANK10 //
     &                         'FIP: '// TRIM(CFIP) // ' and'//
     &                         'pollutant "'// TRIM( ENAM ) //'"'
                        IF( ERRCNT(3) < MXERR ) CALL M3MESG( MESG )
                        ERRCNT(3) = ERRCNT(3) + 1
                        EFLAG = .TRUE.

                    END IF

                END IF

C.............  If this entry is not a combination profile, setup to simply
C               use the single profile specified
            ELSE

                NPCOMBO      = 1
                CCODE  ( 1 ) = SPCODE
                CWEIGHT( 1 ) = 1.0
                USECOMBO     = .FALSE.
                USEFRACS     = .FALSE.

            END IF

C.............  Write speciation profile codes and fractions for this source
C               to the speciation supplemental file (to be used by SMKREPORT).

            WRITE( UDEV, '( 5X, 3A, 2X, A, I9 )' )
     &             QUOTE, TRIM( CSRC ), QUOTE, 'NFRAC=', NPCOMBO
            WRITE( UDEV, '( 5X, 999( 3A, F10.7, : 2X ) )' )
     &             ( QUOTE, TRIM( CCODE(K) ), QUOTE, CWEIGHT(K), K = 1, NPCOMBO )

C.............  Loop over all profiles in combination profile
C.............  Note that for non-combo profiles, NPCOMBO has been set to 1
C                  and CCODE(1) to SPCODE.
            DO NP = 1, NPCOMBO

                VALID = .TRUE.
                K = MAX( INDEX1( CCODE( NP ), NSPROF, SPROFN ),0 )

C.................  If specific profile is not found in set of profiles,
C                   try to apply the default for this pollutant.
C                   If COMBO or fractional, give a different error and
C                   don't look for the default)
                IF( K .LE. 0 ) THEN

                    CALL FMTCSRC( CSRC, NCOUT, BUFFER, L2 )

                    IF ( USECOMBO ) THEN

                        MESG = 'ERROR: Speciation profile "'// CCODE(NP)//
     &                     '" is not in profiles, but it was assigned' //
     &                     CRLF() // BLANK10 // 'as COMBO to source:' //
     &                     CRLF() // BLANK10 // BUFFER( 1:L2 ) //
     &                     ' SCC: ' // TSCCINIT // ' POL: ' // EANAM( V )
                        IF( ERRCNT(5) <= MXERR ) CALL M3MESG( MESG )
                        ERRCNT(5) = ERRCNT(5) + 1
                        EFLAG = .TRUE.
                        CYCLE

                    ELSE IF ( USEFRACS ) THEN

                        MESG = 'ERROR: Speciation profile "'// CCODE(NP)//
     &                     '" is not in profiles, but it was assigned' //
     &                     CRLF() // BLANK10 //
     &                     'as as fractional-profile to source:' //
     &                     CRLF() // BLANK10 // BUFFER( 1:L2 ) //
     &                     ' SCC: ' // TSCCINIT // ' POL: ' // EANAM( V )
                        IF( ERRCNT(5) <= MXERR ) CALL M3MESG( MESG )
                        ERRCNT(5) = ERRCNT(5) + 1
                        EFLAG = .TRUE.
                        CYCLE

                    ELSE

                        MESG = 'WARNING: Speciation profile "'// CCODE(NP)//
     &                     '" is not in profiles, but it was assigned' //
     &                     CRLF() // BLANK10 // 'to source:' //
     &                     CRLF() // BLANK10 // BUFFER( 1:L2 ) //
     &                     ' SCC: ' // TSCCINIT // ' POL: ' // EANAM( V )
                        IF( WARNCNT(1) <= MXWARN ) CALL M3MESG( MESG )
                        WARNCNT(1) = WARNCNT(1) + 1

                        K = MAX( FINDC( CSPT01( V ), NSPROF, SPROFN ), 0 )

                        IF( CSPT01( V ) .NE. EMCMISS3 .AND. K .GT. 0 ) THEN
                            MESG = BLANK5 // 'Using default profile ' //
     &                         CSPT01( V )
                            CALL M3MESG( MESG )

                        ELSE
                            EFLAG = .TRUE.
                            VALID = .FALSE.
                            IF( ERRCNT(4) < MXERR ) CALL REPORT_MISSING_DEFAULT
                            ERRCNT(4) = ERRCNT(4) + 1
                        END IF

                    END IF

                END IF

C.................  Check if pollutant-to-pollutant conversion factor is available
C                   by speciation profile, by checking if their count > 0
                IF( NCNV4 .GT. 0 .AND. VALID ) THEN

                    F1 = FINDC( ADJUSTL( CCODE(NP) ), NCNV4, CNVRT04 )

                    IF( F1 .GT. 0 )THEN
                        CNVFAC = CNVFC04( F1,V )

C.....................  CNVFC00( V ) will equal 1.0 if it has not been set, so
C                       there is no need for error checking
                    ELSE
                        CNVFAC = CNVFC00( V )

                        CALL FMTCSRC( CSRC, NCOUT, BUFFER, L2 )
                        IF( CNVFAC == 1.0 .AND. CNVFLAG( V ) ) THEN
                            MESG='ERROR: Missing assigned speciation profile "'//
     &                          CCODE(NP)//'" from GSCNV input file'//
     &                          CRLF() // BLANK10 // BUFFER( 1:L2 ) //
     &                          ' SCC: ' // TSCCINIT // ' POL: ' // EANAM( V )
                            IF( ERRCNT(6) <= MXERR .AND. CNVERROR ) THEN
                                CALL M3MESG( MESG )
                                ERRCNT(6) = ERRCNT(6) + 1
                                EFLAG = .TRUE.
                            END IF
                        END IF

                    END IF

C.................  Screen for pollutant-to-pollutant conversion factors by checking
C                   if they have been allocated
                ELSE IF( ALLOCATED( CNVRT03 ) .AND. VALID ) THEN

C.....................  If have a combination profile, then this is not a valid
C                       way to set the pollutant-to-pollutant factors
                    IF( USECOMBO ) THEN
                        MESG = 'ERROR: Can not assign pollutant-to-'//
     &                         'pollutant factors using FIPS/SCC '//
     &                         'using combination profiles.'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    ELSE IF ( USEFRACS ) THEN
                        MESG = 'ERROR: Can not assign pollutant-to-'//
     &                   'pollutant factors using  fractional profiles.'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF

C.....................  Try for pollutant-specific FIPS code & SCC match; then
C                           pollutant-specific Cy/st code & SCC match; then
C                           pollutant-specific SCC match
C                           pollutant-specific roadway type match
C                           pollutant-specific vehicle type match

                    F5 = FINDC( CFIP // SCCORIG , NCNV3, CNVRT03 )
                    F4 = FINDC( CSTA // SCCORIG , NCNV2, CNVRT02 )
                    F3 = FINDC( TSCC  , NCNV1, CNVRT01 )
                    F2 = FINDC( CHKRWT, NCNV1, CNVRT01 )
                    F1 = FINDC( SCCORIG( 1:LSCCEND ), NCNV1, CNVRT01 )

                    IF( F5 .GT. 0  ) THEN
                        CNVFAC = CNVFC03( F5,V )

                    ELSE IF( F4 .GT. 0 )THEN
                        CNVFAC = CNVFC02( F4,V )

                    ELSE IF( F3 .GT. 0 )THEN
                        CNVFAC = CNVFC01( F3,V )

                    ELSE IF( F2 .GT. 0 )THEN
                        CNVFAC = CNVFC01( F2,V )

                    ELSE IF( F1 .GT. 0 )THEN
                        CNVFAC = CNVFC01( F1,V )

C.....................  Error check there is no conversion fac found
                    ELSE
                        CNVFAC = CNVFC00( V )
                        CALL FMTCSRC( CSRC, NCOUT, BUFFER, L2 )
                        IF( CNVFAC == 1.0 .AND. CNVFLAG( V ) ) THEN
                            MESG='ERROR: Missing pollutant conversion'//
     &                          ' factor for source:' //
     &                          CRLF() // BLANK10 // BUFFER( 1:L2 ) //
     &                          ' SCC: ' // TSCCINIT // ' POL: ' // EANAM( V )
                            IF( ERRCNT(7) <= MXERR .AND. CNVERROR ) THEN
                                CALL M3MESG( MESG )
                                ERRCNT(7) = ERRCNT(7) + 1
                                EFLAG = .TRUE.
                            END IF
                        END IF

                    END IF

C.................  If they don't exist, simply set the conversion factor to one
                ELSE

                    CNVFAC = 1.

                END IF

C.................  Now that the default profile has been tried, check one last time
C                   for K and then apply speciation factors
                IF( K .GT. 0 ) THEN

C.....................  Get indices to full speciation table
                    ITBL = IDXSPRO ( K )
                    NTBL = NSPECIES( K )

                    I = ITBL - 1
                    DO N = 1, NTBL

                        I = I + 1
                        J = IDXSSPEC( K,N )

                        IF( MASSOUT ) THEN
                            MASSMATX( S,J )= MASSMATX( S,J ) +
     &                            CWEIGHT( NP ) * CNVFAC * MASSFACT( I )
                        END IF

                        IF( MOLEOUT ) THEN
                            MOLEMATX( S,J )= MOLEMATX( S,J ) +
     &                            CWEIGHT( NP ) * CNVFAC * MOLEFACT( I )
                        END IF

                    END DO

                END IF

            END DO  ! Loop over combination profiles, index NP

            END SUBROUTINE SETSOURCE_SMATS

        END SUBROUTINE ASGNSPRO
