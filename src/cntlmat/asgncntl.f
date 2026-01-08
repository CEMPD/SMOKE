
        SUBROUTINE ASGNCNTL( NSRCIN, WDEV, PKTTYP, PNAM, DATSPFLAG, 
     &                       SINDX )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      For each source and pollutant, find the most specific entry in a
C      control data table to apply.  The control data of interest can vary,
C      but no matter the control packet being processed, the control data
C      table index is stored in the ICTL* grouped x-ref tables from MODXREF.
C      The source-based stored index for each control packet comes in
C      as a subroutine argument (SINDX), so this can change for each call to
C      this routine.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 3/99 by M. Houyoux
C
C************************************************************************
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
C*************************************************************************

C...........   MODULES for public variables   
C...........   This module contains the source ararys
        USE M3UTILIO

        USE MODSOURC, ONLY: CSOURC, CSCC, CMACT, CISIC, CSRCTYP

C.........  This module is for cross reference tables
        USE MODXREF, ONLY: ADDPS, TXCNT,
     &          CHRT02A, CHRT02B, CHRT02C, CHRT03, CHRT04,
     &          CHRT05A, CHRT05B, CHRT05C, CHRT06, CHRT07,
     &          CHRT08A, CHRT08B, CHRT08C, CHRT09, CHRT10,
     &          CHRT11, CHRT12, CHRT13, CHRT14, CHRT15, CHRT16,
     &          CHRT26, CHRT27, CHRT28, CHRT29, CHRT30, CHRT31,
     &          CHRT32, CHRT33, CHRT34, CHRT35, CHRT36, CHRT37, CHRT38,
     &          ICTL01, ICTL02A, ICTL02B, ICTL02C, ICTL03, ICTL04,
     &          ICTL05A, ICTL05B, ICTL05C, ICTL06, ICTL07,
     &          ICTL08A, ICTL08B, ICTL08C, ICTL09, ICTL10,
     &          ICTL11, ICTL12, ICTL13, ICTL14, ICTL15, ICTL16,
     &          ICTL26, ICTL27, ICTL28, ICTL29, ICTL30, ICTL31,
     &          ICTL32, ICTL33, ICTL34, ICTL35, ICTL36, ICTL37, ICTL38

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, NCHARS, JSCC, NIPPA, EANAM,
     &                     SCCLEV1, SCCLEV2, SCCLEV3

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  i/o api constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
C       CHARACTER(2)    CRLF
C       LOGICAL         ENVYN
C       INTEGER         FIND1
C       INTEGER         FINDC
C       INTEGER         INDEX1
        LOGICAL         SETSCCTYPE

C        EXTERNAL CRLF, ENVYN, FIND1, FINDC, INDEX1, SETSCCTYPE
        EXTERNAL     SETSCCTYPE

C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NSRCIN         ! number of sources
        INTEGER     , INTENT (IN) :: WDEV           ! warning file
        CHARACTER(*), INTENT (IN) :: PKTTYP         ! packet type of interest
        CHARACTER(*), INTENT (IN) :: PNAM   ! pollutant names
        LOGICAL     , INTENT(OUT) :: DATSPFLAG      ! true: >= 1 data-spec asgn
        INTEGER     , INTENT(OUT) :: SINDX( NSRCIN )! idx to control table

C.........  Other local variables
        INTEGER          I, J, K, L, L2, S, V    !  counters and indices

        INTEGER          F0, F1, F1B, F2, F3, F4, F5, F6  ! tmp find indices
        INTEGER          F7, F8, F9, F10, F11        ! tmp find indices
        INTEGER          IDX    !  tmp index to control data table
        INTEGER          ISTAT  !  tmp indicator for internal subprogram calls
        INTEGER          NCHKCHR          ! position of last non-SCC src char

        LOGICAL       :: EFLAG    = .FALSE. ! true: error detected
        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time subrtn called
        LOGICAL, SAVE :: MACTFLAG = .FALSE. ! true: MACT and src type available in inventory
        LOGICAL, SAVE :: REPDEFLT = .TRUE.  ! true: report when defaults used
        LOGICAL, SAVE :: SICFLAG  = .FALSE. ! true: SIC available in inventory
        LOGICAL, SAVE :: SICFIRST = .TRUE.  ! true: assign by SIC before SCC
        LOGICAL       :: SICXREF  = .FALSE. ! true: SIC assignments in x-ref
        LOGICAL          SCCFLAG            ! true: SCC type is different from previous

        CHARACTER(8), SAVE :: FMTFIP   ! format for writing FIPS code
        CHARACTER(256)     BUFFER   ! source fields buffer
        CHARACTER(256)     MESG     ! message buffer
        CHARACTER(STALEN3) CSTA     ! tmp Country/state code
        CHARACTER(STSLEN3) CSTAS_A  ! tmp Country/state code // level 1 SCC
        CHARACTER(STSLEN3) CSTAS_B  ! tmp Country/state code // level 2 SCC
        CHARACTER(STSLEN3) CSTAS_C  ! tmp Country/state code // level 3 SCC
        CHARACTER(STSLEN3) CSTAS_D  ! tmp Country/state code // all SCC
        CHARACTER(SCCLEN3) TSCC_A   ! tmp level 1 of SCC
        CHARACTER(SCCLEN3) TSCC_B   ! tmp level 2 of SCC
        CHARACTER(SCCLEN3) TSCC_C   ! tmp level 3 of SCC
        CHARACTER(SCCLEN3) TSCC_D   ! tmp level 4 (all) of SCC
        CHARACTER(SCCLEN3) TSCC_BUF ! tmp 10-digit SCC
        CHARACTER(SRCLEN3) CSRC     ! tmp source chars string
        CHARACTER(SS0LEN3) :: CSSC0 = ' ' ! tmp FIPS // Plant // SCC 
        CHARACTER(SS1LEN3) :: CSSC1 = ' ' ! tmp source chars through char1  // SCC 
        CHARACTER(SS2LEN3) :: CSSC2 = ' ' ! tmp source chars through char2  // SCC 
        CHARACTER(SS3LEN3) :: CSSC3 = ' ' ! tmp source chars through char3  // SCC 
        CHARACTER(SS4LEN3) :: CSSC4 = ' ' ! tmp source chars through char4  // SCC 
        CHARACTER(SS5LEN3) :: CSSC5 = ' ' ! tmp source chars through char5  // SCC 
        CHARACTER(SS0LEN3) :: CSRC0 = ' ' ! tmp FIPS // Plant
        CHARACTER(SS1LEN3) :: CSRC1 = ' ' ! tmp source chars through char1
        CHARACTER(SS2LEN3) :: CSRC2 = ' ' ! tmp source chars through char2
        CHARACTER(SS3LEN3) :: CSRC3 = ' ' ! tmp source chars through char3
        CHARACTER(SS4LEN3) :: CSRC4 = ' ' ! tmp source chars through char4
        CHARACTER(SS5LEN3) :: CSRC5 = ' ' ! tmp source chars through char5
        CHARACTER(FIPLEN3) CFIP     ! tmp (character) FIPS code
        CHARACTER(FPSLEN3) CFIPS_A  ! tmp FIPS code // level 1 of SCC
        CHARACTER(FPSLEN3) CFIPS_B  ! tmp FIPS code // level 2 of SCC
        CHARACTER(FPSLEN3) CFIPS_C  ! tmp FIPS code // level 3 of SCC
        CHARACTER(FPSLEN3) CFIPS_D  ! tmp FIPS code // level 4 (all) of SCC
        CHARACTER(SICLEN3) :: CSIC = ' '     ! tmp char SIC
        CHARACTER(SICLEN3) :: CSICL = ' '    ! tmp char left SIC
        CHARACTER(STILEN3) :: CSTASIC = ' '  ! tmp Country/state // char SIC
        CHARACTER(STILEN3) :: CSTASICL = ' ' ! tmp Country/state // char left SIC
        CHARACTER(FPILEN3) :: CFIPSIC = ' '  ! tmp FIPS code // char SIC
        CHARACTER(FPILEN3) :: CFIPSICL = ' ' ! tmp FIPS code // char left SIC
        CHARACTER(MACLEN3) :: CMCT = ' '     ! tmp MACT code
        CHARACTER(MSCLEN3) :: CMSCC = ' '    ! tmp SCC // MACT
        CHARACTER(MSTLEN3) :: CMST = ' '     ! tmp Country/state code // MACT
        CHARACTER(MSSLEN3) :: CMSTSC = ' '   ! tmp Country/state code // SCC // MACT
        CHARACTER(MFPLEN3) :: CMFP = ' '     ! tmp FIPS code // MACT
        CHARACTER(MFSLEN3) :: CMFPSC = ' '   ! tmp FIPS code // SCC // MACT
        CHARACTER(STPLEN3) :: CSTYP = ' '    ! tmp source type code
        CHARACTER(FPMLEN3) :: CFPM = ' '     ! tmp FIPS code // Plant // MACT

        CHARACTER(16) :: PROGNAME = 'ASGNCNTL' ! program name

C***********************************************************************
C   begin body of subroutine ASGNCNTL

C.........  For first time routine is called ...
        IF( FIRSTIME ) THEN

C.............  Retrieve environment variables
            MESG = 'Switch for reporting default control factors'
            REPDEFLT = ENVYN ( 'REPORT_DEFAULTS', MESG, .TRUE., I )

            MESG = 'Match sources by SIC before SCC'
            SICFIRST = ENVYN ( 'XREF_SICOVERSCC', MESG, .TRUE., I )

C.............  Set up format for creating character FIPS code for non-point
            IF( CATEGORY .NE. 'POINT' ) THEN
                WRITE( FMTFIP, 94300 ) '(I', FIPLEN3, '.', FIPLEN3, ')'
            ENDIF

C.............  Figure out if SIC and/or MACT and src type codes are available
            IF ( ASSOCIATED ( CISIC ) ) SICFLAG  = .TRUE.
            IF ( ASSOCIATED ( CMACT ) ) MACTFLAG = .TRUE.

            FIRSTIME = .FALSE.

        ENDIF

C.........  Set flag to indicate if SIC-specific assignments are
C           available in the cross-reference file
        SICXREF = ( MAXVAL( TXCNT( 26:31 ) ) .GT. 0 )

C.........  Initialize SINDX
        SINDX = 0      ! Array

C.........  For each pollutant of interest
        DATSPFLAG = .FALSE.

C.........  Initialize index check
        NCHKCHR = NCHARS
        IF( JSCC .GT. 0 ) NCHKCHR = NCHARS - 1

C.............  Write message, with different message for projections because
C               they are not yet data-specific
            IF( PKTTYP .NE. 'PROJECTION' ) THEN
                MESG = BLANK5 // 'Assigning controls for pollutant "' //
     &                 TRIM( PNAM ) // '"...'
                V = INDEX1( PNAM, NIPPA, EANAM )

            ELSE IF ( PNAM .EQ. 'all' ) THEN
                MESG = BLANK5 // 'Assigning projection factors ' //
     &                 'to all pollutants...'
                V = 1
            ELSE
                MESG = BLANK5 // 'Assigning projection factors ' //
     &                 'for pollutant "' //TRIM( PNAM ) // '"...'
                V = INDEX1( PNAM, NIPPA, EANAM )
            END IF
            CALL M3MSG2( MESG )

            DO S = 1, NSRCIN
   
                CSRC    = CSOURC( S )
                CFIP    = CSRC( 1:FIPLEN3 )
                CSTA    = CFIP( 1:STALEN3 )
                                
C.................  Set type of SCC                
                SCCFLAG = SETSCCTYPE ( CSCC( S ) )
                TSCC_D  = CSCC( S )             ! Level 4 (all)
                TSCC_A  = TSCC_D( 1:SCCLEV1 )   ! Level 1
                TSCC_B  = TSCC_D( 1:SCCLEV2 )   ! Level 2
                TSCC_C  = TSCC_D( 1:SCCLEV3 )   ! Level 3

                CFIPS_A = CFIP // TSCC_A
                CFIPS_B = CFIP // TSCC_B
                CFIPS_C = CFIP // TSCC_C
                CFIPS_D = CFIP // TSCC_D

                CSTAS_A = CSTA // TSCC_A
                CSTAS_B = CSTA // TSCC_B
                CSTAS_C = CSTA // TSCC_C
                CSTAS_D = CSTA // TSCC_D

                IF( SICFLAG ) THEN
                    CSIC = CISIC( S )
                    CSICL    = CSIC( 1:SICEXPLEN3+2 )
                    CSTASICL = CSTA // CSICL
                    CSTASIC  = CSTA // CSIC
                    CFIPSICL = CFIP // CSICL
                    CFIPSIC  = CFIP // CSIC
                END IF

                IF( MACTFLAG ) THEN
                    CMCT   = CMACT( S )
                    CMFPSC = CFIP // TSCC_D // CMCT
                    CMFP   = CFIP // CMCT
                    CMSTSC = CSTA // TSCC_D // CMCT
                    CMST   = CSTA // CMCT
                    CMSCC  = TSCC_D // CMCT
                    CSTYP  = CSRCTYP( S )
                END IF

C.................  Create selection 
                SELECT CASE ( CATEGORY )

                CASE ( 'AREA' ) 

                CASE ( 'MOBILE' )

                CASE ( 'POINT' )

                    CSRC5   = CSRC( 1:PTENDL3( 7 ) ) 
                    CSRC4   = CSRC( 1:PTENDL3( 6 ) ) 
                    CSRC3   = CSRC( 1:PTENDL3( 5 ) ) 
                    CSRC2   = CSRC( 1:PTENDL3( 4 ) ) 
                    CSRC1   = CSRC( 1:PTENDL3( 3 ) ) 
                    CSRC0   = CSRC( 1:PTENDL3( 2 ) )

C...................  Also set search field for SCC-specific assignment with all
C                     source characteristics.
                    CSSC5   = CSRC( 1:PTENDL3( 7 ) ) // TSCC_D
                    CSSC4   = CSRC( 1:PTENDL3( 6 ) ) // TSCC_D
                    CSSC3   = CSRC( 1:PTENDL3( 5 ) ) // TSCC_D
                    CSSC2   = CSRC( 1:PTENDL3( 4 ) ) // TSCC_D
                    CSSC1   = CSRC( 1:PTENDL3( 3 ) ) // TSCC_D
                    CSSC0   = CSRC( 1:PTENDL3( 2 ) ) // TSCC_D

                    IF ( MACTFLAG ) CFPM = CSRC( 1:PTENDL3(2) ) // CMCT
                    
                CASE DEFAULT

                END SELECT


C.................  Try for pollutant-specific CHAR5 non-blank// SCC match; then
C                           pollutant-specific CHAR4 non-blank// SCC or blank match; then
C                           pollutant-specific CHAR3 non-blank// SCC or blank match; then
C                           pollutant-specific CHAR2 non-blank// SCC or blank match; then
C                           pollutant-specific CHAR1 non-blank// SCC or blank match; then
C                           pollutant-specific PLANT non-blank// SCC match; then
C                           pollutant-specific PLANT non-blank// MACT match; then
C                           pollutant-specific PLANT non-blank       match
                F6 = 0
                F5 = 0
                F4 = 0
                F3 = 0
                F2 = 0
                F1 = 0
                F1B= 0
                F0 = 0
                SELECT CASE( NCHKCHR )
                CASE( 7 )
                    F6 = FINDC( CSSC5, TXCNT( 16 ), CHRT16 )
                CASE( 6 )
                    F5 = FINDC( CSSC4, TXCNT( 15 ), CHRT15 )
                CASE( 5 )
                    F4 = FINDC( CSSC3, TXCNT( 14 ), CHRT14 )
                CASE( 4 )
                    F3 = FINDC( CSSC2, TXCNT( 13 ), CHRT13 )
                CASE( 3 )
                    F2 = FINDC( CSSC1, TXCNT( 12 ), CHRT12 )
                END SELECT

                IF( F6 .LE. 0 ) F6 = FINDC( CSRC5, TXCNT( 16 ), CHRT16 )
                IF( F5 .LE. 0 ) F5 = FINDC( CSRC4, TXCNT( 15 ), CHRT15 ) 
                IF( F4 .LE. 0 ) F4 = FINDC( CSRC3, TXCNT( 14 ), CHRT14 ) 
                IF( F3 .LE. 0 ) F3 = FINDC( CSRC2, TXCNT( 13 ), CHRT13 ) 
                IF( F2 .LE. 0 ) F2 = FINDC( CSRC1, TXCNT( 12 ), CHRT12 ) 
                F1 = FINDC( CSSC0, TXCNT( 11 ), CHRT11 )
                F1B= FINDC( CFPM,  TXCNT( 38 ), CHRT38 )  ! FIP, Plant, MACT
                F0 = FINDC( CSRC0, TXCNT( 10 ), CHRT10 )

C..................  Do an additional check for the most detailed source assignment
C                    with SCC as well.

                IF( F6 .GT. 0 .AND. ICTL16(F6,V) .GE. ADDPS ) THEN
                    IDX = ICTL16( F6,V ) - ADDPS
                    CALL SETSOURCE_CONTROL_INDEX
                    DATSPFLAG = .TRUE.
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F5 .GT. 0 .AND. ICTL15(F5,V) .GE. ADDPS ) THEN
                    IDX = ICTL15( F5,V ) - ADDPS
                    CALL SETSOURCE_CONTROL_INDEX
                    DATSPFLAG = .TRUE.
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F4 .GT. 0 .AND. ICTL14(F4,V) .GE. ADDPS ) THEN
                    IDX = ICTL14( F4,V ) - ADDPS
                    CALL SETSOURCE_CONTROL_INDEX
                    DATSPFLAG = .TRUE.
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F3 .GT. 0 .AND. ICTL13(F3,V) .GE. ADDPS ) THEN
                    IDX = ICTL13( F3,V ) - ADDPS
                    CALL SETSOURCE_CONTROL_INDEX
                    DATSPFLAG = .TRUE.
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F2 .GT. 0 .AND. ICTL12(F2,V) .GE. ADDPS ) THEN
                    IDX = ICTL12( F2,V ) - ADDPS
                    CALL SETSOURCE_CONTROL_INDEX
                    DATSPFLAG = .TRUE.
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F1 .GT. 0 .AND. ICTL11(F1,V) .GE. ADDPS ) THEN
                    IDX = ICTL11( F1,V ) - ADDPS
                    CALL SETSOURCE_CONTROL_INDEX
                    DATSPFLAG = .TRUE.
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F1B .GT. 0 .AND. ICTL38(F1B,V) .GE. ADDPS ) THEN
                    IDX = ICTL38( F1B,V ) - ADDPS
                    CALL SETSOURCE_CONTROL_INDEX
                    DATSPFLAG = .TRUE.
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F0 .GT. 0 .AND. ICTL10(F0,V) .GE. ADDPS ) THEN
                    IDX = ICTL10( F0,V ) - ADDPS
                    CALL SETSOURCE_CONTROL_INDEX
                    DATSPFLAG = .TRUE.
                    CYCLE                       !  to end of sources-loop

                END IF

C.................  Try for any CHAR5 non-blank // SCC match; then
C                           any CHAR4 non-blank // SCC match; then
C                           any CHAR3 non-blank // SCC match; then
C                           any CHAR2 non-blank // SCC match; then
C                           any CHAR1 non-blank // SCC match; then
C                           any PLANT non-blank // SCC match; then
C                           any PLANT non-blank        match

                IF( F6 .GT. 0 .AND. ICTL16(F6,V) .NE. IMISS3 ) THEN
                    IDX = ICTL16( F6,V )
                    CALL SETSOURCE_CONTROL_INDEX
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F5 .GT. 0 .AND. ICTL15(F5,V) .NE. IMISS3 ) THEN
                    IDX = ICTL15( F5,V )
                    CALL SETSOURCE_CONTROL_INDEX
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F4 .GT. 0 .AND. ICTL14(F4,V) .NE. IMISS3 ) THEN
                    IDX = ICTL14( F4,V )
                    CALL SETSOURCE_CONTROL_INDEX
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F3 .GT. 0 .AND. ICTL13(F3,V) .NE. IMISS3 ) THEN
                    IDX = ICTL13( F3,V )
                    CALL SETSOURCE_CONTROL_INDEX
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F2 .GT. 0 .AND. ICTL12(F2,V) .NE. IMISS3 ) THEN
                    IDX = ICTL12( F2,V )
                    CALL SETSOURCE_CONTROL_INDEX
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F1 .GT. 0 .AND. ICTL11(F1,V) .NE. IMISS3 ) THEN
                    IDX = ICTL11( F1,V )
                    CALL SETSOURCE_CONTROL_INDEX
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F1B .GT. 0 .AND. ICTL38(F1B,V) .NE. IMISS3 ) THEN
                    IDX = ICTL38( F1B,V )
                    CALL SETSOURCE_CONTROL_INDEX
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F0 .GT. 0 .AND. ICTL10(F0,V) .NE. IMISS3 ) THEN
                    IDX = ICTL10( F0,V )
                    CALL SETSOURCE_CONTROL_INDEX
                    CYCLE                       !  to end of sources-loop

                END IF

C.................  If MACT available in inventory...
                IF ( MACTFLAG ) THEN
                
C.................  Try for pollutant-specific FIPS code, SCC, & MACT match; then
C                           pollutant-specific FIPS code & MACT match; then
C                           pollutant-specific Cy/st code, SCC, & MACT match; then
C                           pollutant-specific Cy/st code & MACT match; then
C                           pollutant-specific SCC & MACT match; then
C                           pollutant-specific MACT match
                    F5 = FINDC( CMFPSC, TXCNT( 37 ), CHRT37 )
                    F4 = FINDC( CMFP  , TXCNT( 36 ), CHRT36 )
                    F3 = FINDC( CMSTSC, TXCNT( 35 ), CHRT35 )
                    F2 = FINDC( CMST  , TXCNT( 34 ), CHRT34 )
                    F1 = FINDC( CMSCC , TXCNT( 33 ), CHRT33 )
                    F0 = FINDC( CMCT  , TXCNT( 32 ), CHRT32 )
                    
                    IF( F5 .GT. 0 .AND. ICTL37( F5,V ) .GE. ADDPS ) THEN
                        IDX = ICTL37( F5,V ) - ADDPS
                        IF( PKTTYP == 'MACT' ) CALL CHECK_SRC_TYPE
                        CALL SETSOURCE_CONTROL_INDEX
                        IF( IDX /= 0 ) DATSPFLAG = .TRUE.
                        CYCLE                       !  to end of sources-loop
                        
                    ELSEIF(F4 .GT. 0 .AND. ICTL36(F4,V) .GE. ADDPS) THEN
                        IDX = ICTL36( F4,V ) - ADDPS
                        IF( PKTTYP == 'MACT' ) CALL CHECK_SRC_TYPE
                        CALL SETSOURCE_CONTROL_INDEX
                        IF( IDX /= 0 ) DATSPFLAG = .TRUE.
                        CYCLE                       !  to end of sources-loop
                        
                    ELSEIF(F3 .GT. 0 .AND. ICTL35(F3,V) .GE. ADDPS) THEN
                        IDX = ICTL35( F3,V ) - ADDPS
                        IF( PKTTYP == 'MACT' ) CALL CHECK_SRC_TYPE
                        CALL SETSOURCE_CONTROL_INDEX
                        IF( IDX /= 0 ) DATSPFLAG = .TRUE.
                        CYCLE                       !  to end of sources-loop
                        
                    ELSEIF(F2 .GT. 0 .AND. ICTL34(F2,V) .GE. ADDPS) THEN
                        IDX = ICTL34( F2,V ) - ADDPS
                        IF( PKTTYP == 'MACT' ) CALL CHECK_SRC_TYPE
                        CALL SETSOURCE_CONTROL_INDEX
                        IF( IDX /= 0 ) DATSPFLAG = .TRUE.
                        CYCLE                       !  to end of sources-loop
                        
                    ELSEIF(F1 .GT. 0 .AND. ICTL33(F1,V) .GE. ADDPS) THEN
                        IDX = ICTL33( F1,V ) - ADDPS
                        IF( PKTTYP == 'MACT' ) CALL CHECK_SRC_TYPE
                        CALL SETSOURCE_CONTROL_INDEX
                        IF( IDX /= 0 ) DATSPFLAG = .TRUE.
                        CYCLE                       !  to end of sources-loop
                        
                    ELSEIF(F0 .GT. 0 .AND. ICTL32(F0,V) .GE. ADDPS) THEN
                        IDX = ICTL32( F0,V ) - ADDPS
                        IF( PKTTYP == 'MACT' ) CALL CHECK_SRC_TYPE
                        CALL SETSOURCE_CONTROL_INDEX
                        IF( IDX /= 0 ) DATSPFLAG = .TRUE.
                        CYCLE                       !  to end of sources-loop
                        
                    END IF

C.................  Try for any FIPS code, SCC, & MACT match; then
C                           any FIPS code & MACT match; then
C                           any Cy/st code, SCC, & MACT match; then
C                           any Cy/st code & MACT match; then
C                           any SCC & MACT match; then
C                           any MACT match 
                    IF( F5 .GT. 0 .AND. ICTL37(F5,V) .NE. IMISS3 ) THEN
                        IDX = ICTL37( F5,V )
                        IF( PKTTYP == 'MACT' ) CALL CHECK_SRC_TYPE
                        CALL SETSOURCE_CONTROL_INDEX
                        CYCLE                       !  to end of sources-loop
                        
                    ELSEIF(F4.GT. 0 .AND. ICTL36(F4,V) .NE. IMISS3) THEN
                        IDX = ICTL36( F4,V )
                        IF( PKTTYP == 'MACT' ) CALL CHECK_SRC_TYPE
                        CALL SETSOURCE_CONTROL_INDEX
                        CYCLE                       !  to end of sources-loop
                        
                    ELSEIF(F3.GT. 0 .AND. ICTL35(F3,V) .NE. IMISS3) THEN
                        IDX = ICTL35( F3,V )
                        IF( PKTTYP == 'MACT' ) CALL CHECK_SRC_TYPE
                        CALL SETSOURCE_CONTROL_INDEX
                        CYCLE                       !  to end of sources-loop
                        
                    ELSEIF(F2.GT. 0 .AND. ICTL34(F2,V) .NE. IMISS3) THEN
                        IDX = ICTL34( F2,V )
                        IF( PKTTYP == 'MACT' ) CALL CHECK_SRC_TYPE
                        CALL SETSOURCE_CONTROL_INDEX
                        CYCLE                       !  to end of sources-loop
                        
                    ELSEIF(F1.GT. 0 .AND. ICTL33(F1,V) .NE. IMISS3) THEN
                        IDX = ICTL33( F1,V )
                        IF( PKTTYP == 'MACT' ) CALL CHECK_SRC_TYPE
                        CALL SETSOURCE_CONTROL_INDEX
                        CYCLE                       !  to end of sources-loop
                        
                    ELSEIF(F0.GT. 0 .AND. ICTL32(F0,V) .NE. IMISS3) THEN
                        IDX = ICTL32( F0,V )
                        IF( PKTTYP == 'MACT' ) CALL CHECK_SRC_TYPE
                        CALL SETSOURCE_CONTROL_INDEX
                        CYCLE                       !  to end of sources-loop
                        
                    END IF
                END IF  ! MACT available in inventory                           

C.................  Depending on SIC/SCC switch, match by SIC first, then SCC
                IF ( SICFIRST ) THEN
                
C.....................  If SIC available in inventory, assign by SIC combinations 1st
                    ISTAT = 0
                    IF ( SICFLAG ) CALL SET_BY_SIC( ISTAT )
                    IF ( ISTAT == 1 ) CYCLE
                    
C.....................  Assign by SCC combinations 2nd
                    ISTAT = 0
                    CALL SET_BY_SCC( ISTAT )
                    IF ( ISTAT == 1 ) CYCLE
                    
                ELSE
                
C.....................  Assign by SCC combinations 1st
                    ISTAT = 0
                    CALL SET_BY_SCC( ISTAT )
                    IF ( ISTAT == 1 ) CYCLE
                    
C.....................  If SIC available in inventory, assign by SIC combinations 2nd
                    ISTAT = 0
                    IF ( SICFLAG ) CALL SET_BY_SIC( ISTAT )
                    IF ( ISTAT == 1 ) CYCLE
                END IF

C.................  Try for any FIPS code match
                F0 = FINDC( CFIP, TXCNT( 7 ), CHRT07 ) 

                IF( F0 .GT. 0 .AND. ICTL07(F0,V) .GE. ADDPS ) THEN
                    IDX = ICTL07( F0,V ) - ADDPS
                    DATSPFLAG = .TRUE.
                    CALL SETSOURCE_CONTROL_INDEX
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F0 .GT. 0 .AND. ICTL07(F0,V) .NE. IMISS3 ) THEN
                    IDX = ICTL07( F0,V )
                    CALL SETSOURCE_CONTROL_INDEX
                    CYCLE                       !  to end of sources-loop

                END IF

C.................  Try for any country/state code match 
                F0 = FINDC( CSTA, TXCNT( 4 ), CHRT04 ) 

                IF( F0 .GT. 0 .AND. ICTL04(F0,V) .GE. ADDPS ) THEN
                    IDX = ICTL04( F0,V ) - ADDPS
                    DATSPFLAG = .TRUE.
                    CALL SETSOURCE_CONTROL_INDEX
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F0 .GT. 0 .AND. ICTL04(F0,V) .NE. IMISS3 ) THEN
                    IDX = ICTL04( F0,V )
                    CALL SETSOURCE_CONTROL_INDEX
                    CYCLE                       !  to end of sources-loop

                END IF

                IF( ICTL01( V ) .NE. IMISS3 .AND. REPDEFLT ) THEN
                    IDX = ICTL01( V )
                    IF( IDX .GT. ADDPS ) THEN
                        IDX = IDX - ADDPS
                        DATSPFLAG = .TRUE.
                    END IF

                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

                    WRITE( MESG,94010 )
     &                     'NOTE: Using default ' // PKTTYP // 
     &                     ' control packet entry for:'//
     &                     CRLF() // BLANK5 // BUFFER( 1:L2 ) //
     &                     CRLF() // BLANK10 // 
     &                     ' SCC: ' // TSCC_D // ' POL: ' // EANAM( V )
                    CALL M3MESG( MESG )

                    CALL SETSOURCE_CONTROL_INDEX

                ELSEIF( ICTL01( V ) .NE. IMISS3 ) THEN
                    IDX = ICTL01( V )
                    IF( IDX .GT. ADDPS ) THEN
                        IDX = IDX - ADDPS
                        DATSPFLAG = .TRUE.
                    END IF

                    CALL SETSOURCE_CONTROL_INDEX

                END IF    !  if default profile code is available or not

            END DO        !  end loop on source, S

        IF( EFLAG ) THEN
            MESG = 'Problem assigning ' // PKTTYP // 
     &             ' controls to sources'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF 

        RETURN
c
C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94300   FORMAT( A, I2.2, A, I2.2, A )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram stores the index of the control
C               data tables for each source and pollutant
            SUBROUTINE SETSOURCE_CONTROL_INDEX

C----------------------------------------------------------------------

            SINDX( S ) = IDX

            END SUBROUTINE SETSOURCE_CONTROL_INDEX

C......................................................................
C......................................................................

C.............  This internal subprogram reports when SCC matches are
C               used instead of SIC matches.
            SUBROUTINE REPORT_SCC_USE( CHARTNUM )

C.............  Internal subprogram arguments.
            CHARACTER(2), INTENT (IN) :: CHARTNUM

C.............  Local variables
            INTEGER   L1, L2, NZ
            CHARACTER(FIPLEN3+SCCLEN3) :: BUFFER
            CHARACTER(256) MESG

C----------------------------------------------------------------------

            SELECT CASE ( CHARTNUM )
            CASE ( '09' )
                L1 = FIPLEN3+1
                L2 = LEN_TRIM( CHRT09( F11 ) )
                BUFFER = CHRT09( F11 )
                NZ = 0
            CASE ( '8C' )
                L1 = FIPLEN3+1
                L2 = LEN_TRIM( CHRT08C( F10 ) )
                BUFFER = CHRT08C( F10 )
                NZ = 2
            CASE ( '8B' )
                L1 = FIPLEN3+1
                L2 = LEN_TRIM( CHRT08B( F9 ) )
                BUFFER = CHRT08B( F9 )
                NZ = 5
            CASE ( '8A' )
                L1 = FIPLEN3+1
                L2 = LEN_TRIM( CHRT08A( F8 ) )
                BUFFER = CHRT08A( F8 )
                NZ = 7
            CASE ( '06' )
                L1 = STALEN3+1
                L2 = LEN_TRIM( CHRT06( F7 ) )
                BUFFER = CHRT06( F7 )
                NZ = 0
            CASE ( '5C' )
                L1 = STALEN3+1
                L2 = LEN_TRIM( CHRT05C( F6 ) )
                BUFFER = CHRT05C( F6 )
                NZ = 2
            CASE ( '5B' )
                L1 = STALEN3+1
                L2 = LEN_TRIM( CHRT05B( F5 ) )
                BUFFER = CHRT05B( F5 )
                NZ = 5
            CASE ( '5A' )
                L1 = STALEN3+1
                L2 = LEN_TRIM( CHRT05A( F4 ) )
                BUFFER = CHRT05A( F4 )
                NZ = 7
            CASE ( '03' )
                L1 = 1
                L2 = LEN_TRIM( CHRT03( F3 ) )
                BUFFER = CHRT03( F3 )
                NZ = 0
            CASE ( '2C' )
                L1 = 1
                L2 = LEN_TRIM( CHRT02C( F2 ) )
                BUFFER = CHRT02C( F2 )
                NZ = 2
            CASE ( '2B' )
                L1 = 1
                L2 = LEN_TRIM( CHRT02B( F1 ) )
                BUFFER = CHRT02B( F1 )
                NZ = 5
            CASE ( '2A' )
                L1 = 1
                L2 = LEN_TRIM( CHRT02A( F0 ) )
                BUFFER = CHRT02A( F0 )
                NZ = 7

            CASE DEFAULT
                MESG = 'INTERNAL ERROR: Case "'// TRIM( CHARTNUM )//
     &                 '" unknown in internal subprogram REPORT_SCC_USE'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

            END SELECT

            IF( CATEGORY .NE. 'POINT' .AND. NZ .GT. 0 ) NZ = NZ + 1

            WRITE( WDEV, 94010 ) 'WARNING: SCC assignment ' //
     &             'from '// TRIM( PKTTYP )// ' entry', IDX,
     &             'for "' // TRIM( PNAM ) // '" using SCC '// 
     &             BUFFER( L1:L2 ) // REPEAT( '0',NZ )

C----------------------  FORMAT  STATEMENTS   ------------------------

C...........   Subprogram Internal buffering formats........ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE REPORT_SCC_USE
            
C......................................................................
C......................................................................

C.............  This internal subprogram checks the source type associated
C               with a MACT packet entry.
            SUBROUTINE CHECK_SRC_TYPE

C.............  This module contains the control packet data and control matrices
            USE MODCNTRL, ONLY: CMACSRCTYP

C----------------------------------------------------------------------

            IF( CMACSRCTYP( IDX ) /= '00' ) THEN
                IF( CMACSRCTYP( IDX ) /= CSTYP ) IDX = 0
            END IF

            END SUBROUTINE CHECK_SRC_TYPE

C......................................................................
C......................................................................

C.............  This internal subprogram assigns the cross-reference
C               information by SIC and all combos. All variables are inherited.
            SUBROUTINE SET_BY_SIC( ISTAT )
            
            INTEGER :: ISTAT
            
C----------------------------------------------------------------------

C.............  Try for pollutant-specific  FIPS code & SIC matches; then
C                  pollutant-specific Cy/st code & SIC matches; then
C                  pollutant-specific SIC matches
            F5 = FINDC( CFIPSIC , TXCNT( 31 ), CHRT31 ) 
            F4 = FINDC( CFIPSICL, TXCNT( 30 ), CHRT30 ) 
            F3 = FINDC( CSTASIC , TXCNT( 29 ), CHRT29 ) 
            F2 = FINDC( CSTASICL, TXCNT( 28 ), CHRT28 ) 
            F1 = FINDC( CSIC    , TXCNT( 27 ), CHRT27 ) 
            F0 = FINDC( CSICL   , TXCNT( 26 ), CHRT26 ) 

            IF( F5 .GT. 0 .AND. ICTL31(F5,V) .GE. ADDPS ) THEN
                IDX = ICTL31( F5,V ) - ADDPS
                CALL SETSOURCE_CONTROL_INDEX
                DATSPFLAG = .TRUE.
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF(F4 .GT. 0 .AND. ICTL30(F4,V) .GE. ADDPS) THEN
                IDX = ICTL30( F4,V ) - ADDPS
                CALL SETSOURCE_CONTROL_INDEX
                DATSPFLAG = .TRUE.
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF(F3 .GT. 0 .AND. ICTL29(F3,V) .GE. ADDPS) THEN
                IDX = ICTL29( F3,V ) - ADDPS
                CALL SETSOURCE_CONTROL_INDEX
                DATSPFLAG = .TRUE.
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF(F2 .GT. 0 .AND. ICTL28(F2,V) .GE. ADDPS) THEN
                IDX = ICTL28( F2,V ) - ADDPS
                CALL SETSOURCE_CONTROL_INDEX
                DATSPFLAG = .TRUE.
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF(F1 .GT. 0 .AND. ICTL27(F1,V) .GE. ADDPS) THEN
                IDX = ICTL27( F1,V ) - ADDPS
                CALL SETSOURCE_CONTROL_INDEX
                DATSPFLAG = .TRUE.
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF(F0 .GT. 0 .AND. ICTL26(F0,V) .GE. ADDPS) THEN
                IDX = ICTL26( F0,V ) - ADDPS
                CALL SETSOURCE_CONTROL_INDEX
                DATSPFLAG = .TRUE.
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            END IF

C...........  Try for any FIPS code & SIC matches; then
C                   any Cy/st code & SIC matches; then
C                   any SIC matches; then
            IF( F5.GT. 0 .AND. ICTL31(F5,V) .NE. IMISS3 ) THEN
                IDX = ICTL31( F5,V )
                CALL SETSOURCE_CONTROL_INDEX
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF(F4.GT. 0 .AND. ICTL30(F4,V) .NE. IMISS3) THEN
                IDX = ICTL30( F4,V )
                CALL SETSOURCE_CONTROL_INDEX
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF(F3.GT. 0 .AND. ICTL29(F3,V) .NE. IMISS3) THEN
                IDX = ICTL29( F3,V )
                CALL SETSOURCE_CONTROL_INDEX
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF(F2.GT. 0 .AND. ICTL28(F2,V) .NE. IMISS3) THEN
                IDX = ICTL28( F2,V )
                CALL SETSOURCE_CONTROL_INDEX
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF(F1.GT. 0 .AND. ICTL27(F1,V) .NE. IMISS3) THEN
                IDX = ICTL27( F1,V )
                CALL SETSOURCE_CONTROL_INDEX
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF(F0.GT. 0 .AND. ICTL26(F0,V) .NE. IMISS3) THEN
                IDX = ICTL26( F0,V )
                CALL SETSOURCE_CONTROL_INDEX
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            END IF

            END SUBROUTINE SET_BY_SIC
            
C......................................................................
C......................................................................

C.............  This internal subprogram assigns the cross-reference
C               information by SCC and all combos. All variables are inherited.

            SUBROUTINE SET_BY_SCC( ISTAT )
            
            INTEGER :: ISTAT
            
C----------------------------------------------------------------------

C.............  Try for pollutant-specific FIPS code & SCC matches; then
C                       pollutant-specific Cy/st code & SCC matches; then
C                       pollutant-specific SCC matches

            F11= FINDC( CFIPS_D, TXCNT( 9 ), CHRT09 ) 
            F10= FINDC( CFIPS_C, TXCNT( 25), CHRT08C ) 
            F9 = FINDC( CFIPS_B, TXCNT( 24), CHRT08B ) 
            F8 = FINDC( CFIPS_A, TXCNT( 23), CHRT08A ) 
            F7 = FINDC( CSTAS_D, TXCNT( 6 ), CHRT06 ) 
            F6 = FINDC( CSTAS_C, TXCNT( 22), CHRT05C ) 
            F5 = FINDC( CSTAS_B, TXCNT( 21), CHRT05B ) 
            F4 = FINDC( CSTAS_A, TXCNT( 20), CHRT05A ) 
            F3 = FINDC( TSCC_D , TXCNT( 3 ), CHRT03 ) 
            F2 = FINDC( TSCC_C , TXCNT( 19), CHRT02C ) 
            F1 = FINDC( TSCC_B , TXCNT( 18), CHRT02B ) 
            F0 = FINDC( TSCC_A , TXCNT( 17), CHRT02A ) 

            IF( F11 .GT. 0 .AND. ICTL09(F11,V) .GE. ADDPS ) THEN
                IDX = ICTL09( F11,V ) - ADDPS
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '09' )
                DATSPFLAG = .TRUE.
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F10 .GT. 0 .AND. ICTL08C(F10,V) .GE. ADDPS )THEN
                IDX = ICTL08C( F10,V ) - ADDPS
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '8C' )
                DATSPFLAG = .TRUE.
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F9 .GT. 0 .AND. ICTL08B(F9,V) .GE. ADDPS ) THEN
                IDX = ICTL08B( F9,V ) - ADDPS
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '8B' )
                DATSPFLAG = .TRUE.
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F8 .GT. 0 .AND. ICTL08A(F8,V) .GE. ADDPS ) THEN
                IDX = ICTL08A( F8,V ) - ADDPS
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '8A' )
                DATSPFLAG = .TRUE.
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F7 .GT. 0 .AND. ICTL06(F7,V) .GE. ADDPS ) THEN
                IDX = ICTL06( F7,V ) - ADDPS
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '06' )
                DATSPFLAG = .TRUE.
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F6 .GT. 0 .AND. ICTL05C(F6,V) .GE. ADDPS ) THEN
                IDX = ICTL05C( F6,V ) - ADDPS
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '5C' )
                DATSPFLAG = .TRUE.
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F5 .GT. 0 .AND. ICTL05B(F5,V) .GE. ADDPS ) THEN
                IDX = ICTL05B( F5,V ) - ADDPS
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '5B' )
                DATSPFLAG = .TRUE.
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F4 .GT. 0 .AND. ICTL05A(F4,V) .GE. ADDPS ) THEN
                IDX = ICTL05A( F4,V ) - ADDPS
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '5A' )
                DATSPFLAG = .TRUE.
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F3 .GT. 0 .AND. ICTL03(F3,V) .GE. ADDPS ) THEN
                IDX = ICTL03( F3,V ) - ADDPS
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '03' )
                DATSPFLAG = .TRUE.
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F2 .GT. 0 .AND. ICTL02C(F2,V) .GE. ADDPS ) THEN
                IDX = ICTL02C( F2,V ) - ADDPS
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '2C' )
                DATSPFLAG = .TRUE.
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F1 .GT. 0 .AND. ICTL02B(F1,V) .GE. ADDPS ) THEN
                IDX = ICTL02B( F1,V ) - ADDPS
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '2B' )
                DATSPFLAG = .TRUE.
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F0 .GT. 0 .AND. ICTL02A(F0,V) .GE. ADDPS ) THEN
                IDX = ICTL02A( F0,V ) - ADDPS
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '2A' )
                DATSPFLAG = .TRUE.
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            END IF

C.............  Try for any FIPS code & SCC matches; then
C                       any Cy/st code & SCC matches; then
C                       any SCC matches; then

            IF( F11 .GT. 0 .AND. ICTL09(F11,V) .NE. IMISS3 ) THEN
                IDX = ICTL09( F11,V )
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '09' )
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F10 .GT. 0 .AND. ICTL08C(F10,V) .NE. IMISS3)THEN
                IDX = ICTL08C( F10,V )
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '8C' )
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F9 .GT. 0 .AND. ICTL08B(F9,V) .NE. IMISS3 ) THEN
                IDX = ICTL08B( F9,V )
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '8B' )
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F8 .GT. 0 .AND. ICTL08A(F8,V) .NE. IMISS3 ) THEN
                IDX = ICTL08A( F8,V )
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '8A' )
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F7 .GT. 0 .AND. ICTL06(F7,V) .NE. IMISS3 ) THEN
                IDX = ICTL06( F7,V )
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '06' )
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F6 .GT. 0 .AND. ICTL05C(F6,V) .NE. IMISS3 ) THEN
                IDX = ICTL05C( F6,V )
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '5C' )
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F5 .GT. 0 .AND. ICTL05B(F5,V) .NE. IMISS3 ) THEN
                IDX = ICTL05B( F5,V )
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '5B' )
                ISTAT = 1
                RETURN                      !  to end of sources-loop
                
            ELSEIF( F4 .GT. 0 .AND. ICTL05A(F4,V) .NE. IMISS3 ) THEN
                IDX = ICTL05A( F4,V )
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '5A' )
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F3 .GT. 0 .AND. ICTL03(F3,V) .NE. IMISS3 ) THEN
                IDX = ICTL03( F3,V )
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '03' )
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F2 .GT. 0 .AND. ICTL02C(F2,V) .NE. IMISS3 ) THEN
                IDX = ICTL02C( F2,V )
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '2C' )
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F1 .GT. 0 .AND. ICTL02B(F1,V) .NE. IMISS3 ) THEN
                IDX = ICTL02B( F1,V )
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '2B' )
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            ELSEIF( F0 .GT. 0 .AND. ICTL02A(F0,V) .NE. IMISS3 ) THEN
                IDX = ICTL02A( F0,V )
                CALL SETSOURCE_CONTROL_INDEX
                IF ( SICXREF ) CALL REPORT_SCC_USE( '2A' )
                ISTAT = 1
                RETURN                      !  to end of sources-loop

            END IF
            
            END SUBROUTINE SET_BY_SCC

        END SUBROUTINE ASGNCNTL
