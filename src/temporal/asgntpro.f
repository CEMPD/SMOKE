
        SUBROUTINE ASGNTPRO( CATEGORY, NSRC, NPOL, NIPOL, PNAM, EINAM,
     &                       TREFFMT )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      For each source and pollutant, find the most specific temporal profile
C      that applies to that source. Do this using the grouped tables of
C      temporal cross references from RDTREF.  The hierarchical order is
C      defined in this subroutine, and can be determined from the in-source
C      comments below. Once a profile code has been identified, search for this
C      code in the temporal profile tables (from RDTPROF) and save the index
C      to these tables for each source and pollutant.
C
C  PRECONDITIONS REQUIRED:
C Expects tables to have IMISS3 where they are undefined
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 1/99 by M. Houyoux
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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
C***************************************************************************

C...........   MODULES for public variables   
C...........   This module contains the source ararys
        USE MODSOURC

C...........   This module contains the cross-reference tables
        USE MODXREF

C...........   This module contains the temporal profile tables
        USE MODTPRO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        LOGICAL         ENVYN
        INTEGER         FIND1
        INTEGER         FINDC
        INTEGER         INDEX1

        EXTERNAL CRLF, ENVYN, FIND1, FINDC, INDEX1

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: CATEGORY       ! source category
        INTEGER     , INTENT (IN) :: NSRC           ! number of sources
        INTEGER     , INTENT (IN) :: NPOL           ! number of pols in group
        INTEGER     , INTENT (IN) :: NIPOL          ! total number of pollutants
        CHARACTER(*), INTENT (IN) :: PNAM ( NPOL )  ! pollutant names in group
        CHARACTER(*), INTENT (IN) :: EINAM( NIPOL ) ! all pollutant names
        CHARACTER(*), INTENT (IN) :: TREFFMT        ! temporal xref format

C.........  Other local variables
        INTEGER          I, J, L2, S, V    !  counters and indices

        INTEGER          F0, F1, F2, F3, F4, F5, F6  ! tmp find indices
        INTEGER          MREF    !  tmp monthly profile code
        INTEGER          WREF    !  tmp weekly  profile code
        INTEGER          DREF    !  tmp diurnal profile code

        LOGICAL       :: EFLAG    = .FALSE.
        LOGICAL, SAVE :: FIRSTIME = .TRUE.
        LOGICAL          MFLAG   !  use monthly profiles
        LOGICAL          WFLAG   !  use weekly  profiles
        LOGICAL, SAVE :: REPDEFLT = .TRUE.

        CHARACTER*8            FMTFIP   ! format for writing FIPS code
        CHARACTER*300          BUFFER   ! source fields buffer
        CHARACTER*300          LINE     ! line buffer
        CHARACTER*300          MESG     ! message buffer
        CHARACTER(LEN=STALEN3) CSTA     ! tmp Country/state code
        CHARACTER(LEN=STSLEN3) CSTASCC  ! tmp Country/state code // SCC
        CHARACTER(LEN=STSLEN3) CSTASL5  ! tmp Country/state code // left SCC
        CHARACTER(LEN=SCLLEN3) TSCCL5   ! tmp left digits of TSCC
        CHARACTER(LEN=SRCLEN3) CSRC     ! tmp source chars string
        CHARACTER(LEN=SS0LEN3) CSSC0    ! tmp FIPS // Plant // SCC
        CHARACTER(LEN=SS1LEN3) CSSC1    ! tmp source chars through char1 // SCC
        CHARACTER(LEN=SS2LEN3) CSSC2    ! tmp source chars through char2 // SCC
        CHARACTER(LEN=SS3LEN3) CSSC3    ! tmp source chars through char3 // SCC
        CHARACTER(LEN=SS4LEN3) CSSC4    ! tmp source chars through char4 // SCC
        CHARACTER(LEN=SS5LEN3) CSSC5    ! tmp source chars through char5 // SCC
        CHARACTER(LEN=FIPLEN3) CFIP     ! tmp (character) FIPS code
        CHARACTER(LEN=FPLLEN3) CFIPPLT  ! tmp FIPS code // plant id
        CHARACTER(LEN=FPSLEN3) CFIPSCC  ! tmp FIPS code // SCC
        CHARACTER(LEN=FPSLEN3) CFIPSL5  ! tmp FIPS code // left SCC
        CHARACTER(LEN=SCCLEN3) TSCC     ! tmp 10-digit SCC

        CHARACTER*16 :: PROGNAME = 'ASGNTPRO' ! program name

C***********************************************************************
C   begin body of subroutine ASGNTPRO

C.........  For list-formatted temporal cross-reference (one entry per source)
C           from EMS-95 files, the profiles are not applied per pollutant.  So,
C           we can set these for the first group of pollutants used when
C           calling this subroutine and then use them for all pollutants.
        IF( FIRSTIME .AND. TREFFMT .EQ. 'LIST' ) THEN

C.............  Set for first pollutant in group  
            J = 1
            DO S = 1, NSRC

C.................  Set MFLAG to true for using monthly temporal adjustments
                MFLAG = ( MOD( TPFLAG( S ), MTPRFAC ) .EQ. 0 )

C.................  Set WFLAG to trur for using weekly temporal adjustments
                WFLAG = ( MOD( TPFLAG( S ), WTPRFAC ) .EQ. 0 .OR.
     &                    MOD( TPFLAG( S ), WDTPFAC ) .EQ. 0      )

                MREF = MPRNA( S )
                WREF = WPRNA( S )
                DREF = DPRNA( S )
                CALL SETSOURCE_TPROFS

            ENDDO

C.............  Set for remaining pollutants in group  
            DO J = 2, NPOL
                MDEX( :,J ) = MDEX( :,1 )
                WDEX( :,J ) = WDEX( :,1 )
                DDEX( :,J ) = DDEX( :,1 )
                EDEX( :,J ) = EDEX( :,1 )
            ENDDO

       ENDIF

C.........  For first time routine is called in all cases,
        IF( FIRSTIME ) THEN

C.............  Retrieve environment variables
            MESG = 'Switch for reporting default temporal profiles'
            REPDEFLT = ENVYN ( 'REPORT_DEFAULTS', MESG, .TRUE., I )

C.............  Set up format for creating character FIPS code for non-point
            IF( CATEGORY .NE. 'POINT' ) THEN
                WRITE( FMTFIP, 94300 ) '(I', FIPLEN3, '.', FIPLEN3, ')'
            ENDIF

            FIRSTIME = .FALSE.

        ENDIF

C.........  Exit subroutine for list-formatted temporal x-ref because we
C           do not have a heirarchial application of temporal profiles
C           to worry about.
        IF( TREFFMT .EQ. 'LIST' ) RETURN

        DO J = 1, NPOL

C.............  Find index in complete list of pollutants
            V = INDEX1( PNAM( J ), NIPOL, EINAM )

            DO S = 1, NSRC

C.................  Set MFLAG to true for using monthly temporal adjustments
                MFLAG = ( MOD( TPFLAG( S ), MTPRFAC ) .EQ. 0 )

C.................  Set WFLAG to trur for using weekly temporal adjustments
                WFLAG = ( MOD( TPFLAG( S ), WTPRFAC ) .EQ. 0 .OR.
     &                    MOD( TPFLAG( S ), WDTPFAC ) .EQ. 0      )

C.................  Create selection 
                SELECT CASE ( CATEGORY )

                CASE ( 'AREA' ) 
                    WRITE( CFIP, FMTFIP ) IFIP( S )

                    CSTA    = CFIP( 1:STALEN3 )
                    TSCC    = CSCC( S )
                    TSCCL5  = TSCC( 1:SCLLEN3 )
                    CFIPSCC = CFIP // TSCC
                    CFIPSL5 = CFIP // TSCCL5
                    CSTASCC = CSTA // TSCC
                    CSTASL5 = CSTA // TSCCL5

                CASE ( 'MOBILE' )

                CASE ( 'POINT' )

                    CSRC    = CSOURC( S )
                    TSCC    = CSCC( S )
                    CSSC5   = CSRC( 1:PTENDL3( 7 ) ) // TSCC
                    CSSC4   = CSRC( 1:PTENDL3( 6 ) ) // TSCC
                    CSSC3   = CSRC( 1:PTENDL3( 5 ) ) // TSCC
                    CSSC2   = CSRC( 1:PTENDL3( 4 ) ) // TSCC
                    CSSC1   = CSRC( 1:PTENDL3( 3 ) ) // TSCC
                    CSSC0   = CSRC( 1:PTENDL3( 2 ) ) // TSCC
                    CFIPPLT = CSRC( 1:PTENDL3( 2 ) )
                    CFIP    = CSRC( 1:FIPLEN3 )
                    CSTA    = CSRC( 1:STALEN3 )
                    TSCCL5  = TSCC( 1:SCLLEN3 )
                    CFIPSCC = CFIP // TSCC
                    CFIPSL5 = CFIP // TSCCL5
                    CSTASCC = CSTA // TSCC
                    CSTASL5 = CSTA // TSCCL5
                    
                CASE DEFAULT

                END SELECT

C NOTE: what happens when TXCNT( 10 ) = 0 ?? Do I need a test for > 0 ?

C.................  In the tables used in the following heirarchy, a pollutant-
C                   specific cross-reference entry has not been use as the
C                   default for all pollutants.  So the monthly profile number
C                   tables (MPRT*) are checked to ensure the pollutant has
C                   been defined for a level of matching of interest.

C.................  Try for pollutant-specific CHAR5 non-blank// SCC match; then
C                           pollutant-specific CHAR4 non-blank// SCC match; then
C                           pollutant-specific CHAR3 non-blank// SCC match; then
C                           pollutant-specific CHAR2 non-blank// SCC match; then
C                           pollutant-specific CHAR1 non-blank// SCC match; then
C                           pollutant-specific PLANT non-blank// SCC match; then
C                           pollutant-specific PLANT non-blank       match

                F6 = FINDC( CSSC5  , TXCNT( 16 ), CHRT16 ) 
                F5 = FINDC( CSSC4  , TXCNT( 15 ), CHRT15 ) 
                F4 = FINDC( CSSC3  , TXCNT( 14 ), CHRT14 ) 
                F3 = FINDC( CSSC2  , TXCNT( 13 ), CHRT13 ) 
                F2 = FINDC( CSSC1  , TXCNT( 12 ), CHRT12 ) 
                F1 = FINDC( CSSC0  , TXCNT( 11 ), CHRT11 ) 
                F0 = FINDC( CFIPPLT, TXCNT( 10 ), CHRT10 ) 

                IF( F6 .GT. 0 .AND. DPRT16(F6,V) .GE. ADDPS ) THEN
                    MREF = MPRT16( F6,V )
                    WREF = WPRT16( F6,V )
                    DREF = DPRT16( F6,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F5 .GT. 0 .AND. DPRT15(F5,V) .GE. ADDPS ) THEN
                    MREF = MPRT15( F5,V )
                    WREF = WPRT15( F5,V )
                    DREF = DPRT15( F5,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F4 .GT. 0 .AND. DPRT14(F4,V) .GE. ADDPS ) THEN
                    MREF = MPRT14( F4,V )
                    WREF = WPRT14( F4,V )
                    DREF = DPRT14( F4,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F3 .GT. 0 .AND. DPRT13(F3,V) .GE. ADDPS ) THEN
                    MREF = MPRT13( F3,V )
                    WREF = WPRT13( F3,V )
                    DREF = DPRT13( F3,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F2 .GT. 0 .AND. DPRT12(F2,V) .GE. ADDPS ) THEN
                    MREF = MPRT12( F2,V )
                    WREF = WPRT12( F2,V )
                    DREF = DPRT12( F2,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F1 .GT. 0 .AND. DPRT11(F1,V) .GE. ADDPS ) THEN
                    MREF = MPRT11( F1,V )
                    WREF = WPRT11( F1,V )
                    DREF = DPRT11( F1,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F0 .GT. 0 .AND. DPRT10(F0,V) .GE. ADDPS ) THEN
                    MREF = MPRT10( F0,V )
                    WREF = WPRT10( F0,V )
                    DREF = DPRT10( F0,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                END IF

C.................  Try for any CHAR5 non-blank // SCC match; then
C                           any CHAR4 non-blank // SCC match; then
C                           any CHAR3 non-blank // SCC match; then
C                           any CHAR2 non-blank // SCC match; then
C                           any CHAR1 non-blank // SCC match; then
C                           any PLANT non-blank // SCC match; then
C                           any PLANT non-blank        match

                IF( F6 .GT. 0 .AND. DPRT16(F6,V) .NE. IMISS3 ) THEN
                    MREF = MPRT16( F6,V )
                    WREF = WPRT16( F6,V )
                    DREF = DPRT16( F6,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F5 .GT. 0 .AND. DPRT15(F5,V) .NE. IMISS3 ) THEN
                    MREF = MPRT15( F5,V )
                    WREF = WPRT15( F5,V )
                    DREF = DPRT15( F5,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F4 .GT. 0 .AND. DPRT14(F4,V) .NE. IMISS3 ) THEN
                    MREF = MPRT14( F4,V )
                    WREF = WPRT14( F4,V )
                    DREF = DPRT14( F4,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F3 .GT. 0 .AND. DPRT13(F3,V) .NE. IMISS3 ) THEN
                    MREF = MPRT13( F3,V )
                    WREF = WPRT13( F3,V )
                    DREF = DPRT13( F3,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F2 .GT. 0 .AND. DPRT12(F2,V) .NE. IMISS3 ) THEN
                    MREF = MPRT12( F2,V )
                    WREF = WPRT12( F2,V )
                    DREF = DPRT12( F2,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F1 .GT. 0 .AND. DPRT11(F1,V) .NE. IMISS3 ) THEN
                    MREF = MPRT11( F1,V )
                    WREF = WPRT11( F1,V )
                    DREF = DPRT11( F1,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F0 .GT. 0 .AND. DPRT10(F0,V) .NE. IMISS3 ) THEN
                    MREF = MPRT10( F0,V )
                    WREF = WPRT10( F0,V )
                    DREF = DPRT10( F0,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                END IF

C.................  Try for pollutant-specific FIPS code & SCC match; then
C                           pollutant-specific FIPS code & left SCC match; then
C                           pollutant-specific Cy/st code & SCC match; then
C                           pollutant-specific Cy/st code & left SCC match; then
C                           pollutant-specific SCC match; then
C                           pollutant-specific left SCC match

                F5 = FINDC( CFIPSCC, TXCNT( 9 ), CHRT09 ) 
                F4 = FINDC( CFIPSL5, TXCNT( 8 ), CHRT08 ) 
                F3 = FINDC( CSTASCC, TXCNT( 6 ), CHRT06 ) 
                F2 = FINDC( CSTASL5, TXCNT( 5 ), CHRT05 ) 
                F1 = FINDC( TSCC   , TXCNT( 3 ), CHRT03 ) 
                F0 = FINDC( TSCCL5 , TXCNT( 2 ), CHRT02 ) 

                IF( F5 .GT. 0 .AND. DPRT09(F5,V) .GE. ADDPS ) THEN
                    MREF = MPRT09( F5,V )
                    WREF = WPRT09( F5,V )
                    DREF = DPRT09( F5,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F4 .GT. 0 .AND. DPRT08(F4,V) .GE. ADDPS ) THEN
                    MREF = MPRT08( F4,V )
                    WREF = WPRT08( F4,V )
                    DREF = DPRT08( F4,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F3 .GT. 0 .AND. DPRT06(F3,V) .GE. ADDPS ) THEN
                    MREF = MPRT06( F3,V )
                    WREF = WPRT06( F3,V )
                    DREF = DPRT06( F3,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F2 .GT. 0 .AND. DPRT05(F2,V) .GE. ADDPS ) THEN
                    MREF = MPRT05( F2,V )
                    WREF = WPRT05( F2,V )
                    DREF = DPRT05( F2,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F1 .GT. 0 .AND. DPRT03(F1,V) .GE. ADDPS ) THEN
                    MREF = MPRT03( F1,V )
                    WREF = WPRT03( F1,V )
                    DREF = DPRT03( F1,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F0 .GT. 0 .AND. DPRT02(F0,V) .GE. ADDPS ) THEN
                    MREF = MPRT02( F0,V )
                    WREF = WPRT02( F0,V )
                    DREF = DPRT02( F0,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                END IF

C.................  Try for any FIPS code & SCC match; then
C                           any FIPS code & left SCC match; then
C                           any Cy/st code & SCC match; then
C                           any Cy/st code & left SCC match; then
C                           any SCC match; then
C                           any left SCC match

                IF( F5 .GT. 0 .AND. DPRT09(F5,V) .NE. IMISS3 ) THEN
                    MREF = MPRT09( F5,V ) 
                    WREF = WPRT09( F5,V )
                    DREF = DPRT09( F5,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F4 .GT. 0 .AND. DPRT08(F4,V) .NE. IMISS3 ) THEN
                    MREF = MPRT08( F4,V ) 
                    WREF = WPRT08( F4,V )
                    DREF = DPRT08( F4,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F3 .GT. 0 .AND. DPRT06(F3,V) .NE. IMISS3 ) THEN
                    MREF = MPRT06( F3,V ) 
                    WREF = WPRT06( F3,V )
                    DREF = DPRT06( F3,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F2 .GT. 0 .AND. DPRT05(F2,V) .NE. IMISS3 ) THEN
                    MREF = MPRT05( F2,V ) 
                    WREF = WPRT05( F2,V )
                    DREF = DPRT05( F2,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F1 .GT. 0 .AND. DPRT03(F1,V) .NE. IMISS3 ) THEN
                    MREF = MPRT03( F1,V ) 
                    WREF = WPRT03( F1,V )
                    DREF = DPRT03( F1,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F0 .GT. 0 .AND. DPRT02(F0,V) .NE. IMISS3 ) THEN
                    MREF = MPRT02( F0,V ) 
                    WREF = WPRT02( F0,V )
                    DREF = DPRT02( F0,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                END IF

C.................  Try for any FIPS code match
                F0 = FINDC( CFIP, TXCNT( 7 ), CHRT07 ) 

                IF( F0 .GT. 0 ) THEN
                    MREF = MPRT07( F0 ) 
                    WREF = WPRT07( F0 )
                    DREF = DPRT07( F0 )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop
                END IF

C.................  Try for any country/state code match (not, pol-specific)
                F0 = FINDC( CSTA, TXCNT( 4 ), CHRT04 ) 

                IF( F0 .GT. 0 ) THEN
                    MREF = MPRT04( F0 ) 
                    WREF = WPRT04( F0 )
                    DREF = DPRT04( F0 )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop
                END IF

                IF( MPRT01 .NE. IMISS3 .AND. REPDEFLT ) THEN
                    MREF = MPRT01
                    WREF = WPRT01
                    DREF = DPRT01
                    
                    CALL FMTCSRC( CSRC, 7, BUFFER, L2 )

                    WRITE( MESG,94010 )
     &                     'NOTE: Using default temporal profile for:'//
     &                     CRLF() // BLANK5 // BUFFER( 1:L2 ) //
     &                     CRLF() // BLANK10 // 
     &                     ' SCC: ' // TSCC // ' POL: ' // EINAM( V )
                    CALL M3MESG( MESG )

                    CALL SETSOURCE_TPROFS

                ELSEIF( MPRT01 .NE. IMISS3 ) THEN
                    MREF = MPRT01
                    WREF = WPRT01
                    DREF = DPRT01
                    CALL SETSOURCE_TPROFS

                ELSE
                    EFLAG = .TRUE.

                    WRITE( MESG,94010 )
     &                     'ERROR: No temporal cross-reference ' //
     &                     'available (and no default) for:' //
     &                     CRLF() // BLANK5 // BUFFER( 1:L2 )

                    CALL M3MESG( MESG )

                END IF    !  if default profile code is available or not

            END DO        !  end loop on source, S

        END DO            !  end loop on pollutant, V

        IF( EFLAG ) THEN
            MESG = 'Problem assigning temporal profiles to sources'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF 

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94300   FORMAT( A, I2.2, A, I2.2, A )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram stores the index of the temporal 
C               profile codes in the temporal profile tables for each source.
C.............  All variables are defined through host association.
            SUBROUTINE SETSOURCE_TPROFS

C----------------------------------------------------------------------

            IF( MFLAG ) THEN

                MDEX( S,J ) = MAX( FIND1( MREF, NMON, MONREF ), 0 )

                IF( MDEX( S,J ) .EQ. 0 ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                     'ERROR: Monthly profile', MREF, 
     &                     'is not in profiles, but was assigned' //
     &                     CRLF() // BLANK5 // 'to source:' //
     &                     CRLF() // BLANK5 // BUFFER( 1:L2 )
                    CALL M3MESG( MESG )

                END IF

            END IF  ! If monthly profiles are being used or not

            IF( WFLAG ) THEN

                WDEX( S,J ) = MAX( FIND1( WREF, NWEK, WEKREF ), 0 )

                IF( WDEX( S,J ) .EQ. 0 ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                     'ERROR: Weekly profile', WREF, 
     &                     'is not in profiles, but was assigned' //
     &                     CRLF() // BLANK5 // 'to source:' //
     &                     CRLF() // BLANK5 // BUFFER( 1:L2 )
                    CALL M3MESG( MESG )
                END IF

            END IF

            DDEX( S,J ) = MAX( FIND1( DREF, NDIU, DIUREF ), 0 )

            IF( NEND .GT. 0 ) THEN

                EDEX( S,J ) = MAX( FIND1( DREF, NEND, ENDREF ), 0 )

            END IF  ! If there are weekend diurnal profiles

            IF( DDEX( S,J ) .EQ. 0 .AND. EDEX( S,J ) .EQ. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &                 'ERROR: Diurnal profile', DREF, 
     &                 'is not in profiles, but was assigned' //
     &                 CRLF() // BLANK5 // 'to source:' //
     &                 CRLF() // BLANK5 // BUFFER( 1:L2 )
                CALL M3MESG( MESG )

            END IF

C------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE SETSOURCE_TPROFS

        END SUBROUTINE ASGNTPRO
