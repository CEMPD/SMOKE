
        SUBROUTINE ASGNSPRO( MASSOUT, MOLEOUT, REPORT, NSRCIN, PNAM, 
     &                       MASSMATX, MOLEMATX )

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
C...........   This module contains the source arrays
        USE MODSOURC

C...........   This module contains the cross-reference tables
        USE MODXREF

C...........   This module contains the speciation profile tables
        USE MODSPRO

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        LOGICAL         ENVYN
        INTEGER         FINDC
        INTEGER         INDEX1

        EXTERNAL CRLF, ENVYN, FINDC, INDEX1

C.........  SUBROUTINE ARGUMENTS
        LOGICAL     , INTENT    (IN) :: MASSOUT        ! true: create mass-based
        LOGICAL     , INTENT    (IN) :: MOLEOUT        ! true: create mole-based
        LOGICAL     , INTENT    (IN) :: REPORT         ! true: rep defaults
        INTEGER     , INTENT    (IN) :: NSRCIN         ! number of sources
        CHARACTER(*), INTENT    (IN) :: PNAM           ! pol name of interest
        REAL        , INTENT(IN OUT) :: MASSMATX( NSRCIN,* )! mass spec matx
        REAL        , INTENT(IN OUT) :: MOLEMATX( NSRCIN,* )! mole spec matx

C.........  Other local variables
        INTEGER          L2, LV, S, V    !  counters and indices

        INTEGER          F0, F1, F2, F3, F4, F5, F6  ! tmp find indices
        INTEGER          IOS         ! i/o status
        INTEGER          NCOUT       ! no. output source chars for mesgs

        REAL             CNVFAC      ! tmp pol-to-pol conversion factor

        LOGICAL       :: EFLAG    = .FALSE.
        LOGICAL, SAVE :: FIRSTIME = .TRUE.
        LOGICAL, SAVE :: REPDEFLT = .TRUE.

        CHARACTER*8               FMTFIP   ! format for writing FIPS code
        CHARACTER*300             BUFFER   ! source fields buffer
        CHARACTER*300             MESG     ! message buffer
        CHARACTER(LEN=FIPLEN3)    CFIP     ! tmp (character) FIPS code
        CHARACTER(LEN=FPLLEN3) :: CFIPPLT = ' ' ! tmp FIPS code // plant id
        CHARACTER(LEN=FPSLEN3)    CFIPSCC  ! tmp FIPS code // SCC
        CHARACTER(LEN=FPSLEN3)    CFIPSL   ! tmp FIPS code // left SCC
        CHARACTER(LEN=STALEN3)    CSTA     ! tmp Country/state code
        CHARACTER(LEN=STSLEN3)    CSTASCC  ! tmp Country/state code // SCC
        CHARACTER(LEN=STSLEN3)    CSTASL   ! tmp Country/state code // left SCC
        CHARACTER(LEN=SRCLEN3)    CSRC     ! tmp source chars string
        CHARACTER(LEN=SS0LEN3) :: CSSC0 = ' '! tmp FIPS // Plant // SCC
        CHARACTER(LEN=SS1LEN3) :: CSSC1 = ' '! tmp src chars -> char1 // SCC
        CHARACTER(LEN=SS2LEN3) :: CSSC2 = ' '! tmp src chars -> char2 // SCC
        CHARACTER(LEN=SS3LEN3) :: CSSC3 = ' '! tmp src chars -> char3 // SCC
        CHARACTER(LEN=SS4LEN3) :: CSSC4 = ' '! tmp src chars -> char4 // SCC
        CHARACTER(LEN=SS5LEN3) :: CSSC5 = ' '! tmp src chars -> char5 // SCC
        CHARACTER(LEN=SCCLEN3)    TSCC     ! tmp 10-digit SCC
        CHARACTER(LEN=SCCLEN3)    TSCCL    ! tmp left digits of TSCC
        CHARACTER(LEN=SPNLEN3)    SPCODE   ! tmp speciation profile code

        CHARACTER*16 :: PROGNAME = 'ASGNSPRO' ! program name

C***********************************************************************
C   begin body of subroutine ASGNSPRO

C.........  For first time routine is called in all cases,
        IF( FIRSTIME ) THEN

C.............  Retrieve environment variables
            MESG = 'Switch for reporting default speciation profiles'
            REPDEFLT = ENVYN ( 'REPORT_DEFAULTS', MESG, .TRUE., IOS )

C.............  Set up format for creating character FIPS code for non-point
            IF( CATEGORY .NE. 'POINT' ) THEN
                WRITE( FMTFIP, 94300 ) '(I', FIPLEN3, '.', FIPLEN3, ')'
            ENDIF

            FIRSTIME = .FALSE.

        ENDIF

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

C.........  Find index in complete list of pollutants and set length of name
        V  = INDEX1( PNAM, NIPOL, EINAM ) 
        LV = LEN_TRIM( EINAM( V ) )

        DO S = 1, NSRCIN

C.............  Create selection 
            SELECT CASE ( CATEGORY )

            CASE ( 'AREA' ) 

                CSRC    = CSOURC( S )
                CFIP    = CSRC( 1:FIPLEN3 )
                CSTA    = CFIP( 1:STALEN3 )
                TSCC    = CSCC( S )
                TSCCL   = TSCC( 1:LSCCEND )
                CFIPSCC = CFIP // TSCC
                CFIPSL  = CFIP // TSCCL 
                CSTASCC = CSTA // TSCC
                CSTASL  = CSTA // TSCCL 

            CASE ( 'MOBILE' )

                CSRC    = CSOURC( S )
                CFIP    = CSRC( 1:FIPLEN3 )
                CSTA    = CFIP( 1:STALEN3 )
c note: add here when ready for mobile

            CASE ( 'POINT' )

                CSRC    = CSOURC( S )
                CFIP    = CSRC( 1:FIPLEN3 )
                CSTA    = CSRC( 1:STALEN3 )
                TSCC    = CSCC( S )
                CSSC5   = CSRC( 1:PTENDL3( 7 ) ) // TSCC
                CSSC4   = CSRC( 1:PTENDL3( 6 ) ) // TSCC
                CSSC3   = CSRC( 1:PTENDL3( 5 ) ) // TSCC
                CSSC2   = CSRC( 1:PTENDL3( 4 ) ) // TSCC
                CSSC1   = CSRC( 1:PTENDL3( 3 ) ) // TSCC
                CSSC0   = CSRC( 1:PTENDL3( 2 ) ) // TSCC
                CFIPPLT = CSRC( 1:PTENDL3( 2 ) )
                TSCCL   = TSCC( 1:LSCCEND )
                CFIPSCC = CFIP // TSCC
                CFIPSL  = CFIP // TSCCL
                CSTASCC = CSTA // TSCC
                CSTASL  = CSTA // TSCCL
                    
            CASE DEFAULT

            END SELECT

C.........................................................................
C.............  Initialize speciation matrices using pollutant-to-pollutant
C               conversion factors, if they exist
C.........................................................................

C.............  Screen for pollutant-to-pollutant conversion factors by checking
C               if they have been allocated
            IF( ALLOCATED( CNVRT03 ) ) THEN

C.................  Try for pollutant-specific FIPS code & SCC match; then
C                           pollutant-specific Cy/st code & SCC match; then
C                           pollutant-specific SCC match

                F3 = FINDC( CFIPSCC, NCNV3, CNVRT03 ) 
                F2 = FINDC( CSTASCC, NCNV2, CNVRT02 ) 
                F1 = FINDC( TSCC   , NCNV1, CNVRT01 ) 

        	IF( F3 .GT. 0 .AND. CNVFC03(F3,V) .NE. AMISS3 ) THEN
                    CNVFAC = CNVFC03( F3,V )

        	ELSE IF( F2 .GT. 0 .AND. CNVFC02(F2,V) .NE. AMISS3 ) THEN
                    CNVFAC = CNVFC02( F2,V )

        	ELSE IF( F1 .GT. 0 .AND. CNVFC01(F1,V) .NE. AMISS3 ) THEN
                    CNVFAC = CNVFC01( F1,V )

C.................  CNVFC00( V ) will equal 1.0 if it has not been set, so 
C                   there is no need for error checking
        	ELSE
                    CNVFAC = CNVFC00( V )

        	END IF

C.............  If they don't exist, simply set the conversion factor to one
            ELSE

                CNVFAC = 1.

            END IF

C.........................................................................
C.............  Now find and apply speciation profiles data 
C.........................................................................

C.............  In the tables used in the following heirarchy, all cross-
C               reference entries are by definition, pollutant- specific.  
C               The cross-reference tables (e.g,, CHRT02 come from MODXREF)

C.............  Try for pollutant-specific CHAR5 non-blank// SCC match; then
C                       pollutant-specific CHAR4 non-blank// SCC match; then
C                       pollutant-specific CHAR3 non-blank// SCC match; then
C                       pollutant-specific CHAR2 non-blank// SCC match; then
C                       pollutant-specific CHAR1 non-blank// SCC match; then
C                       pollutant-specific PLANT non-blank// SCC match; then
C                       pollutant-specific PLANT non-blank       match

            F6 = FINDC( CSSC5  , TXCNT( 16 ), CHRT16 ) 
            F5 = FINDC( CSSC4  , TXCNT( 15 ), CHRT15 ) 
            F4 = FINDC( CSSC3  , TXCNT( 14 ), CHRT14 ) 
            F3 = FINDC( CSSC2  , TXCNT( 13 ), CHRT13 ) 
            F2 = FINDC( CSSC1  , TXCNT( 12 ), CHRT12 ) 
            F1 = FINDC( CSSC0  , TXCNT( 11 ), CHRT11 ) 
            F0 = FINDC( CFIPPLT, TXCNT( 10 ), CHRT10 ) 

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

C.............  Try for pollutant-specific FIPS code & SCC match; then
C                       pollutant-specific FIPS code & left SCC match; then
C                       pollutant-specific Cy/st code & SCC match; then
C                       pollutant-specific Cy/st code & left SCC match; then
C                       pollutant-specific SCC match; then
C                       pollutant-specific left SCC match

            F5 = FINDC( CFIPSCC, TXCNT( 9 ), CHRT09 ) 
            F4 = FINDC( CFIPSL , TXCNT( 8 ), CHRT08 ) 
            F3 = FINDC( CSTASCC, TXCNT( 6 ), CHRT06 ) 
            F2 = FINDC( CSTASL , TXCNT( 5 ), CHRT05 ) 
            F1 = FINDC( TSCC   , TXCNT( 3 ), CHRT03 ) 
            F0 = FINDC( TSCCL  , TXCNT( 2 ), CHRT02 ) 

            IF( F5 .GT. 0 .AND. CSPT09(F5,V) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT09( F5,V ) 
                CALL SETSOURCE_SMATS
                CYCLE                       !  to end of sources-loop

            ELSEIF( F4 .GT. 0 .AND. CSPT08(F4,V) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT08( F4,V ) 
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

            ELSEIF( F1 .GT. 0 .AND. CSPT03(F1,V) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT03( F1,V ) 
                CALL SETSOURCE_SMATS
                CYCLE                       !  to end of sources-loop

            ELSEIF( F0 .GT. 0 .AND. CSPT02(F0,V) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT02( F0,V ) 
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

            IF( CSPT01( V ) .NE. EMCMISS3 .AND. 
     &          REPDEFLT .AND. REPORT           ) THEN
                SPCODE = CSPT01( V )
                    
                CALL FMTCSRC( CSRC, NCOUT, BUFFER, L2 )

C.................  For default speciation profile, make sure that it has been
C                   defined for the current pollutant
                IF( SPCODE .EQ. ' ' ) THEN

                    EFLAG = .TRUE.
                    MESG = 'ERROR: Default speciation profile not ' //
     &                     'defined for POL: '// EINAM(V)(1:LV)// ','//
     &                     CRLF() // BLANK10 // 'but needed for:'//
     &                     CRLF() // BLANK10 // BUFFER( 1:L2 ) //
     &                     CRLF() // BLANK10 // 
     &                     'SCC: ' // TSCC
                    CALL M3MESG( MESG )

C.................  It profile defined, try to apply it
                ELSE

                    MESG = 'NOTE: Using default speciation profile "' //
     &                     SPCODE // '" for:'//
     &                     CRLF() // BLANK10 // BUFFER( 1:L2 ) //
     &                     CRLF() // BLANK10 // 
     &                     'SCC: ' // TSCC // ' POL: ' // EINAM( V )
                    CALL M3MESG( MESG )

                    CALL SETSOURCE_SMATS

                 END IF

            ELSEIF( CSPT01( V ) .NE. EMCMISS3 ) THEN
                SPCODE = CSPT01( V )
                CALL SETSOURCE_SMATS

            ELSE
                CALL FMTCSRC( CSRC, NCOUT, BUFFER, L2 )

                EFLAG = .TRUE.
                MESG = 'ERROR: No speciation cross-reference ' //
     &                 'available (and no default) for:' //
     &                 CRLF() // BLANK10 // BUFFER( 1:L2 ) //
     &                 CRLF() // BLANK10 // 
     &                 'SCC: ' // TSCC // ' POL: ' // EINAM( V )

                CALL M3MESG( MESG )

            END IF    !  if default profile code is available or not

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

C.............  This internal subprogram searches for the speciation profile
C               code in the abriged list (from MODSPRO) and if found, applies
C               the speciation factors for all species in that profile to 
C               the speciation matrices.
C.............  Most variables are defined through host association.
            SUBROUTINE SETSOURCE_SMATS

C.............  Local variables
            INTEGER   I, J, K, N  ! indices and counters

            INTEGER   ITBL        ! position in full table of current profile
            INTEGER   NTBL        ! number of species of current profile

C----------------------------------------------------------------------

            K = MAX( FINDC( SPCODE, NSPROF, SPROFN ), 0 )

            IF( K .EQ. 0 ) THEN

                CALL FMTCSRC( CSRC, NCOUT, BUFFER, L2 )

                EFLAG = .TRUE.
                MESG = 'ERROR: Speciation profile "' // SPCODE // 
     &                 '" is not in profiles, but it was assigned' //
     &                 CRLF() // BLANK10 // 'to source:' //
     &                 CRLF() // BLANK10 // BUFFER( 1:L2 ) //
     &                 CRLF() // BLANK10 // 
     &                 'SCC: ' // TSCC // ' POL: ' // EINAM( V )
                CALL M3MESG( MESG )

            ELSE

C.................  Get indices to full speciation table
                ITBL = IDXSPRO ( K )
                NTBL = NSPECIES( K )

                I = ITBL - 1
                DO N = 1, NTBL

                    I = I + 1
                    J = IDXSSPEC( K,N )

                    IF( MASSOUT ) THEN
                        MASSMATX( S,J )= CNVFAC * MASSFACT( I )
                    END IF

                    IF( MOLEOUT ) THEN
                        MOLEMATX( S,J )= CNVFAC * MOLEFACT( I )
                    END IF

                END DO 

            END IF

C------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

            END SUBROUTINE SETSOURCE_SMATS

        END SUBROUTINE ASGNSPRO
