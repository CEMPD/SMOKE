
        SUBROUTINE ASGNPSI( NACT, ANAM )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C     This routine assigns the parameter scheme index (PSI) for each activity
C     to all sources and 24 hours.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 6/99 by M. Houyoux
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

C.........  This module contains the information about the source category
        USE MODINFO

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
        INTEGER     , INTENT (IN) :: NACT         ! number of activities 
        CHARACTER(*), INTENT (IN) :: ANAM( NACT ) ! activity names

C.........  Other local variables
        INTEGER          I, J, L2, S, V    !  counters and indices

        INTEGER          F0, F1, F2, F3, F4, F5, F6  ! tmp find indices
        INTEGER          EREF    !  tmp diurnal profile code

        LOGICAL       :: EFLAG    = .FALSE.
        LOGICAL, SAVE :: FIRSTIME = .TRUE.
        LOGICAL, SAVE :: REPDEFLT = .TRUE.

        CHARACTER*10           RWTFMT   ! roadway type format
        CHARACTER*300          BUFFER   ! source fields buffer
        CHARACTER*300          MESG     ! message buffer
        CHARACTER(LEN=STALEN3) CSTA     ! tmp Country/state code
        CHARACTER(LEN=STSLEN3) CSTASCC  ! tmp Country/state code // SCC
        CHARACTER(LEN=STSLEN3) CSTASL   ! tmp Country/state code // left SCC
        CHARACTER(LEN=SCCLEN3) TSCCL    ! tmp left digits of TSCC
        CHARACTER(LEN=SRCLEN3) CSRC     ! tmp source chars string
        CHARACTER(LEN=SS0LEN3) CSSC0    ! tmp FIPS // Plant // SCC
        CHARACTER(LEN=SS1LEN3) CSSC1    ! tmp source chars through char1 // SCC
        CHARACTER(LEN=SS2LEN3) CSSC2    ! tmp source chars through char2 // SCC
        CHARACTER(LEN=SS3LEN3) CSSC3    ! tmp source chars through char3 // SCC
        CHARACTER(LEN=SS4LEN3) CSSC4    ! tmp source chars through char4 // SCC
        CHARACTER(LEN=SS5LEN3) CSSC5    ! tmp source chars through char5 // SCC
        CHARACTER(LEN=FIPLEN3) CFIP     ! tmp (character) FIPS code
        CHARACTER(LEN=FPLLEN3) CFIPPLS  ! tmp FIPS code // plant id or rdtype
        CHARACTER(LEN=FPSLEN3) CFIPSCC  ! tmp FIPS code // SCC
        CHARACTER(LEN=FPSLEN3) CFIPSL   ! tmp FIPS code // left SCC
        CHARACTER(LEN=SCCLEN3) TSCC     ! tmp 10-digit SCC

        CHARACTER*16 :: PROGNAME = 'ASGNPSI' ! program name

C***********************************************************************
C   begin body of subroutine ASGNPSI

C.........  For first time routine is called in all cases,
        IF( FIRSTIME ) THEN

C.............  Retrieve environment variables
            MESG = 'Switch for reporting default emission factors'
            REPDEFLT = ENVYN ( 'REPORT_DEFAULTS', MESG, .TRUE., I )

            FIRSTIME = .FALSE.

        ENDIF

C.........  Set up roadway type format
        WRITE( RWTFMT, '("(I",I2.2")")' ) RWTLEN3

C.........  Initialize search fields to blank
        CSRC    = ' '
        TSCC    = ' '
        CSSC5   = ' '
        CSSC4   = ' '
        CSSC3   = ' '
        CSSC2   = ' '
        CSSC1   = ' '
        CSSC0   = ' '
        CFIPPLS = ' '
        CFIP    = ' '
        CSTA    = ' '
        TSCCL   = ' '
        CFIPSCC = ' '
        CFIPSL  = ' '
        CSTASCC = ' '
        CSTASL  = ' '

C.........  Loop through activities
        DO J = 1, NACT

C.............  Find index in complete list of activities
            V = INDEX1( ANAM( J ), NIACT, ACTVTY )

            DO S = 1, NSRC

C.................  Create selection 
                SELECT CASE ( CATEGORY )

                CASE ( 'AREA' ) 
                    CFIP    = CSRC( 1:FIPLEN3 )
                    CSTA    = CFIP( 1:STALEN3 )
                    TSCC    = CSCC( S )
                    TSCCL   = TSCC( 1:LSCCEND )
                    CFIPSCC = CFIP // TSCC
                    CFIPSL  = CFIP // TSCCL
                    CSTASCC = CSTA // TSCC
                    CSTASL  = CSTA // TSCCL

C.................  For PSIs, SCC is not used for matching, but roadway type
C                   and link IDs are used, so adjust fields accordingly
                CASE ( 'MOBILE' )

                    CSRC    = CSOURC( S )
                    WRITE( TSCC, RWTFMT ) IRCLAS( S )
                    CALL PADZERO( TSCC )

                    CFIP    = CSRC( 1:FIPLEN3 )
                    CSTA    = CSRC( 1:STALEN3 )
                    TSCCL   = TSCC( 1:LSCCEND )
                    CFIPSCC = CFIP // TSCC
                    CFIPSL  = CFIP // TSCCL 
                    CSTASCC = CSTA // TSCC
                    CSTASL  = CSTA // TSCCL
                    CFIPPLS = CSRC( 1:MBENDL3( 2 ) )
                    CSSC0   = CSRC( 1:MBENDL3( 2 ) ) // TSCC  ! w/ road class
                    CSSC1   = CSRC( 1:MBENDL3( 3 ) ) // TSCC  ! w/ link
                    CSSC2   = CSRC( 1:MBENDL3( 4 ) ) // TSCC  ! w/ vehicle type

                CASE ( 'POINT' )

                    CSRC    = CSOURC( S )
                    TSCC    = CSCC( S )
                    CSSC5   = CSRC( 1:PTENDL3( 7 ) ) // TSCC
                    CSSC4   = CSRC( 1:PTENDL3( 6 ) ) // TSCC
                    CSSC3   = CSRC( 1:PTENDL3( 5 ) ) // TSCC
                    CSSC2   = CSRC( 1:PTENDL3( 4 ) ) // TSCC
                    CSSC1   = CSRC( 1:PTENDL3( 3 ) ) // TSCC
                    CSSC0   = CSRC( 1:PTENDL3( 2 ) ) // TSCC
                    CFIPPLS = CSRC( 1:PTENDL3( 2 ) )
                    CFIP    = CSRC( 1:FIPLEN3 )
                    CSTA    = CSRC( 1:STALEN3 )
                    TSCCL   = TSCC( 1:LSCCEND )
                    CFIPSCC = CFIP // TSCC
                    CFIPSL  = CFIP // TSCCL 
                    CSTASCC = CSTA // TSCC
                    CSTASL  = CSTA // TSCCL
                    
                CASE DEFAULT

                END SELECT

C NOTE: what happens when TXCNT( 10 ) = 0 ?? Do we need a test for > 0 ?

C.................  In the tables used in the following heirarchy, all matches
C                   must be activity-specific.  There is no such thing as
C                   an emission factor that is not activity-specific.  The 
C                   emission factor index tables are checked to make sure
C                   that the activity was specified in the cross-reference.

C.................  Try for activity-specific CHAR5 non-blank// SCC match; then
C                           activity-specific CHAR4 non-blank// SCC match; then
C                           activity-specific CHAR3 non-blank// SCC match; then
C                           activity-specific CHAR2 non-blank// SCC match; then
C                           activity-specific CHAR1 non-blank// SCC match; then
C                           activity-specific PLANT non-blank// SCC match; then
C                           activity-specific PLANT non-blank       match

                F6 = FINDC( CSSC5  , TXCNT( 16 ), CHRT16 ) 
                F5 = FINDC( CSSC4  , TXCNT( 15 ), CHRT15 ) 
                F4 = FINDC( CSSC3  , TXCNT( 14 ), CHRT14 ) 
                F3 = FINDC( CSSC2  , TXCNT( 13 ), CHRT13 ) 
                F2 = FINDC( CSSC1  , TXCNT( 12 ), CHRT12 ) 
                F1 = FINDC( CSSC0  , TXCNT( 11 ), CHRT11 ) 
                F0 = FINDC( CFIPPLS, TXCNT( 10 ), CHRT10 ) 

                IF( F6 .GT. 0 .AND. IEFS16(F6,V) .GE. ADDPS ) THEN
                    EREF = IEFS16( F6,V ) - ADDPS
                    CALL SETSOURCE_EFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F5 .GT. 0 .AND. IEFS15(F5,V) .GE. ADDPS ) THEN
                    EREF = IEFS15( F5,V ) - ADDPS
                    CALL SETSOURCE_EFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F4 .GT. 0 .AND. IEFS14(F4,V) .GE. ADDPS ) THEN
                    EREF = IEFS14( F4,V ) - ADDPS
                    CALL SETSOURCE_EFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F3 .GT. 0 .AND. IEFS13(F3,V) .GE. ADDPS ) THEN
                    EREF = IEFS13( F3,V ) - ADDPS
                    CALL SETSOURCE_EFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F2 .GT. 0 .AND. IEFS12(F2,V) .GE. ADDPS ) THEN
                    EREF = IEFS12( F2,V ) - ADDPS
                    CALL SETSOURCE_EFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F1 .GT. 0 .AND. IEFS11(F1,V) .GE. ADDPS ) THEN
                    EREF = IEFS11( F1,V ) - ADDPS
                    CALL SETSOURCE_EFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F0 .GT. 0 .AND. IEFS10(F0,V) .GE. ADDPS ) THEN
                    EREF = IEFS10( F0,V ) - ADDPS
                    CALL SETSOURCE_EFS
                    CYCLE                       !  to end of sources-loop

                END IF

C.................  Try for activity-specific FIPS code & SCC match; then
C                           activity-specific FIPS code & left SCC match; then
C                           activity-specific Cy/st code & SCC match; then
C                           activity-specific Cy/st code & left SCC match; then
C                           activity-specific SCC match; then
C                           activity-specific left SCC match

                F5 = FINDC( CFIPSCC, TXCNT( 9 ), CHRT09 ) 
                F4 = FINDC( CFIPSL , TXCNT( 8 ), CHRT08 ) 
                F3 = FINDC( CSTASCC, TXCNT( 6 ), CHRT06 ) 
                F2 = FINDC( CSTASL , TXCNT( 5 ), CHRT05 ) 
                F1 = FINDC( TSCC   , TXCNT( 3 ), CHRT03 ) 
                F0 = FINDC( TSCCL  , TXCNT( 2 ), CHRT02 ) 

                IF( F5 .GT. 0 .AND. IEFS09(F5,V) .GE. ADDPS ) THEN
                    EREF = IEFS09( F5,V ) - ADDPS
                    CALL SETSOURCE_EFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F4 .GT. 0 .AND. IEFS08(F4,V) .GE. ADDPS ) THEN
                    EREF = IEFS08( F4,V ) - ADDPS
                    CALL SETSOURCE_EFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F3 .GT. 0 .AND. IEFS06(F3,V) .GE. ADDPS ) THEN
                    EREF = IEFS06( F3,V ) - ADDPS
                    CALL SETSOURCE_EFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F2 .GT. 0 .AND. IEFS05(F2,V) .GE. ADDPS ) THEN
                    EREF = IEFS05( F2,V ) - ADDPS
                    CALL SETSOURCE_EFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F1 .GT. 0 .AND. IEFS03(F1,V) .GE. ADDPS ) THEN
                    EREF = IEFS03( F1,V ) - ADDPS
                    CALL SETSOURCE_EFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F0 .GT. 0 .AND. IEFS02(F0,V) .GE. ADDPS ) THEN
                    EREF = IEFS02( F0,V ) - ADDPS
                    CALL SETSOURCE_EFS
                    CYCLE                       !  to end of sources-loop

                END IF

C.................  Try for activity-specific FIPS code match
                F0 = FINDC( CFIP, TXCNT( 7 ), CHRT07 ) 

                IF( F0 .GT. 0 .AND. IEFS10(F0,V) .GE. ADDPS ) THEN
                    EREF = IEFS07( F0,V )
                    CALL SETSOURCE_EFS
                    CYCLE                       !  to end of sources-loop
                END IF

C.................  Try for activity-specific country/state code match
                F0 = FINDC( CSTA, TXCNT( 4 ), CHRT04 ) 

                IF( F0 .GT. 0 .AND. IEFS04(F0,V) .GE. ADDPS ) THEN
                    EREF = IEFS04( F0,V )
                    CALL SETSOURCE_EFS
                    CYCLE                       !  to end of sources-loop
                END IF

                IF( IEFS01( V ) .NE. IMISS3 .AND. REPDEFLT ) THEN
                    EREF = IEFS01( V )
                    
                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

                    WRITE( MESG,94010 )
     &                     'NOTE: Using default parameter scheme ' //
     &                     'indices for:'// CRLF() // BLANK5 // 
     &                     BUFFER( 1:L2 ) // CRLF() // BLANK10 // 
     &                     ' SCC: ' // TSCC // ' ACT: ' // ACTVTY( V )
                    CALL M3MESG( MESG )

                    CALL SETSOURCE_EFS

                ELSEIF( IEFS01( V ) .NE. IMISS3 ) THEN
                    EREF = IEFS01( V )
                    CALL SETSOURCE_EFS

                ELSE
                    EFLAG = .TRUE.

                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

                    WRITE( MESG,94010 )
     &                     'ERROR: No parameter scheme indices ' //
     &                     'available (and no default) for:' //
     &                     CRLF() // BLANK5 // BUFFER( 1:L2 ) //
     &                     CRLF() // BLANK10 // 
     &                     ' SCC: ' // TSCC // ' ACT: ' // ACTVTY( V )

                    CALL M3MESG( MESG )

                END IF    !  if default profile code is available or not

            END DO        !  end loop on source, S

        END DO            !  end loop on activity, V

        IF( EFLAG ) THEN
            MESG = 'Problem assigning emission factors to sources'
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

C.............  This internal subprogram stores the index to the parameter
C               scheme index unsorted xref table.
C.............  This is a subroutine so that if more complex steps need
C               to be added later, they can be; and to be consistent with the
C               other ASGN* routines.
C.............  All variables are defined through host association.
            SUBROUTINE SETSOURCE_EFS

C----------------------------------------------------------------------

            EFSIDX( S,J ) = EREF

            RETURN

C------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE SETSOURCE_EFS

        END SUBROUTINE ASGNPSI
