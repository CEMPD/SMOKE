
        SUBROUTINE ASGNPSI( NACT, ANAM, ISTAT )

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
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
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
        INTEGER     , INTENT (IN) :: NACT          ! number of activities 
        CHARACTER(*), INTENT (IN) :: ANAM ( NACT ) ! activity names
        INTEGER     , INTENT (IN) :: ISTAT( NACT ) ! 0=don't use; >0=use

C.........  Other local variables
        INTEGER          I, J, L2, S, V    !  counters and indices

        INTEGER          F0, F1, F2, F3, F4, F5, F6  ! tmp find indices
        INTEGER       :: F0B = 0 ! extra find index for mobile
        INTEGER       :: F2B = 0 ! extra find index for mobile
        INTEGER       :: F4B = 0 ! extra find index for mobile
        INTEGER          EREF           !  tmp diurnal profile code

        LOGICAL       :: EFLAG    = .FALSE.
        LOGICAL, SAVE :: FIRSTIME = .TRUE.
        LOGICAL, SAVE :: REPDEFLT = .TRUE.

        CHARACTER*10             RWTFMT   ! roadway type format
        CHARACTER*300            BUFFER   ! source fields buffer
        CHARACTER*300            MESG     ! message buffer
        CHARACTER(LEN=FIPLEN3):: CFIP     ! tmp FIPS code
        CHARACTER(LEN=STALEN3):: CSTA     ! tmp st/county code
        CHARACTER(LEN=LNKLEN3):: CLNK     ! temporary link code
        CHARACTER(LEN=RWTLEN3):: CRWT     ! buffer for roadway type
        CHARACTER(LEN=SRCLEN3):: CSRC =' '! tmp source chars string
        CHARACTER(LEN=VIDLEN3):: CVID =' '! buffer for vehicle type number
        CHARACTER(LEN=SS5LEN3):: CHK16=' '! tmp source chars through char5// SCC
        CHARACTER(LEN=SS4LEN3):: CHK15=' '! tmp source chars through char4// SCC
        CHARACTER(LEN=SS3LEN3):: CHK14=' '! tmp source chars through char3// SCC
        CHARACTER(LEN=SS2LEN3):: CHK13=' '! tmp source chars through char2// SCC
        CHARACTER(LEN=SS1LEN3):: CHK12=' '! tmp source chars through char1// SCC
        CHARACTER(LEN=SS0LEN3):: CHK11=' '! tmp FIPS // Plant // SCC
        CHARACTER(LEN=FPLLEN3):: CHK10=' '! tmp FIPS code // plant id
        CHARACTER(LEN=FPSLEN3):: CHK09=' '! tmp FIPS code // SCC
        CHARACTER(LEN=FPSLEN3):: CHK08=' '! tmp FIPS code // left SCC
        CHARACTER(LEN=FPSLEN3):: CHK08B=' '! tmp FIPS code // veh ID SCC
        CHARACTER(LEN=STSLEN3):: CHK06=' '! tmp Country/state code // SCC
        CHARACTER(LEN=STSLEN3):: CHK05=' '! tmp Country/state code // left SCC
        CHARACTER(LEN=STSLEN3):: CHK05B=' '! tmp Country/state code// veh ID SCC
        CHARACTER(LEN=SCCLEN3):: CHK02B=' '! tmp veh ID SCC
        CHARACTER(LEN=SCCLEN3)   CHKRWT   ! tmp roadway type only SCC
        CHARACTER(LEN=SCCLEN3)   CHKVID   ! tmp vehicle-type only SCC
        CHARACTER(LEN=RWTLEN3):: RWTZERO=' '! zero roadway type
        CHARACTER(LEN=SCCLEN3):: TSCC =' '! tmp 10-digit SCC
        CHARACTER(LEN=SCCLEN3):: TSCCL=' '! tmp left digits of TSCC
        CHARACTER(LEN=VIDLEN3):: VIDHOLD  ! placeholder vehicle type
        CHARACTER(LEN=VIDLEN3):: VIDZERO  ! zero vehicle type

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
        WRITE( RWTFMT, '("(I",I2.2,")")' ) RWTLEN3

C.........  Set up roadway type and vehicle types with all zeros
        RWTZERO = REPEAT( '0', RWTLEN3 )
        VIDZERO = REPEAT( '0', VIDLEN3 )
        VIDHOLD = REPEAT( '-', VIDLEN3 )

C.........  Initialize search fields to blank
        CSRC    = ' '
        TSCC    = ' '
        CFIP    = ' '
        CSTA    = ' '
        TSCCL   = ' '

C.........  Loop through activities
        DO J = 1, NACT

C.............  Skip activity if status is to not use
            IF( ISTAT( J ) .EQ. 0 ) CYCLE

C.............  Find index in complete list of activities
            V = INDEX1( ANAM( J ), NIACT, ACTVTY )

            DO S = 1, NSRC

                CSRC  = CSOURC( S )
                TSCC  = CSCC( S )
                TSCCL = TSCC( 1:LSCCEND )
                CFIP  = CSRC( 1:FIPLEN3 )
                CSTA  = CFIP( 1:STALEN3 )
                CHK09 = CFIP // TSCC
                CHK08 = CFIP // TSCCL
                CHK06 = CSTA // TSCC
                CHK05 = CSTA // TSCCL

C.................  Create selection 
                SELECT CASE ( CATEGORY )

                CASE ( 'AREA' ) 

C.................  For PSIs, SCC is not used for matching, but roadway type
C                   and link IDs are used, so adjust fields accordingly
                CASE ( 'MOBILE' )

                    CRWT    = CSRC( MBBEGL3( 2 ):MBENDL3( 2 ) )
                    CLNK    = CSRC( MBBEGL3( 3 ):MBENDL3( 3 ) )
                    CVID    = CSRC( MBBEGL3( 4 ):MBENDL3( 4 ) )
                    CALL PADZERO( CVID )

                    TSCC = CRWT // CVID
                    CALL PADZERO( TSCC )
                    TSCCL= TSCC( 1:LSCCEND )

                    CHKVID = RWTZERO // CVID
                    CALL PADZERO( CHKVID )

                    CHKRWT = CRWT // VIDHOLD
                    CALL PADZERO( CHKRWT )

                    CHK09  = CFIP // TSCC
                    CHK08  = CFIP // TSCCL 
                    CHK08B = CFIP // CHKVID                  ! County// VTP
                    CHK06  = CSTA // TSCC
                    CHK05  = CSTA // TSCCL
                    CHK05B = CSTA // CHKVID                   ! State // vtype
                    CHK02B = CHKVID                           ! Vehicle type
                    CHK12  = CSRC( 1:MBENDL3( 3 ) ) // CHKRWT ! County//rcl//lnk
                    CHK13  = CSRC( 1:MBENDL3( 4 ) ) // TSCC   ! w/ vehicle type

                CASE ( 'POINT' )

                    CHK16   = CSRC( 1:PTENDL3( 7 ) ) // TSCC
                    CHK15   = CSRC( 1:PTENDL3( 6 ) ) // TSCC
                    CHK14   = CSRC( 1:PTENDL3( 5 ) ) // TSCC
                    CHK13   = CSRC( 1:PTENDL3( 4 ) ) // TSCC
                    CHK12   = CSRC( 1:PTENDL3( 3 ) ) // TSCC
                    CHK11   = CSRC( 1:PTENDL3( 2 ) ) // TSCC 
                    CHK10   = CSRC( 1:PTENDL3( 2 ) )           ! County // plant
                    
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

                F6 = FINDC( CHK16, TXCNT( 16 ), CHRT16 ) 
                F5 = FINDC( CHK15, TXCNT( 15 ), CHRT15 ) 
                F4 = FINDC( CHK14, TXCNT( 14 ), CHRT14 ) 
                F3 = FINDC( CHK13, TXCNT( 13 ), CHRT13 ) 
                F2 = FINDC( CHK12, TXCNT( 12 ), CHRT12 ) 
                F1 = FINDC( CHK11, TXCNT( 11 ), CHRT11 ) 
                F0 = FINDC( CHK10, TXCNT( 10 ), CHRT10 )

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

                F5 = FINDC( CHK09, TXCNT( 9 ), CHRT09 ) 
                F4 = FINDC( CHK08, TXCNT( 8 ), CHRT08 ) 
                F3 = FINDC( CHK06, TXCNT( 6 ), CHRT06 ) 
                F2 = FINDC( CHK05, TXCNT( 5 ), CHRT05 ) 
                F1 = FINDC( TSCC , TXCNT( 3 ), CHRT03 ) 
                F0 = FINDC( TSCCL, TXCNT( 2 ), CHRT02 )

C................. Check for mobile-specific matches that use a TSCC with
C                  road class of zero and vehicle type. The assignment of
C                  temporal profile based on  a vehicle type and no road class
C                  comes after the road class only match (or TSCCL in CHRT08,
C                  for example) but the match uses the full TSCC (or CHRT09, for
C                  example).
                IF( CATEGORY .EQ. 'MOBILE' ) THEN
                    F4B = FINDC( CHK08B, TXCNT( 9 ), CHRT09 )
                    F2B = FINDC( CHK05B, TXCNT( 6 ), CHRT06 )
                    F0B = FINDC( CHK02B, TXCNT( 3 ), CHRT03 )
                END IF

                IF( F5 .GT. 0 .AND. IEFS09(F5,V) .GE. ADDPS ) THEN
                    EREF = IEFS09( F5,V ) - ADDPS
                    CALL SETSOURCE_EFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F4 .GT. 0 .AND. IEFS08(F4,V) .GE. ADDPS ) THEN
                    EREF = IEFS08( F4,V ) - ADDPS
                    CALL SETSOURCE_EFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( F4B .GT. 0 .AND. IEFS08(F4B,V) .GE. ADDPS ) THEN
                    EREF = IEFS08( F4B,V ) - ADDPS
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

                ELSEIF( F2B .GT. 0 .AND. IEFS05(F2B,V) .GE. ADDPS ) THEN
                    EREF = IEFS05( F2B,V ) - ADDPS
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

                ELSEIF( F0B .GT. 0 .AND. IEFS03(F0B,V) .GE. ADDPS ) THEN
                    EREF = IEFS03( F0B,V ) - ADDPS
                    CALL SETSOURCE_EFS
                    CYCLE                       !  to end of sources-loop

                END IF

C.................  Try for activity-specific FIPS code match
                F0 = FINDC( CFIP, TXCNT( 7 ), CHRT07 ) 

                IF( F0 .GT. 0 .AND. IEFS07(F0,V) .GE. ADDPS ) THEN
                    EREF = IEFS07( F0,V ) - ADDPS
                    CALL SETSOURCE_EFS
                    CYCLE                       !  to end of sources-loop
                END IF

C.................  Try for activity-specific country/state code match
                F0 = FINDC( CSTA, TXCNT( 4 ), CHRT04 ) 

                IF( F0 .GT. 0 .AND. IEFS04(F0,V) .GE. ADDPS ) THEN
                    EREF = IEFS04( F0,V ) - ADDPS
                    CALL SETSOURCE_EFS
                    CYCLE                       !  to end of sources-loop
                END IF

                IF( IEFS01( V ) .NE. IMISS3 .AND. REPDEFLT ) THEN
                    EREF = IEFS01( V ) - ADDPS
                    
                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

                    WRITE( MESG,94010 )
     &                     'NOTE: Using default parameter scheme ' //
     &                     'indices for:'// CRLF() // BLANK5 // 
     &                     BUFFER( 1:L2 ) // CRLF() // BLANK10 // 
     &                     ' SCC: ' // TSCC // ' ACT: ' // ACTVTY( V )
                    CALL M3MESG( MESG )

                    CALL SETSOURCE_EFS

                ELSEIF( IEFS01( V ) .NE. IMISS3 ) THEN
                    EREF = IEFS01( V ) - ADDPS
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
