
        SUBROUTINE ASGNTPRO( NGSZ, ANAM, TREFFMT )

C***********************************************************************
C  subroutine body starts at line
C
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
C     Expects tables to have IMISS3 where they are undefined
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 1/99 by M. Houyoux
C
C****************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         ENVINT
        LOGICAL         ENVYN
        INTEGER         FIND1
        INTEGER         FINDC
        INTEGER         INDEX1

        EXTERNAL CRLF, ENVINT, ENVYN, FIND1, FINDC, INDEX1

C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN):: NGSZ          ! no. pols/emis-types in group
        CHARACTER(*), INTENT (IN):: ANAM( NGSZ )  ! group pol names
        CHARACTER(*), INTENT (IN):: TREFFMT       ! temporal x-ref format

C.........  Arrays for evaluating x-ref cases
        LOGICAL          STAT  ( 14 )  ! matching and pollutant status

C.........  Other local variables
        INTEGER          I, J, L2, S, V    !  counters and indices

        INTEGER          ERRCNT  !  count of errors
        INTEGER          F0, F1, F2, F3, F4, F5, F6  ! tmp find indices
        INTEGER       :: F0B = 0 ! extra find index for mobile
        INTEGER       :: F2B = 0 ! extra find index for mobile
        INTEGER       :: F4B = 0 ! extra find index for mobile
        INTEGER          MREF    !  tmp monthly profile code
        INTEGER          WREF    !  tmp weekly  profile code
        INTEGER          DREF    !  tmp diurnal profile code
        INTEGER          MXERR   !  max error messages to output
        INTEGER          MXWARN  !  max warning messages to output
        INTEGER          WRNCNT  !  count of warnings

        LOGICAL       :: EFLAG    = .FALSE. ! true: error found
        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time routine called
        LOGICAL          MFLAG              ! true: use monthly profiles
        LOGICAL          WFLAG              ! true: use weekly  profiles
        LOGICAL, SAVE :: REPDEFLT = .TRUE.  ! true: report default x-ref applied

        CHARACTER*10             RWTFMT   ! fmt to write roadway type to string
        CHARACTER*10             VIDFMT   ! format to write veh ID to string
        CHARACTER*300            BUFFER   ! source fields buffer
        CHARACTER*300            MESG     ! message buffer
        CHARACTER(LEN=SRCLEN3)   CSRC     ! tmp source chars string
        CHARACTER(LEN=FIPLEN3)   CFIP     ! tmp (character) FIPS code
        CHARACTER(LEN=STALEN3)   CSTA     ! tmp Country/state code
        CHARACTER(LEN=SCCLEN3)   TSCC     ! tmp 10-digit SCC
        CHARACTER(LEN=SCCLEN3)   TSCCL    ! tmp left digits of TSCC
        CHARACTER(LEN=SCCLEN3)   TSCCSAV  ! TSCC saved for msg (mb: resets TSCC)
        CHARACTER(LEN=SCCLEN3)   CHKRWT   ! tmp roadway type only SCC
        CHARACTER(LEN=SCCLEN3)   CHKVID   ! tmp vehicle-type only SCC
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
        CHARACTER(LEN=RWTLEN3)   CRWT     ! tmp char roadway type
        CHARACTER(LEN=RWTLEN3)   RWTZERO  ! zero roadway type
        CHARACTER(LEN=VIDLEN3)   CVID     ! tmp vehicle type
        CHARACTER(LEN=VIDLEN3)   VIDZERO  ! zero vehicle type

        CHARACTER*16 :: PROGNAME = 'ASGNTPRO' ! program name

C***********************************************************************
C   begin body of subroutine ASGNTPRO

C.........  For list-formatted temporal cross-reference (one entry per source)
C           from EMS-95 files, the profiles are not applied per pollutant.  So,
C           we can set these for the first group of pollutants used when
C           calling this subroutine and then use them for all pollutants.
        IF( FIRSTIME .AND. TREFFMT .EQ. 'SOURCE' ) THEN

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
                CALL SETSOURCE_TPROFS  ! Sets MDEX, WDEX, DDEX, EDEX

            ENDDO

C.............  Set for remaining pollutants in group  
            DO J = 2, NGSZ
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

C.............  Get error and warning limits from the environment
            MXERR  = ENVINT( ERRSET , ' ', 100, I )
            MXWARN = ENVINT( WARNSET, ' ', 100, I )

            FIRSTIME = .FALSE.

        ENDIF

C.........  Set up roadway type format
        WRITE( RWTFMT, '("(I",I2.2,".",I2.2,")")' ) RWTLEN3, RWTLEN3
        WRITE( VIDFMT, '("(I",I2.2,".",I2.2,")")' ) VIDLEN3, VIDLEN3

C.........  Set up roadway type and vehicle types with all zeros
        RWTZERO = REPEAT( '0', RWTLEN3 )
        VIDZERO = REPEAT( '0', VIDLEN3 )

C.........  Exit subroutine for list-formatted temporal x-ref because we
C           do not have a heirarchial application of temporal profiles
C           to worry about.
        IF( TREFFMT .EQ. 'SOURCE' ) RETURN

        ERRCNT = 0
        WRNCNT = 0
        DO J = 1, NGSZ

C.............  Find index in complete list of pollutants
            V = INDEX1( ANAM( J ), NIPPA, EANAM )

            DO S = 1, NSRC

C.................  Set MFLAG to true for using monthly temporal adjustments
                MFLAG = ( MOD( TPFLAG( S ), MTPRFAC ) .EQ. 0 )

C.................  Set WFLAG to trur for using weekly temporal adjustments
                WFLAG = ( MOD( TPFLAG( S ), WTPRFAC ) .EQ. 0 .OR.
     &                    MOD( TPFLAG( S ), WDTPFAC ) .EQ. 0      )

C.................  Retrieve local variables for source characteristics
                CSRC    = CSOURC( S )
                TSCC    = CSCC( S )
                TSCCL   = TSCC( 1:LSCCEND )
                CFIP    = CSRC( 1:FIPLEN3 )
                CSTA    = CFIP( 1:STALEN3 )
                TSCCSAV = TSCC
                CHK09   = CFIP // TSCC                       ! County // SCC
                CHK08   = CFIP // TSCCL                 ! County // left SCC
                CHK06   = CSTA // TSCC                ! Country/state // SCC
                CHK05   = CSTA // TSCCL          ! Country/state // left SCC

C.................  Set category-specific source characteristic combinations
                SELECT CASE ( CATEGORY )

                CASE ( 'AREA' )   ! Already set above

                CASE ( 'MOBILE' )
         	    WRITE( CRWT, RWTFMT ) IRCLAS( S )
                    WRITE( CVID, VIDFMT ) IVTYPE( S )

                    TSCC = CRWT // CVID
                    CALL PADZERO( TSCC )
                    TSCCL= TSCC( 1:LSCCEND )

                    CHKVID = RWTZERO // CVID
                    CALL PADZERO( CHKVID )

                    CHKRWT = CRWT // VIDZERO
                    CALL PADZERO( CHKRWT )

                    CHK13  = CSRC( 1:MBENDL3(4) )// TSCC   ! Cnty//RWT//LNK//VTP
                    CHK12  = CSRC( 1:MBENDL3(3) )// CHKRWT    ! Cnty// RWT// LNK
                    CHK09  = CFIP // TSCC                   ! County// RWT// VTP
                    CHK08  = CFIP // TSCCL                        ! County// RWT
                    CHK08B = CFIP // CHKVID                       ! County// VTP
                    CHK06  = CSTA // TSCC                   ! State // RWT// VTP
                    CHK05  = CSTA // TSCCL                  ! State // road type
                    CHK05B = CSTA // CHKVID                  ! State // veh type
                    CHK02B = CHKVID                               ! Vehicle type

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

C.................  In the tables used in the following heirarchy, a pollutant-
C                   specific cross-reference entry has not been use as the
C                   default for all pollutants.  So the diurnal profile number
C                   tables (DPRT*) are checked to ensure the pollutant has
C                   been defined for a level of matching of interest.  This is
C                   why DPRT* arrays are compared to IMISS3

C.................  Try to find source characteristic combinations for the
C                   first seven types of matches.  These depend on source
C                   category.

                F6 = FINDC( CHK16, TXCNT( 16 ), CHRT16 ) 
                F5 = FINDC( CHK15, TXCNT( 15 ), CHRT15 ) 
                F4 = FINDC( CHK14, TXCNT( 14 ), CHRT14 ) 
                F3 = FINDC( CHK13, TXCNT( 13 ), CHRT13 ) 
                F2 = FINDC( CHK12, TXCNT( 12 ), CHRT12 ) 
                F1 = FINDC( CHK11, TXCNT( 11 ), CHRT11 ) 
                F0 = FINDC( CHK10, TXCNT( 10 ), CHRT10 )

C.................  Initialize status for all comparisons in first group
                STAT = .FALSE.    ! array

C.................  Evaluate x-ref cases for pollutant/emistype-specific
                STAT(1)= (F6 .GT. 0 .AND. DPRT16(F6,V) .GE. ADDPS)
                STAT(2)= (F5 .GT. 0 .AND. DPRT15(F5,V) .GE. ADDPS)
                STAT(3)= (F4 .GT. 0 .AND. DPRT14(F4,V) .GE. ADDPS)
                STAT(4)= (F3 .GT. 0 .AND. DPRT13(F3,V) .GE. ADDPS)
                STAT(5)= (F2 .GT. 0 .AND. DPRT12(F2,V) .GE. ADDPS)
                STAT(6)= (F1 .GT. 0 .AND. DPRT11(F1,V) .GE. ADDPS)
                STAT(7)= (F0 .GT. 0 .AND. DPRT10(F0,V) .GE. ADDPS)

C.................  Based on evaluation of cases, store reference information
C                   for pollutant-specific
                IF( STAT( 1 ) ) THEN
                    MREF = MPRT16( F6,V )
                    WREF = WPRT16( F6,V )
                    DREF = DPRT16( F6,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 2 ) ) THEN
                    MREF = MPRT15( F5,V )
                    WREF = WPRT15( F5,V )
                    DREF = DPRT15( F5,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 3 ) ) THEN
                    MREF = MPRT14( F4,V )
                    WREF = WPRT14( F4,V )
                    DREF = DPRT14( F4,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 4 ) ) THEN
                    MREF = MPRT13( F3,V )
                    WREF = WPRT13( F3,V )
                    DREF = DPRT13( F3,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 5 ) ) THEN
                    MREF = MPRT12( F2,V )
                    WREF = WPRT12( F2,V )
                    DREF = DPRT12( F2,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 6 ) ) THEN
                    MREF = MPRT11( F1,V )
                    WREF = WPRT11( F1,V )
                    DREF = DPRT11( F1,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 7 ) ) THEN
                    MREF = MPRT10( F0,V )
                    WREF = WPRT10( F0,V )
                    DREF = DPRT10( F0,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                END IF

C.................  Look at the same x-ref cases for no pollutant/emistype
                STAT(1)= (F6 .GT. 0 .AND. DPRT16(F6,V) .NE. IMISS3)
                STAT(2)= (F5 .GT. 0 .AND. DPRT15(F5,V) .NE. IMISS3)
                STAT(3)= (F4 .GT. 0 .AND. DPRT14(F4,V) .NE. IMISS3)
                STAT(4)= (F3 .GT. 0 .AND. DPRT13(F3,V) .NE. IMISS3)
                STAT(5)= (F2 .GT. 0 .AND. DPRT12(F2,V) .NE. IMISS3)
                STAT(6)= (F1 .GT. 0 .AND. DPRT11(F1,V) .NE. IMISS3)
                STAT(7)= (F0 .GT. 0 .AND. DPRT10(F0,V) .NE. IMISS3)

C.................  Continue to evaluate cases and store reference information
                IF( STAT( 1 ) ) THEN
                    MREF = MPRT16( F6,V )
                    WREF = WPRT16( F6,V )
                    DREF = DPRT16( F6,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 2 ) ) THEN
                    MREF = MPRT15( F5,V )
                    WREF = WPRT15( F5,V )
                    DREF = DPRT15( F5,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 3 ) ) THEN
                    MREF = MPRT14( F4,V )
                    WREF = WPRT14( F4,V )
                    DREF = DPRT14( F4,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 4 ) ) THEN
                    MREF = MPRT13( F3,V )
                    WREF = WPRT13( F3,V )
                    DREF = DPRT13( F3,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 5 ) ) THEN
                    MREF = MPRT12( F2,V )
                    WREF = WPRT12( F2,V )
                    DREF = DPRT12( F2,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 6 ) ) THEN
                    MREF = MPRT11( F1,V )
                    WREF = WPRT11( F1,V )
                    DREF = DPRT11( F1,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 7 ) ) THEN
                    MREF = MPRT10( F0,V )
                    WREF = WPRT10( F0,V )
                    DREF = DPRT10( F0,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                END IF

C.................  Try to find source characteristic combinations for the
C                   next six types of matches.
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

C.................  Initialize status for all comparisons in second group
                STAT = .FALSE.    ! array

C.................  Make second round of comparisons for pollutant/emistype
C                   specific cases
                STAT(1)= (F5 .GT. 0 .AND. DPRT09(F5,V) .GE. ADDPS)
                STAT(2)= (F4 .GT. 0 .AND. DPRT08(F4,V) .GE. ADDPS)
                STAT(4)= (F3 .GT. 0 .AND. DPRT06(F3,V) .GE. ADDPS)
                STAT(5)= (F2 .GT. 0 .AND. DPRT05(F2,V) .GE. ADDPS)
                STAT(7)= (F1 .GT. 0 .AND. DPRT03(F1,V) .GE. ADDPS)
                STAT(8)= (F0 .GT. 0 .AND. DPRT02(F0,V) .GE. ADDPS)

C.................  Evaluate mobile-specific cases
                IF( CATEGORY .EQ. 'MOBILE' ) THEN
                    STAT(3)= (F4B .GT. 0 .AND. DPRT09(F4B,V) .GE. ADDPS)
                    STAT(6)= (F2B .GT. 0 .AND. DPRT06(F2B,V) .GE. ADDPS)
                    STAT(9)= (F0B .GT. 0 .AND. DPRT03(F0B,V) .GE. ADDPS)
                END IF

C.................  Continue to evaluate cases and store reference information
                IF( STAT( 1 ) ) THEN
                    MREF = MPRT09( F5,V )
                    WREF = WPRT09( F5,V )
                    DREF = DPRT09( F5,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 2 ) ) THEN
                    MREF = MPRT08( F4,V )
                    WREF = WPRT08( F4,V )
                    DREF = DPRT08( F4,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 3 ) ) THEN
                    MREF = MPRT09( F4B,V )
                    WREF = WPRT09( F4B,V )
                    DREF = DPRT09( F4B,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 4 ) ) THEN
                    MREF = MPRT06( F3,V )
                    WREF = WPRT06( F3,V )
                    DREF = DPRT06( F3,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 5 ) ) THEN
                    MREF = MPRT05( F2,V )
                    WREF = WPRT05( F2,V )
                    DREF = DPRT05( F2,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 6 ) ) THEN
                    MREF = MPRT06( F2B,V )
                    WREF = WPRT06( F2B,V )
                    DREF = DPRT06( F2B,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 7 ) ) THEN
                    MREF = MPRT03( F1,V )
                    WREF = WPRT03( F1,V )
                    DREF = DPRT03( F1,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 8 ) ) THEN
                    MREF = MPRT02( F0,V )
                    WREF = WPRT02( F0,V )
                    DREF = DPRT02( F0,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 9 ) ) THEN
                    MREF = MPRT03( F0B,V )
                    WREF = WPRT03( F0B,V )
                    DREF = DPRT03( F0B,V ) - ADDPS
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                END IF

C.................  Evaluate remainder of x-ref cases
                STAT(1)= (F5 .GT. 0 .AND. DPRT09(F5,V) .NE. IMISS3)
                STAT(2)= (F4 .GT. 0 .AND. DPRT08(F4,V) .NE. IMISS3)
                STAT(4)= (F3 .GT. 0 .AND. DPRT06(F3,V) .NE. IMISS3)
                STAT(5)= (F2 .GT. 0 .AND. DPRT05(F2,V) .NE. IMISS3)
                STAT(7)= (F1 .GT. 0 .AND. DPRT03(F1,V) .NE. IMISS3)
                STAT(8)= (F0 .GT. 0 .AND. DPRT02(F0,V) .NE. IMISS3)

C.................  Remainder of mobile-specific evaluations
                IF( CATEGORY .EQ. 'MOBILE' ) THEN
                    STAT(3)=(F4B .GT. 0 .AND. DPRT09(F4B,V) .NE. IMISS3)
                    STAT(6)=(F2B .GT. 0 .AND. DPRT06(F2B,V) .NE. IMISS3)
                    STAT(9)=(F0B .GT. 0 .AND. DPRT03(F0B,V) .NE. IMISS3)
                END IF

C.................  Continue to evaluate cases and store reference information
C                   for non-pollutant/emission-type specific
C.................  No "ADDPS" used here, because it is 0 in all cases
                IF( STAT( 1 ) ) THEN
                    MREF = MPRT09( F5,V ) 
                    WREF = WPRT09( F5,V )
                    DREF = DPRT09( F5,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 2 ) ) THEN
                    MREF = MPRT08( F4,V ) 
                    WREF = WPRT08( F4,V )
                    DREF = DPRT08( F4,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 3 ) ) THEN
                    MREF = MPRT09( F4B,V )
                    WREF = WPRT09( F4B,V )
                    DREF = DPRT09( F4B,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 4 ) ) THEN
                    MREF = MPRT06( F3,V ) 
                    WREF = WPRT06( F3,V )
                    DREF = DPRT06( F3,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 5 ) ) THEN
                    MREF = MPRT05( F2,V ) 
                    WREF = WPRT05( F2,V )
                    DREF = DPRT05( F2,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 6 ) ) THEN
                    MREF = MPRT06( F2B,V )
                    WREF = WPRT06( F2B,V )
                    DREF = DPRT06( F2B,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 7 ) ) THEN
                    MREF = MPRT03( F1,V ) 
                    WREF = WPRT03( F1,V )
                    DREF = DPRT03( F1,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 8 ) ) THEN
                    MREF = MPRT02( F0,V ) 
                    WREF = WPRT02( F0,V )
                    DREF = DPRT02( F0,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                ELSEIF( STAT( 9 ) ) THEN
                    MREF = MPRT03( F0B,V )
                    WREF = WPRT03( F0B,V )
                    DREF = DPRT03( F0B,V )
                    CALL SETSOURCE_TPROFS
                    CYCLE                       !  to end of sources-loop

                END IF

C.................  Try for county or state match
C.................  NOTE - there is no longer a reason to have this as an
C                   internal subprogram, but no need to change it back either
                CALL COUNTY_OR_STATE
                IF( F0 .GT. 0 ) CYCLE

C.................  Check for and apply ultimate defaults
                IF( MPRT01( V ) .NE. IMISS3 .AND. REPDEFLT .AND.
     &              WRNCNT      .LE. MXWARN                      ) THEN

                    WRNCNT = WRNCNT + 1
                    MREF = MPRT01( V )
                    WREF = WPRT01( V )
                    DREF = DPRT01( V )
                    
                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

                    WRITE( MESG,94010 )
     &                     'NOTE: Using default temporal profile for:'//
     &                     CRLF() // BLANK5 // BUFFER( 1:L2 ) //
     &                     CRLF() // BLANK10 // 
     &                     ' SCC: ' // TSCCSAV // ' POL: ' // ANAM( V )
                    CALL M3MESG( MESG )

                    CALL SETSOURCE_TPROFS

                ELSEIF( MPRT01( V ) .NE. IMISS3 ) THEN
                    MREF = MPRT01( V )
                    WREF = WPRT01( V )
                    DREF = DPRT01( V )
                    CALL SETSOURCE_TPROFS

                ELSE IF( ERRCNT .LE. MXERR ) THEN
                    EFLAG = .TRUE.
                    ERRCNT = ERRCNT + 1

                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

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

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94300   FORMAT( A, I2.2, A, I2.2, A )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram checks charts 04 and 07 for county
C               and state matches 
            SUBROUTINE COUNTY_OR_STATE

C----------------------------------------------------------------------

C.................  Try for any FIPS code match
                F0 = 0
                F0 = FINDC( CFIP, TXCNT( 7 ), CHRT07 ) 

                IF( F0 .GT. 0 ) THEN
                    MREF = MPRT07( F0 ) 
                    WREF = WPRT07( F0 )
                    DREF = DPRT07( F0 )
                    CALL SETSOURCE_TPROFS
                    RETURN                       !  to end of sources-loop
                END IF

C.................  Try for any country/state code match (not, pol-specific)
                F0 = FINDC( CSTA, TXCNT( 4 ), CHRT04 ) 

                IF( F0 .GT. 0 ) THEN
                    MREF = MPRT04( F0 ) 
                    WREF = WPRT04( F0 )
                    DREF = DPRT04( F0 )
                    CALL SETSOURCE_TPROFS
                    RETURN                       !  to end of sources-loop
                END IF

            END SUBROUTINE COUNTY_OR_STATE

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram stores the index of the temporal 
C               profile codes in the temporal profile tables for each source.
C.............  All variables are defined through host association.
            SUBROUTINE SETSOURCE_TPROFS

C----------------------------------------------------------------------

            IF( MFLAG ) THEN

                MDEX( S,J ) = MAX( FIND1( MREF, NMON, MONREF ), 0 )

                IF( MDEX( S,J ) .EQ. 0 ) THEN

                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

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

                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                     'ERROR: Weekly profile', WREF, 
     &                     'is not in profiles, but was assigned' //
     &                     CRLF() // BLANK5 // 'to source:' //
     &                     CRLF() // BLANK5 // BUFFER( 1:L2 )
                    CALL M3MESG( MESG )
                END IF

            END IF

            DDEX( S,J ) = MAX( FIND1( DREF, NWKD, WKDREF ), 0 )

            IF( NEND .GT. 0 ) THEN

                EDEX( S,J ) = MAX( FIND1( DREF, NEND, ENDREF ), 0 )

            END IF  ! If there are weekend diurnal profiles

            IF( DDEX( S,J ) .EQ. 0 .AND. EDEX( S,J ) .EQ. 0 ) THEN

                CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &                 'ERROR: Diurnal profile', DREF, 
     &                 'is not in profiles, but was assigned' //
     &                 CRLF() // BLANK5 // 'to source:' //
     &                 CRLF() // BLANK5 // BUFFER( 1:L2 )
                CALL M3MESG( MESG )

            END IF

            RETURN

C------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE SETSOURCE_TPROFS

        END SUBROUTINE ASGNTPRO
