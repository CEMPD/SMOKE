
        SUBROUTINE ASGNUNCERT
        
C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      For each source, find the most specific entry in the uncertainty
C      table to apply.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines:  I/O API subroutines, CHECKMEM
C      Functions:    CRLF, FIND1, FINDC, INDEX1
C
C  REVISION  HISTORY:
C     Created 9/2001 by A. Holland
C
C************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
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
C*************************************************************************

C...........   MODULES for public variables   
C...........   This module contains the source ararys
        USE MODSOURC

C.........  This module is for cross reference tables
        USE MODXREF

C.........  This module contains the information about the source category
        USE MODINFO
        
C.........  This module contains uncertainty-specific settings        
        USE MODUNCERT

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         FIND1
        INTEGER         FINDC
        INTEGER         INDEX1

        EXTERNAL CRLF, FIND1, FINDC, INDEX1


C.........  Local variables
        INTEGER          I, J, K, L, L2              !  counters and indices
        INTEGER          S, V, IOS, N                !  counters and indices
        INTEGER          F0, F1, F2, F3, F4, F5, F6  ! tmp find indices
        INTEGER          IDX             !  tmp index to uncert data table
        INTEGER       :: UINDX( NSRC, UNVAR )        ! index to uncert table
        
        REAL             RIMISS3            ! real value of integer missing

        LOGICAL       :: EFLAG    = .FALSE. ! true: error detected
        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time subrtn called

        CHARACTER*8            FMTFIP   ! format for writing FIPS code
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
        CHARACTER(LEN=FPLLEN3) CFIPPLT  ! tmp FIPS code // plant id
        CHARACTER(LEN=FPSLEN3) CFIPSCC  ! tmp FIPS code // SCC
        CHARACTER(LEN=FPSLEN3) CFIPSL   ! tmp FIPS code // left SCC
        CHARACTER(LEN=SCCLEN3) TSCC     ! tmp 10-digit SCC
        
        CHARACTER*6         :: PKTTYP = 'UNCERT'
        CHARACTER*1            CHARMETH ! local method
        CHARACTER*1            CTYPE    ! local parametric/empirical type
        CHARACTER*2            CAPRCH   ! local approach
        CHARACTER(LEN=IOVLEN3) CBUF     ! tmp pollutant name

        CHARACTER*16 :: PROGNAME = 'ASGNUNCERT' ! program name

C***********************************************************************
C   begin body of subroutine ASGNUNCERT

C.........  Compute real value of integer missing
        RIMISS3 = REAL( IMISS3 )

C.........  For first time routine is called ...
        IF( FIRSTIME ) THEN

C.............  Set up format for creating character FIPS code for non-point
            IF( CATEGORY .NE. 'POINT' ) THEN
                WRITE( FMTFIP, 94300 ) '(I', FIPLEN3, '.', FIPLEN3, ')'
            ENDIF

            FIRSTIME = .FALSE.

        ENDIF

C.........  Initialize UINDX and USTAT
	ALLOCATE( USTAT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'USTAT', PROGNAME )
        
        UINDX = 0        ! Array
        USTAT = .FALSE.  ! Array

C.........  For each pollutant of interest
        DO J = 1, UNVAR
        
C.............  Find index in complete list of pollutants
            CBUF = UNAMES( J )
            V = INDEX1( CBUF, UNVAR, UNAMES )

            

            L = LEN_TRIM( CBUF )
            MESG = BLANK5 // 'Assigning uncertainty to ' //
     &             'pollutant "' // CBUF( 1:L ) // '"...'
            CALL M3MSG2( MESG )                
                                                

            DO S = 1, NSRC
	    
                CSRC    = CSOURC( S )
                CFIP    = CSRC( 1:FIPLEN3 )
                CSTA    = CFIP( 1:STALEN3 )
                TSCC    = CSCC( S )
                TSCCL   = TSCC( 1:LSCCEND )
                CFIPSCC = CFIP // TSCC
                CFIPSL  = CFIP // TSCCL
                CSTASCC = CSTA // TSCC
                CSTASL  = CSTA // TSCCL

C.................  Create selection 
                SELECT CASE ( CATEGORY )

                CASE ( 'AREA' ) 
                    CRL     = 'A'
                    CATDESC = 'Area'

                CASE ( 'MOBILE' )
                    CRL     = 'M'
                    CATDESC = 'Mobile'

                CASE ( 'POINT' )
                    CRL     = 'P'
                    CATDESC = 'Point'

                    CSSC5   = CSRC( 1:PTENDL3( 7 ) ) // TSCC
                    CSSC4   = CSRC( 1:PTENDL3( 6 ) ) // TSCC
                    CSSC3   = CSRC( 1:PTENDL3( 5 ) ) // TSCC
                    CSSC2   = CSRC( 1:PTENDL3( 4 ) ) // TSCC
                    CSSC1   = CSRC( 1:PTENDL3( 3 ) ) // TSCC
                    CSSC0   = CSRC( 1:PTENDL3( 2 ) ) // TSCC
                    CFIPPLT = CSRC( 1:PTENDL3( 2 ) )
                    
                CASE DEFAULT

                END SELECT

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

            IF( F6 .GT. 0 .AND. ICTL16(F6,V) .GE. ADDPS ) THEN
                IDX = ICTL16( F6,V ) - ADDPS
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F5 .GT. 0 .AND. ICTL15(F5,V) .GE. ADDPS ) THEN
                IDX = ICTL15( F5,V ) - ADDPS
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F4 .GT. 0 .AND. ICTL14(F4,V) .GE. ADDPS ) THEN
                IDX = ICTL14( F4,V ) - ADDPS
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F3 .GT. 0 .AND. ICTL13(F3,V) .GE. ADDPS ) THEN
                IDX = ICTL13( F3,V ) - ADDPS
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F2 .GT. 0 .AND. ICTL12(F2,V) .GE. ADDPS ) THEN
                IDX = ICTL12( F2,V ) - ADDPS
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F1 .GT. 0 .AND. ICTL11(F1,V) .GE. ADDPS ) THEN
                IDX = ICTL11( F1,V ) - ADDPS
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F0 .GT. 0 .AND. ICTL10(F0,V) .GE. ADDPS ) THEN
                IDX = ICTL10( F0,V ) - ADDPS
                CALL SETSOURCE_UNCERT_INDEX
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
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F5 .GT. 0 .AND. ICTL15(F5,V) .NE. IMISS3 ) THEN
                IDX = ICTL15( F5,V )
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F4 .GT. 0 .AND. ICTL14(F4,V) .NE. IMISS3 ) THEN
                IDX = ICTL14( F4,V )
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F3 .GT. 0 .AND. ICTL13(F3,V) .NE. IMISS3 ) THEN
                IDX = ICTL13( F3,V )
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F2 .GT. 0 .AND. ICTL12(F2,V) .NE. IMISS3 ) THEN
                IDX = ICTL12( F2,V )
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F1 .GT. 0 .AND. ICTL11(F1,V) .NE. IMISS3 ) THEN
                IDX = ICTL11( F1,V )
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F0 .GT. 0 .AND. ICTL10(F0,V) .NE. IMISS3 ) THEN
                IDX = ICTL10( F0,V )
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            END IF

C.................  Try for pollutant-specific FIPS code & SCC match; then
C                           pollutant-specific FIPS code & left SCC match; then
C                           pollutant-specific Cy/st code & SCC match; then
C                           pollutant-specific Cy/st code & left SCC match; then
C                           pollutant-specific SCC match; then
C                           pollutant-specific left SCC match

            F5 = FINDC( CFIPSCC, TXCNT( 9 ), CHRT09 ) 
            F4 = FINDC( CFIPSL , TXCNT( 8 ), CHRT08 ) 
            F3 = FINDC( CSTASCC, TXCNT( 6 ), CHRT06 ) 
            F2 = FINDC( CSTASL , TXCNT( 5 ), CHRT05 ) 
            F1 = FINDC( TSCC   , TXCNT( 3 ), CHRT03 ) 
            F0 = FINDC( TSCCL  , TXCNT( 2 ), CHRT02 ) 

            IF( F5 .GT. 0 .AND. ICTL09(F5,V) .GE. ADDPS ) THEN
                IDX = ICTL09( F5,V ) - ADDPS
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F4 .GT. 0 .AND. ICTL08(F4,V) .GE. ADDPS ) THEN
                IDX = ICTL08( F4,V ) - ADDPS
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop
                
            ELSEIF( F3 .GT. 0 .AND. ICTL06(F3,V) .GE. ADDPS ) THEN
                IDX = ICTL06( F3,V ) - ADDPS
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F2 .GT. 0 .AND. ICTL05(F2,V) .GE. ADDPS ) THEN
                IDX = ICTL05( F2,V ) - ADDPS
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F1 .GT. 0 .AND. ICTL03(F1,V) .GE. ADDPS ) THEN
                IDX = ICTL03( F1,V ) - ADDPS
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F0 .GT. 0 .AND. ICTL02(F0,V) .GE. ADDPS ) THEN
                IDX = ICTL02( F0,V ) - ADDPS
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            END IF

C.................  Try for any FIPS code & SCC match; then
C                           any FIPS code & left SCC match; then
C                           any Cy/st code & SCC match; then
C                           any Cy/st code & left SCC match; then
C                           any SCC match; then
C                           any left SCC match

            IF( F5 .GT. 0 .AND. ICTL09(F5,V) .NE. IMISS3 ) THEN
                IDX = ICTL09( F5,V )
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F4 .GT. 0 .AND. ICTL08(F4,V) .NE. IMISS3 ) THEN
                IDX = ICTL08( F4,V )
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop
                    
            ELSEIF( F3 .GT. 0 .AND. ICTL06(F3,V) .NE. IMISS3 ) THEN
                IDX = ICTL06( F3,V )
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F2 .GT. 0 .AND. ICTL05(F2,V) .NE. IMISS3 ) THEN
                IDX = ICTL05( F2,V )
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F1 .GT. 0 .AND. ICTL03(F1,V) .NE. IMISS3 ) THEN
                IDX = ICTL03( F1,V )
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F0 .GT. 0 .AND. ICTL02(F0,V) .NE. IMISS3 ) THEN
                IDX = ICTL02( F0,V )
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            END IF

C.................  Try for any FIPS code match
            F0 = FINDC( CFIP, TXCNT( 7 ), CHRT07 ) 

            IF( F0 .GT. 0 .AND. ICTL07(F0,V) .GE. ADDPS ) THEN
                IDX = ICTL07( F0,V ) - ADDPS
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F0 .GT. 0 .AND. ICTL07(F0,V) .NE. IMISS3 ) THEN
                IDX = ICTL07( F0,V )
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            END IF

C.................  Try for any country/state code match 
            F0 = FINDC( CSTA, TXCNT( 4 ), CHRT04 ) 

            IF( F0 .GT. 0 .AND. ICTL04(F0,V) .GE. ADDPS ) THEN
                IDX = ICTL04( F0,V )
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            ELSEIF( F0 .GT. 0 .AND. ICTL04(F0,V) .NE. IMISS3 ) THEN
                IDX = ICTL04( F0,V )
                CALL SETSOURCE_UNCERT_INDEX
                CYCLE                       !  to end of sources-loop

            END IF

c            IF( ICTL01( V ) .NE. IMISS3 .AND. REPDEFLT ) THEN
c                IDX = ICTL01( V )
                    
c                CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

c                WRITE( MESG,94010 )
c     &                 'NOTE: Using default ' // PKTTYP // 
c     &                 ' control packet entry for:'//
c     &                 CRLF() // BLANK5 // BUFFER( 1:L2 ) //
c     &                 CRLF() // BLANK10 // 
c     &                 ' SCC: ' // TSCC // ' POL: ' // EANAM( V )
c                CALL M3MESG( MESG )

c                CALL SETSOURCE_UNCERT_INDEX

            IF( ICTL01( V ) .NE. IMISS3 ) THEN
                IDX = ICTL01( V )
                CALL SETSOURCE_UNCERT_INDEX

            END IF    !  if default profile code is available or not

            END DO        !  end loop on source, S
            
        END DO            !  end loop on pollutant, J/V


        IF( EFLAG ) THEN
            MESG = 'Problem assigning uncertainty to sources'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF 
        
C.........  Get count of all effected sources
        DO S = 1, NSRC
            IF( USTAT( S ) ) UCOUNT = UCOUNT + 1
        END DO

C.........  Allocate and initialize output variables        
        ALLOCATE( SRCNUM( UCOUNT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCNUM', PROGNAME )
        ALLOCATE( METHOD( UCOUNT, UNVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'METHOD', PROGNAME )
        ALLOCATE( EPTYP( UCOUNT, UNVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EPTYP', PROGNAME )
        ALLOCATE( APRCH( UCOUNT, UNVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'APRCH', PROGNAME )
        ALLOCATE( NUMEP( UCOUNT, UNVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NUMEP', PROGNAME )
        ALLOCATE( UNCIDX( UCOUNT, UNVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UNCIDX', PROGNAME )
        ALLOCATE( PARMS( MXPARDAT, NPPCKT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PARMS', PROGNAME )
        ALLOCATE( EMFVAL( MXEMPDAT, NEPCKT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMFVAL', PROGNAME )
        ALLOCATE( PROBVAL( MXEMPDAT, NEPCKT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PROBVAL', PROGNAME )
        
        SRCNUM  = IMISS3
        METHOD  = IMISS3
        EPTYP   = IMISS3
        APRCH   = IMISS3
        NUMEP   = IMISS3
        UNCIDX  = IMISS3
        PARMS   = RIMISS3
        EMFVAL  = RIMISS3
        PROBVAL = RIMISS3
        
C.........  Get source numbers for effected sources        
        N = 0
        
        DO S = 1, NSRC
            IF( USTAT( S ) ) THEN
                N = N + 1
                SRCNUM( N ) = S
            END IF
        END DO
        
C.........  Fill output arrays for effected sources and pollutants
        N = 0
        
        DO J = 1, UNVAR          
        
            DO S = 1, NSRC
            
                IF( USTAT( S ) ) N = N + 1
            
                IF( UINDX( S, J ) .GT. 0 ) THEN
                
                    IDX = UINDX( S, J )
                
                    CHARMETH = FAPCKT( IDX )%METH
                                
                    IF( CHARMETH .EQ. 'E' ) THEN
                        METHOD( N, J ) = 0
                    ELSE
                        METHOD( N, J ) = 1
                    END IF
                    
                    CAPRCH = FAPCKT( IDX )%APRCH
                    
                    IF( CAPRCH .EQ. 'S'  ) APRCH( N, J ) = 1
                    IF( CAPRCH .EQ. 'I'  ) APRCH( N, J ) = 2
                    IF( CAPRCH .EQ. 'ST' ) APRCH( N, J ) = 3
                    IF( CAPRCH .EQ. 'IT' ) APRCH( N, J ) = 4
                
                    IF( METHOD( N, J ) .EQ. 1 ) THEN
                
                        CTYPE = PARPCKT( FAPCKT( IDX )%INDX )%PTYPE
                    
                        IF( CTYPE .EQ. 'N' ) EPTYP( N, J ) = 1 
                        IF( CTYPE .EQ. 'L' ) EPTYP( N, J ) = 2
                        IF( CTYPE .EQ. 'G' ) EPTYP( N, J ) = 3
                        IF( CTYPE .EQ. 'W' ) EPTYP( N, J ) = 4
                        IF( CTYPE .EQ. 'B' ) EPTYP( N, J ) = 5
                    
                        NUMEP( N, J ) = PARPCKT( FAPCKT
     &                                  ( IDX )%INDX )%NUMP
     
                        UNCIDX( N, J ) = FAPCKT(IDX)%INDX
                    
                    
                    ELSE 
                    
                        CTYPE = EMPPCKT( FAPCKT( IDX )%INDX )%ETYPE
                    
                        IF( CTYPE .EQ. 'S' ) EPTYP( N, J ) = 6
                        IF( CTYPE .EQ. 'L' ) EPTYP( N, J ) = 7
                    
                        NUMEP( N, J ) = EMPPCKT( FAPCKT( IDX )
     &                                  %INDX )%EPKTENT
     
                        UNCIDX( N, J ) = FAPCKT( IDX )%INDX

                                                                    
                    END IF
                    
                END IF                    

            END DO
            N = 0
        END DO
        
        DO I = 1, MXPARDAT
            
            DO J = 1, NPPCKT
            
                PARMS( I, J ) = PARAMET( J, I )
     
            END DO
            
        END DO
             
        DO I = 1, MXEMPDAT
        
            DO J = 1, NEPCKT
            
                EMFVAL( I, J ) = EFVAL( J, I )

                PROBVAL( I, J ) = PROB( J, I )
                
            END DO
            
        END DO


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
            SUBROUTINE SETSOURCE_UNCERT_INDEX

C----------------------------------------------------------------------

            UINDX( S,J ) = IDX

            IF( IDX .GT. 0 ) THEN
                USTAT( S )  = .TRUE.   ! true indicates source is affected
            ELSE
                USTAT( S )  = .FALSE.
            END IF

            END SUBROUTINE SETSOURCE_UNCERT_INDEX

        END SUBROUTINE ASGNUNCERT
