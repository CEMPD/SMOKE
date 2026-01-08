
        SUBROUTINE ASGNSURG

C***********************************************************************
C  subroutine body starts at line 109
C
C  DESCRIPTION:
C      For each source, find the most specific gridding surrogate code
C      that applies to that source. Do this using the grouped tables of
C      gridding cross references from RDGREF.  The hierarchical order is
C      defined in this subroutine, and can be determined from the in-source
C      comments below. Once a surrogate code has been identified, search for 
C      this code in the gridding surrogates tables (from RDSRG) and save the 
C      index to these tables for each source.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 4/99 by M. Houyoux
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
C***************************************************************************

C...........   MODULES for public variables   
C...........   This module contains the source ararys
        USE M3UTILIO

        USE MODSOURC, ONLY: CIFIP, CSOURC, CSCC, CLINK, IRCLAS, IVTYPE

C...........   This module contains the cross-reference tables
        USE MODXREF, ONLY: CHRT02, CHRT03, CHRT04, CHRT05, 
     &                     CHRT06, CHRT07, CHRT08, CHRT09,
     &                     ISRG01, ISRG02, ISRG03, ISRG04, ISRG05,
     &                     ISRG06, ISRG07, ISRG08, ISRG09,
     &                     TXCNT, ASRGID, SRGFIPIDX

C...........   This module contains the gridding surrogates tables
        USE MODSURG, ONLY: NSRGS, SRGLIST

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, LSCCEND, NCHARS, NSRC

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  i/o api constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
C       CHARACTER(2)    CRLF
C       LOGICAL         ENVYN
C       INTEGER         FIND1
C       INTEGER         FINDC
        LOGICAL         SETSCCTYPE

C        EXTERNAL    CRLF, ENVYN, FIND1, FINDC, SETSCCTYPE
        EXTERNAL     SETSCCTYPE

C.........  Other local variables
        INTEGER          I, II, J, L2, S    !  counters and indices

        INTEGER          F0, F1, F2, F3, F4, F5  ! tmp find indices
        INTEGER          ISRG    !  tmp surrogate code

        LOGICAL       :: EFLAG    = .FALSE.
        LOGICAL, SAVE :: FIRSTIME = .TRUE.
        LOGICAL, SAVE :: REPDEFLT = .TRUE.
        LOGICAL          SCCFLAG           ! true: SCC type is different from previous

        CHARACTER(10)         RWTFMT   ! formt to write rdway type to string
        CHARACTER(10)         VIDFMT   ! format to write veh ID to string
        CHARACTER(300)        BUFFER   ! source fields buffer
        CHARACTER(300)        MESG     ! message buffer
        CHARACTER(LNKLEN3) :: CLNK= ' '! tmp link ID
        CHARACTER(STALEN3)    CSTA     ! tmp Country/state code
        CHARACTER(STSLEN3)    CSTASCC  ! tmp Country/state code // SCC
        CHARACTER(STSLEN3)    CSTASL   ! tmp Country/state code // left SCC
        CHARACTER(SCCLEN3)    TSCCL    ! tmp left digits of TSCC
        CHARACTER(SRCLEN3)    CSRC     ! tmp source chars string
        CHARACTER(RWTLEN3)    CRWT     !  buffer for roadway type
        CHARACTER(FIPLEN3)    CFIP     ! tmp (character) FIPS code
        CHARACTER(FPLLEN3)    CFIPPLT  ! tmp FIPS code // plant id
        CHARACTER(FPSLEN3)    CFIPSCC  ! tmp FIPS code // SCC
        CHARACTER(FPSLEN3)    CFIPSL   ! tmp FIPS code // left SCC
        CHARACTER(SCCLEN3)    TSCC     ! tmp 10-digit SCC
        CHARACTER(VIDLEN3)    CVID     ! buffer for vehicle type ID

        CHARACTER(16) :: PROGNAME = 'ASGNSURG' ! program name

C***********************************************************************
C   begin body of subroutine ASGNSURG

C.........  For first time routine is called in all cases,
        IF( FIRSTIME ) THEN

C.............  Retrieve environment variables
            MESG = 'Switch for reporting default gridding surrogates'
            REPDEFLT = ENVYN ( 'REPORT_DEFAULTS', MESG, .TRUE., I )

C.............  Allocate memory for surrogate ID position

            FIRSTIME = .FALSE.

        ENDIF

C.........  Set up formats
        WRITE( RWTFMT, '("(I",I2.2,".",I2.2,")")' ) RWTLEN3, RWTLEN3
        WRITE( VIDFMT, '("(I",I2.2,".",I2.2,")")' ) VIDLEN3, VIDLEN3
        
C.........  Loop through the sources
        DO S = 1, NSRC

C.............  Create selection 
            SELECT CASE ( CATEGORY )

            CASE ( 'AREA' )
                CSRC    = CSOURC( S )
                CFIP    = CSRC( 1:FIPLEN3 )
                CSTA    = CFIP( 1:STALEN3 )
                TSCC    = CSCC( S )
                
C.................  Set type of SCC                
                SCCFLAG = SETSCCTYPE ( TSCC )
                TSCCL   = TSCC( 1:LSCCEND )
                
                CFIPSCC = CFIP // TSCC
                CFIPSL  = CFIP // TSCCL
                CSTASCC = CSTA // TSCC
                CSTASL  = CSTA // TSCCL

            CASE ( 'BIOG' )

c note: insert here when needed

            CASE ( 'MOBILE' )

                CSRC    = CSOURC( S )
                CFIP    = CSRC( 1:FIPLEN3 )
                CSTA    = CFIP( 1:STALEN3 )
                TSCC    = CSCC( S )

C.................  Set type of SCC                
                SCCFLAG = SETSCCTYPE ( TSCC )
                TSCCL   = TSCC( 1:LSCCEND )

                CFIPSCC = CFIP // TSCC
                CFIPSL  = CFIP // TSCCL
                CSTASCC = CSTA // TSCC
                CSTASL  = CSTA // TSCCL

            END SELECT

C.............  Skip finding a surrogate if current source is a link source
            IF( CLNK .NE. ' ' ) CYCLE

C.................  Try for FIPS code & SCC match; then
C                           FIPS code & left SCC match; then
C                           Cy/st code & SCC match; then
C                           Cy/st code & left SCC match; then
C                           SCC match; then
C                           left SCC match

            F5 = FINDC( CFIPSCC, TXCNT( 9 ), CHRT09 ) 
            F4 = FINDC( CFIPSL , TXCNT( 8 ), CHRT08 ) 
            F3 = FINDC( CSTASCC, TXCNT( 6 ), CHRT06 ) 
            F2 = FINDC( CSTASL , TXCNT( 5 ), CHRT05 ) 
            F1 = FINDC( TSCC   , TXCNT( 3 ), CHRT03 ) 
            F0 = FINDC( TSCCL  , TXCNT( 2 ), CHRT02 )

            IF( F5 .GT. 0 ) THEN
                ISRG = ISRG09( F5 ) 
                CALL SETSOURCE_GSURG
                CYCLE                       !  to end of sources-loop

            ELSEIF( F4 .GT. 0 ) THEN
                ISRG = ISRG08( F4 ) 
                CALL SETSOURCE_GSURG
                CYCLE                       !  to end of sources-loop

            ELSEIF( F3 .GT. 0 ) THEN
                ISRG = ISRG06( F3 ) 
                CALL SETSOURCE_GSURG
                CYCLE                       !  to end of sources-loop

            ELSEIF( F2 .GT. 0 ) THEN
                ISRG = ISRG05( F2 ) 
                CALL SETSOURCE_GSURG
                CYCLE                       !  to end of sources-loop

            ELSEIF( F1 .GT. 0 ) THEN
                ISRG = ISRG03( F1 ) 
                CALL SETSOURCE_GSURG
                CYCLE                       !  to end of sources-loop

            ELSEIF( F0 .GT. 0 ) THEN
                ISRG = ISRG02( F0 ) 
                CALL SETSOURCE_GSURG
                CYCLE                       !  to end of sources-loop

            END IF

C.............  Try for any FIPS code match
            F0 = FINDC( CFIP, TXCNT( 7 ), CHRT07 ) 

            IF( F0 .GT. 0 ) THEN
                ISRG = ISRG07( F0 ) 
                CALL SETSOURCE_GSURG
                CYCLE                       !  to end of sources-loop
            END IF

C.............  Try for any country/state code match (not, pol-specific)
            F0 = FINDC( CSTA, TXCNT( 4 ), CHRT04 ) 

            IF( F0 .GT. 0 ) THEN
                ISRG = ISRG04( F0 ) 
                CALL SETSOURCE_GSURG
                CYCLE                       !  to end of sources-loop
            END IF

            IF( ISRG01 .NE. IMISS3 .AND. REPDEFLT ) THEN
                ISRG = ISRG01
                    
                CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

                WRITE( MESG,94010 )
     &                 'WARNING: Using default gridding ' //
     &                 'cross-reference code of', ISRG, 'for:' //
     &                 CRLF() // BLANK10 // BUFFER( 1:L2 )
                CALL M3MESG( MESG )

                CALL SETSOURCE_GSURG

            ELSEIF( ISRG01 .NE. IMISS3 ) THEN
                ISRG = ISRG01
                CALL SETSOURCE_GSURG

            ELSE
                EFLAG = .TRUE.

                CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

                WRITE( MESG,94010 )
     &                 'ERROR: No gridding cross-reference ' //
     &                 'available (and no default) for:' //
     &                 CRLF() // BLANK10 // BUFFER( 1:L2 )

                CALL M3MESG( MESG )

            END IF    !  if default profile code is available or not

        END DO        !  end loop on source, S

        IF( EFLAG ) THEN
            MESG = 'Problem assigning gridding surrogates to sources'
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

C.............  This internal subprogram stores the index of the surrogate 
C               codes from the surrogates file for each source.
            SUBROUTINE SETSOURCE_GSURG

C.............  Local variables
            INTEGER          ISRGPOS  ! position of surrogate code in list

C----------------------------------------------------------------------
            SRGFIPIDX( S ) = S
            
            ISRGPOS = MAX( FIND1( ISRG, NSRGS, SRGLIST ), 0 )

            IF( ISRGPOS .EQ. 0 ) THEN

                CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &                 'ERROR: Gridding surrogate code', ISRG, 
     &                 'is not in surrogates file, but was ' //
     &                 CRLF() // BLANK5 // 'assigned to source:' //
     &                 CRLF() // BLANK10 // BUFFER( 1:L2 )
                CALL M3MESG( MESG )

            END IF
            
c            SRGIDPOS( S ) = ISRGPOS
            ASRGID( S ) = ISRG

            RETURN

C------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE SETSOURCE_GSURG

        END SUBROUTINE ASGNSURG
