
        INTEGER FUNCTION GETVMIX( CFIP, CRWT, CLNK )

C***********************************************************************
C  function body starts at line 
C
C  DESCRIPTION:
C     This function assigns the vehicle mix to a source.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 2/2000 by M. Houyoux
C
C**************************************************************************
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

C.........  MODULES for public variables
C...........   This module is for cross reference tables
        USE MODXREF, ONLY: CHRT02, CHRT04, CHRT05, CHRT07, CHRT08,
     &                     CHRT10, CHRT11,
     &                     IMVS01, IMVS02, IMVS04, IMVS05, IMVS07,
     &                     IMVS08, IMVS10, IMVS11,
     &                     TXCNT

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: LSCCEND, NCHARS

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        LOGICAL         ENVYN
        INTEGER         FINDC

        EXTERNAL  CRLF, ENVYN, FINDC

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: CFIP   ! co/st/cy code
        CHARACTER(*), INTENT (IN) :: CRWT   ! roadway type
        CHARACTER(*), INTENT (IN) :: CLNK   ! link code

C...........   Local parameters

C...........   Other local variables
        INTEGER         IOS          ! i/o status
        INTEGER         L2
        INTEGER         F6, F5, F4, F3, F2, F1, F0

        LOGICAL     :: FIRSTIME = .TRUE.  ! true: first time routine is called
        LOGICAL     :: REPDEFLT = .TRUE.  ! true: report default applied

        CHARACTER*300            BUFFER   !  source characteristics
        CHARACTER*300            MESG     !  message buffer

        CHARACTER(LEN=SRCLEN3)   CSRC     ! tmp source chars string
        CHARACTER(LEN=STALEN3)   CSTA     ! tmp Country/state code
        CHARACTER(LEN=SCCLEN3)   TSCC     ! tmp 10-digit SCC
        CHARACTER(LEN=SCCLEN3)   TSCCL    ! tmp left digits of TSCC
        CHARACTER(LEN=SS0LEN3):: CHK11=' '! tmp FIPS // Plant // SCC
        CHARACTER(LEN=FPLLEN3):: CHK10=' '! tmp FIPS code // plant id
        CHARACTER(LEN=FPSLEN3):: CHK08=' '! tmp FIPS code // left SCC
        CHARACTER(LEN=STSLEN3):: CHK05=' '! tmp Country/state code // left SCC
        CHARACTER(LEN=VIDLEN3), SAVE :: VIDZERO  ! zero vehicle type

        CHARACTER*16 :: PROGNAME = 'GETVMIX' ! program name

C***********************************************************************
C   begin body of subroutine GETVMIX

C.........  The first time routine is called...
        IF( FIRSTIME ) THEN

C.............  Initialize zero vehicle type ID
            VIDZERO = REPEAT( '0', VIDLEN3 )

C.............  Retrieve environment variables
            MESG = 'Switch for reporting default temporal profiles'
            REPDEFLT = ENVYN ( 'REPORT_DEFAULTS', MESG, .TRUE., IOS )

            FIRSTIME = .FALSE.

        END IF

C.........  Adjust roadway type to internal SCC format for matching
        TSCC = CRWT // VIDZERO
        CALL PADZERO( TSCC )

C.........  Create source characteristics full string
        CALL BLDCSRC( CFIP, RWTBLNK3, CLNK, CHRBLNK3,
     &                    CHRBLNK3, CHRBLNK3, CHRBLNK3,
     &                    POLBLNK3, CSRC )

C.........  Define various source characteristic combinations
C.........  NOTE - the full-SCC matches are the same as the left-SCC matches
C           because the vehicle mix fields are not used to assign the VMT mix
C           (by definition of the data in the VMT mix file).
        TSCCL   = TSCC( 1:LSCCEND )

        CHK11  = CSRC( 1:MBENDL3(3) )// TSCC          ! Cnty// RWT// LNK
        CHK10  = CSRC( 1:MBENDL3(3) )                 ! Cnty// LNK
        CHK08  = CFIP // TSCCL                        ! County// RWT
        CHK05  = CSTA // TSCCL                        ! State // road type

C.........  Try for FIPS code & roadway type & link match; then
C                   FIPS code & link match; then
C                   FIPS code & roadway type match; then
C                   Cy/st code & roadway type; then
C                   FIPS code; then
C                   Cy/st code; then
C                   roadway type
        
        F6 = FINDC( CHK11, TXCNT( 11 ), CHRT11 ) 
        F5 = FINDC( CHK10, TXCNT( 10 ), CHRT10 ) 
        F4 = FINDC( CHK08, TXCNT( 8  ), CHRT08 ) 
        F3 = FINDC( CHK05, TXCNT( 5  ), CHRT05 ) 
        F2 = FINDC( CFIP , TXCNT( 7  ), CHRT07 ) 
        F1 = FINDC( CSTA , TXCNT( 4  ), CHRT04 ) 
        F0 = FINDC( TSCCL, TXCNT( 2  ), CHRT02 )

        IF( F6 .GT. 0 ) THEN
            GETVMIX = IMVS11( F6 ) 
            RETURN 

        ELSEIF( F5 .GT. 0 ) THEN
            GETVMIX = IMVS10( F5 ) 
            RETURN 

        ELSEIF( F4 .GT. 0 ) THEN
            GETVMIX = IMVS08( F4 ) 
            RETURN

        ELSEIF( F3 .GT. 0 ) THEN
            GETVMIX = IMVS05( F3 ) 
            RETURN

        ELSEIF( F2 .GT. 0 ) THEN
            GETVMIX = IMVS07( F2 ) 
            RETURN

        ELSEIF( F1 .GT. 0 ) THEN
            GETVMIX = IMVS04( F1 ) 
            RETURN

        ELSEIF( F0 .GT. 0 ) THEN
            GETVMIX = IMVS02( F0 ) 
            RETURN

        END IF

C.........  Apply default and report it
        IF( IMVS01 .NE. IMISS3 .AND. REPDEFLT ) THEN
            GETVMIX = IMVS01

            CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

            WRITE( MESG,94010 )
     &                 'NOTE: Using default gridding ' //
     &                 'cross-reference for:' //
     &                 CRLF() // BLANK10 // BUFFER( 1:L2 )
            CALL M3MESG( MESG )

C.........  Apply default with no report
        ELSEIF( IMVS01 .NE. IMISS3 ) THEN
            GETVMIX = IMVS01            

C.........  Report that no match was found - this is an error
        ELSE
            GETVMIX = -1

            CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

            WRITE( MESG,94010 )
     &             'ERROR: No VMT mix profile ' //
     &             'available (and no default) for:' //
     &             CRLF() // BLANK10 // BUFFER( 1:L2 )

            CALL M3MESG( MESG )

        END IF    !  if default profile code is available or not

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END FUNCTION GETVMIX
