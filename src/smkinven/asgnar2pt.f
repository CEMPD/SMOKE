
        SUBROUTINE ASGNAR2PT

C***********************************************************************
C  subroutine body starts at line 109
C
C  DESCRIPTION:
C      For each source, find the most specific area-to-point adjustment
C      that applies to that source. Do this using the grouped tables of
C      gridding cross references from RDAR2PT (and stored in MODLISTS).  
C      Only one group (full FIPS & full SCC) has been defined at this time
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 11/02 by M. Houyoux
C     09/2025 by HT UNC-IE:  Use M3UTILIO
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
C***************************************************************************
        USE M3UTILIO

C...........   MODULES for public variables   
C...........   This module contains the source ararys
        USE MODSOURC, ONLY: CSOURC, CSCC

C...........   This module contains the cross-reference tables
        USE MODXREF, ONLY: AR2PTTBL, AR2PTIDX, AR2PTCNT, TXCNT,
     &                     CHRT09, CHRT08, ARPT09, ARPT08

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, NSRC, LSCCEND

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C       INCLUDE 'PARMS3.EXT'    !  i/o api constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
C       CHARACTER(2)    CRLF
C       INTEGER         FINDC
C       LOGICAL         SETSCCTYPE

C       EXTERNAL        CRLF, FINDC, SETSCCTYPE
        LOGICAL, EXTERNAL :: SETSCCTYPE

C...........   SUBROUTINE ARGUMENTS

C.........  Other local variables
        INTEGER          S           !  counters and indices

        INTEGER          IOS         ! i/o status
        INTEGER          F4, F5      ! tmp find indices

        LOGICAL       :: EFLAG    = .FALSE.
        LOGICAL          SCCFLAG           ! true: SCC type is different from previous

        CHARACTER(256)        MESG     ! message buffer
        CHARACTER(SRCLEN3)    CSRC     ! tmp source chars string
        CHARACTER(FIPLEN3)    CFIP     ! tmp (character) FIPS code
        CHARACTER(FPSLEN3)    CFIPSCC  ! tmp FIPS code // SCC
        CHARACTER(FPSLEN3)    CFIPSL   ! tmp FIPS code // left SCC
        CHARACTER(SCCLEN3)    TSCC     ! tmp 10-digit SCC
        CHARACTER(SCCLEN3)    TSCCL    ! tmp left digits of TSCC

        CHARACTER(16) :: PROGNAME = 'ASGNAR2PT' ! program name

C***********************************************************************
C   begin body of subroutine ASGNAR2PT

C.........  Allocate memory for source-based arrays
        ALLOCATE( AR2PTTBL( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AR2PTTBL', PROGNAME )
        ALLOCATE( AR2PTIDX( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AR2PTIDX', PROGNAME )
        ALLOCATE( AR2PTCNT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AR2PTCNT', PROGNAME )
        
C.........  Initialize arrays
        AR2PTTBL = 0  ! array
        AR2PTIDX = 0  ! array
        AR2PTCNT = 0  ! array

C.........  Loop through the sorted sources
        DO S = 1, NSRC

c NOTE: Perhaps change this to be a generic routine for both cross-reference need?  The
C N: only problem with this is that CHRT* will be shared by both steps :(

C.............  Create selection 
            SELECT CASE ( CATEGORY )

            CASE ( 'AREA', 'MOBILE' )
                CSRC    = CSOURC( S )
                TSCC    = CSCC( S )
                
C.................  Set type of SCC                
                SCCFLAG = SETSCCTYPE( TSCC )
                TSCCL   = TSCC( 1:LSCCEND )
                
                CFIP    = CSRC( 1:FIPLEN3 )
                CFIPSCC = CFIP // TSCC
                CFIPSL  = CFIP // TSCCL
            END SELECT

C.................  Try for FIPS code & SCC match
            F5 = FINDC( CFIPSCC, TXCNT( 9 ), CHRT09 ) 
            F4 = FINDC( CFIPSL , TXCNT( 8 ), CHRT08 ) 

            IF( F5 .GT. 0 ) THEN
                AR2PTTBL( S ) = ARPT09( F5,1 )
                AR2PTIDX( S ) = ARPT09( F5,2 )
                AR2PTCNT( S ) = ARPT09( F5,3 )

            ELSE IF( F4 .GT. 0 ) THEN
                AR2PTTBL( S ) = ARPT08( F4,1 )
                AR2PTIDX( S ) = ARPT08( F4,2 )
                AR2PTCNT( S ) = ARPT08( F4,3 )

            END IF

        END DO        !  end loop on sources

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE ASGNAR2PT
