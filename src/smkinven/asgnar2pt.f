
        SUBROUTINE ASGNAR2PT( NRAWBP )

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
C
C************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC  
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
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
        INTEGER         FINDC

        EXTERNAL        CRLF, FINDC

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: NRAWBP  ! no. raw records by pollutant

C.........  Other local variables
        INTEGER          I, J, LS, S    !  counters and indices

        INTEGER          IOS         ! i/o status
        INTEGER          F4, F5      ! tmp find indices

        LOGICAL       :: EFLAG    = .FALSE.

        CHARACTER*256             MESG     ! message buffer
        CHARACTER(LEN=SRCLEN3)    CSRC     ! tmp source chars string
        CHARACTER(LEN=FIPLEN3)    CFIP     ! tmp (character) FIPS code
        CHARACTER(LEN=FPSLEN3)    CFIPSCC  ! tmp FIPS code // SCC
        CHARACTER(LEN=FPSLEN3)    CFIPSL   ! tmp FIPS code // left SCC
        CHARACTER(LEN=SCCLEN3)    TSCC     ! tmp 10-digit SCC
        CHARACTER(LEN=SCCLEN3)    TSCCL    ! tmp left digits of TSCC

        CHARACTER*16 :: PROGNAME = 'ASGNAR2PT' ! program name

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

C.........  Loop through the unsorted sources
        LS = 0
        DO I = 1, NRAWBP

            J = INDEXA( I )
            S = SRCIDA( I )

            IF ( S .EQ. LS ) CYCLE

c NOTE: Perhaps change this to be a generic routine for both cross-reference need?  The
C N: only problem with this is that CHRT* will be shared by both steps :(

C.............  Create selection 
            SELECT CASE ( CATEGORY )

            CASE ( 'AREA', 'MOBILE' )
                CSRC    = CSOURCA( J )
                TSCC    = CSCCA  ( J )
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

            LS = S

        END DO        !  end loop on source x pollutants

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE ASGNAR2PT
