
        SUBROUTINE FILLCNTL( PKTTYP, JTMAX, JXMAX, 
     &                       PKTINFO, JPOL, JT, JX )

C***********************************************************************
C  subroutine body starts at line 
C      This routine increments JX and JT and stores the control data tables 
C      and ungrouped control x-ref tables
C
C  DESCRIPTION:
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 3/99 by M. Houyoux
C
C****************************************************************************/
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
C.........  This module is for cross reference tables
        USE MODXREF, ONLY: INDXTA, ISPTA, CSCCTA, CSRCTA, MPRNA, CMACTA,
     &                     CISICA

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'CPKTDAT.EXT'   !  control packet contents

C...........   EXTERNAL FUNCTIONS:
        INTEGER       STR2INT
        EXTERNAL      STR2INT

C...........   SUBROUTINE ARGUMENTS:
        INTEGER     , INTENT (IN) :: PKTTYP ! packet type number
        INTEGER     , INTENT (IN) :: JTMAX  ! max allowed JT
        INTEGER     , INTENT (IN) :: JXMAX  ! max allowed JX
        TYPE( CPACKET ),INTENT(IN):: PKTINFO! packet information
        INTEGER     , INTENT (IN) :: JPOL   ! pollutant position
        INTEGER  , INTENT(IN OUT) :: JT     ! idx to control data tables
        INTEGER  , INTENT(IN OUT) :: JX     ! idx to ungrouped cntl x-ref tables

C...........   Other local variables
        INTEGER         N               !  counters and indices
        INTEGER         SIC             !  tmp SIC

        CHARACTER(SCCLEN3) TMPSCC   !  tmp SCC
        CHARACTER(ALLLEN3) CSRCALL  !  buffer for source char, incl pol

        CHARACTER(16) :: PROGNAME = 'FILLCNTL'   ! program name

C***********************************************************************
C   begin body of subroutine FILLCNTL

C.........  Increment data table counter
        JT = JT + 1

C.........  Store packet information in temporary variables
        TMPSCC = PKTINFO%TSCC
        IF( PKTINFO%CSIC /= REPEAT( '0', SICLEN3 ) ) THEN
            TMPSCC = PKTINFO%CSIC
        END IF

        JX = JX + 1

        IF( JX .GT. JXMAX ) RETURN  ! to next iteration

C.........  Store unsorted x-ref table entries
        INDXTA( JX ) = JX
        ISPTA ( JX ) = JPOL

C.........  Parse the line of data into segments based on the rules
C.........  Ensure that pollutant is in master list of pollutants or
C               skip the pollutant-specific entry
        CSCCTA( JX ) = PKTINFO%TSCC
        CMACTA( JX ) = PKTINFO%CMCT
        CISICA( JX ) = PKTINFO%CSIC
        MPRNA ( JX ) = JT   ! Position in data table

C.........  Store sorting criteria as right-justified in fields
        CSRCALL = ' '
        SELECT CASE( CATEGORY )

        CASE( 'AREA' )
            CALL BLDCSRC( PKTINFO%CFIP, PLTBLNK3, CHRBLNK3,
     &                    CHRBLNK3, CHRBLNK3, CHRBLNK3, 
     &                    CHRBLNK3, POLBLNK3, CSRCALL )

            CSRCTA( JX ) = CSRCALL( 1:SRCLEN3 ) // TMPSCC // 
     &                     PKTINFO%CMCT // PKTINFO%CPOS

        CASE( 'MOBILE' )
            CALL BLDCSRC( PKTINFO%CFIP, PLTBLNK3, CHRBLNK3,
     &                    CHRBLNK3, CHRBLNK3, CHRBLNK3, 
     &                    CHRBLNK3, POLBLNK3, CSRCALL )

            CSRCTA( JX ) = CSRCALL( 1:SRCLEN3 ) // TMPSCC // 
     &                     PKTINFO%CPOS

        CASE( 'POINT' )

            CALL BLDCSRC( PKTINFO%CFIP, PKTINFO%PLT, PKTINFO%CHAR1,
     &                    PKTINFO%CHAR2, PKTINFO%CHAR3, 
     &                    PKTINFO%CHAR4, PKTINFO%CHAR5, POLBLNK3, 
     &                    CSRCALL )

            CSRCTA( JX ) = CSRCALL( 1:SRCLEN3 ) // TMPSCC // 
     &                     PKTINFO%CMCT // PKTINFO%CPOS

        END SELECT

C.........  Now store the packet information in the packet tables...

C.........  Double check for memory allocation
        IF( JT .GT. JTMAX ) RETURN

C.........  Store control data table, depending on type of packet
        CALL FILLCDAT( PKTLIST( PKTTYP ), JT, PKTINFO )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010  FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE FILLCNTL
