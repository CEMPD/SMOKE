
        SUBROUTINE FILLCNTL( PKTTYP, JTMAX, JXMAX, NSTART, NEND, 
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

C.........  MODULES for public variables
C.........  This module is for cross reference tables
        USE MODXREF

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'CPKTDAT.EXT'   !  control packet contents

C...........   SUBROUTINE ARGUMENTS:

        CHARACTER(*), INTENT (IN) :: PKTTYP ! packet type 
        INTEGER     , INTENT (IN) :: JTMAX  ! max allowed JT
        INTEGER     , INTENT (IN) :: JXMAX  ! max allowed JX
        INTEGER     , INTENT (IN) :: NSTART ! loop beg for expding SIC-based rec
        INTEGER     , INTENT (IN) :: NEND   ! loop end for expding SIC-based rec
        TYPE( CPACKET ),INTENT(IN):: PKTINFO! packet information
        INTEGER     , INTENT (IN) :: JPOL   ! pollutant position
        INTEGER  , INTENT(IN OUT) :: JT     ! idx to control data tables
        INTEGER  , INTENT(IN OUT) :: JX     ! idx to ungrouped cntl x-ref tables

C...........   Other local variables
        INTEGER         N               !  counters and indices

        CHARACTER(LEN=SCCLEN3) TMPSCC   !  tmp SCC
        CHARACTER(LEN=ALLLEN3) CSRCALL  !  buffer for source char, incl pol

        CHARACTER*16 :: PROGNAME = 'FILLCNTL'   ! program name

C***********************************************************************
C   begin body of subroutine FILLCNTL

C.........  Increment data table counter
        JT = JT + 1

C.........  Loop through expansion records and store all SCCs for SIC, or
C           the same SCC if no expansion
        DO N = NSTART, NEND

            JX = JX + 1

            IF( JX .GT. JXMAX ) CYCLE  ! to next iteration

            TMPSCC = INVSCC( N )

C.............  Store unsorted x-ref table entries
            INDXTA( JX ) = JX
            ISPTA ( JX ) = JPOL

C.........  Parse the line of data into segments based on the rules
C.........  Ensure that pollutant is in master list of pollutants or
C           skip the pollutant-specific entry
            CSCCTA( JX ) = TMPSCC
            MPRNA ( JX ) = JT   ! Position in data table

C.............  Store sorting criteria as right-justified in fields
            CSRCALL = ' '
            CALL BLDCSRC( PKTINFO%CFIP, PKTINFO%PLT, PKTINFO%CHAR1,
     &                    PKTINFO%CHAR2, PKTINFO%CHAR3, PKTINFO%CHAR4,
     &                    PKTINFO%CHAR5, POLBLNK3, CSRCALL )

            CSRCTA( JX ) = CSRCALL( 1:SRCLEN3 )// TMPSCC// PKTINFO%CPOS

        END DO  ! End loop through the expansion records

C.........  Now store the packet information in the packet tables...

C.........  Double check for memory allocation
        IF( JT .GT. JTMAX ) RETURN

C.........  Store control data table, depending on type of packet
        CALL FILLCDAT( PKTTYP, JT, PKTINFO )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010  FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE FILLCNTL
