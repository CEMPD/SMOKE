
        SUBROUTINE ERRPKTS( PKTTYP, JT, JX, SKIPPOL, JTMAX, JXMAX, 
     &                      EFLAG )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      Check dimensions with maximum dimensions for packets, and reset 
C      maximum to actual in case any records were skipped.  Also, produce
C      messages if pollutant was skipped or if the number of usable packets
C      was zero.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 3/99 by M. Houyoux
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

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)    CRLF

        EXTERNAL CRLF

C...........   SUBROUTINE ARGUMENTS:

        CHARACTER(*), INTENT (IN) :: PKTTYP ! packet type 
        INTEGER     , INTENT (IN) :: JT     ! idx to control data tables
        INTEGER     , INTENT (IN) :: JX     ! idx to ungrouped cntl x-ref tables
        LOGICAL     , INTENT (IN) :: SKIPPOL! true: skipped rec is pol-spcfic
        INTEGER  , INTENT(IN OUT) :: JTMAX  ! max allowed JT
        INTEGER  , INTENT(IN OUT) :: JXMAX  ! max allowed JX
        LOGICAL     , INTENT(OUT) :: EFLAG  ! true: error occurred

C...........   Other local variables
        INTEGER         L        !  counters and indices

        CHARACTER(256)  MESG     ! message buffer

        CHARACTER(16) :: PROGNAME = 'ERRPKTS' ! program name

C***********************************************************************
C   begin body of subroutine ERRPKTS

C.........  Get length of packet name for messages
        L = LEN_TRIM( PKTTYP )

C.........  Write warning message for pollutants in file that are
C           not in master list
        IF( SKIPPOL ) THEN
            MESG = 'Pollutant-specific entries in the ' //
     &             PKTTYP // ' packet' // CRLF() // 
     &             BLANK10 // 'have been skipped.'
            CALL M3WARN( PROGNAME, 0, 0, MESG )
        END IF

C.........  Error for overflow of cross-reference information
        IF( JX .GT. JXMAX ) THEN

            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'INTERNAL ERROR: dimension for ' //
     &             'storing cross-reference for the ' //
     &             PKTTYP // CRLF() // BLANK10 // 
     &             'packet was', JXMAX, 'but actually needed', JX
            CALL M3MSG2( MESG )

C.........  Otherwise, store final count
        ELSE
            JXMAX = JX

        END IF

C.........  Error for overflow of control table information
        IF( JT .EQ. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'No usable ' // TRIM( PKTTYP ) // 
     &             ' control packet entries!'
            CALL M3MSG2( MESG )

        ELSE IF( JT .GT. JTMAX ) THEN

            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'INTERNAL ERROR: dimension for ' //
     &             'storing control data for the ' //
     &             PKTTYP // CRLF() // BLANK10 // 
     &             'packet was', JTMAX, 'but actually needed', JT
            CALL M3MSG2( MESG )

C.........  Otherwise, store final count
        ELSE
            JTMAX = JT

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE ERRPKTS
