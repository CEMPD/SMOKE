
        SUBROUTINE WRDAYMSG( JDATE, MESG )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      Writes a text message to stdout and log file about which day is
C      being processed
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 9/99 by M. Houyoux
C
C**************************************************************************
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

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS:
        CHARACTER*14    MMDDYY
        INTEGER         WKDAY

        EXTERNAL        MMDDYY, WKDAY

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: JDATE    ! Julian date
        CHARACTER(*), INTENT(OUT) :: MESG     ! message buffer

C...........   Local variables
        INTEGER         DAY          !  day of week number
        INTEGER         L            !  length of day name

        CHARACTER*16 :: PROGNAME = 'WRDAYMSG' ! program name

C***********************************************************************
C   begin body of subroutine WRDAYMSG

        DAY = WKDAY( JDATE )

        L = LEN_TRIM( DAYS( DAY ) )
        MESG= 'Processing '// DAYS( DAY )( 1:L )// ' '// MMDDYY( JDATE )
        CALL M3MSG2( MESG )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

        END SUBROUTINE WRDAYMSG

