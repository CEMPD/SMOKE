
        SUBROUTINE WRSHOUR( NSRC, JDATE, JTIME, FNAME, VHOUR )
   
C***********************************************************************
C  subroutine WRSHOUR body starts at line < >
C
C  DESCRIPTION:
C      Write per-source hourly temperature data
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION HISTORY:
C
C***************************************************************************
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
C****************************************************************************

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NSRC           ! no. sources
        INTEGER     , INTENT (IN) :: JDATE          ! julian date
        INTEGER     , INTENT (IN) :: JTIME          ! time HHMMSS
        CHARACTER(*), INTENT (IN) :: FNAME          ! logical file name
        REAL        , INTENT (IN) :: VHOUR( NSRC,0:23 ) ! hourly value

C...........   Local variables
        INTEGER         I       ! index variable

        CHARACTER*300   MESG    ! message buffer

        CHARACTER*16 :: PROGNAME = 'WRSHOUR' ! program name

C***********************************************************************
C   begin body of subroutine WRSHOUR

C...........   Loop through 24 hours
        DO I = 0, 23

C.................  Write out temperature information for the current hour

            IF( .NOT. WRITE3( FNAME, 'TKHOUR', JDATE, JTIME, 
     &                        VHOUR(:,I) ) ) THEN 

                MESG = 'Could not write hourly data to "' //
     &                  FNAME( 1:LEN_TRIM( FNAME ) ) //  '".'

                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

            END IF

C..................  Increment output time
            CALL NEXTIME( JDATE, JTIME, 10000 )

        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

        END SUBROUTINE WRSHOUR
