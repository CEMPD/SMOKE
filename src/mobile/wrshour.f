
        SUBROUTINE WRSHOUR( FNAME, JDATE, JTIME, NCNTY, ARRAYPOS,
     &                      CNTYCODES, HOURTEMP )
   
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
        CHARACTER(*), INTENT (IN) :: FNAME          ! logical file name
        INTEGER     , INTENT (IN) :: JDATE          ! julian date
        INTEGER     , INTENT (IN) :: JTIME          ! time HHMMSS
        INTEGER     , INTENT (IN) :: NCNTY          ! no. counties
        INTEGER     , INTENT (IN) :: ARRAYPOS       ! position in temperature array
        INTEGER     , INTENT (IN) :: CNTYCODES( NCNTY )   ! county FIPS codes
        REAL        , INTENT (IN) :: HOURTEMP( NCNTY,24 ) ! hourly values

C...........   Local variables
        INTEGER         I       ! index variable

        CHARACTER*300   MESG    ! message buffer

        CHARACTER*16 :: PROGNAME = 'WRSHOUR' ! program name

C***********************************************************************
C   begin body of subroutine WRSHOUR

C.........  Write county codes to file
        IF( .NOT. WRITE3( FNAME, 'COUNTIES', JDATE, JTIME,
     &                    CNTYCODES ) ) THEN       
     	    MESG = 'Could not write county codes to "' //
     &              FNAME( 1:LEN_TRIM( FNAME ) ) // '".'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
        END IF

C.........  Write one hour of temperatures to file        
        IF( .NOT. WRITE3( FNAME, 'TKCOUNTY', JDATE, JTIME,  
     &                    HOURTEMP( :,ARRAYPOS ) ) ) THEN 
            MESG = 'Could not write hourly data to "' //
     &              FNAME( 1:LEN_TRIM( FNAME ) ) //  '".'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

        END SUBROUTINE WRSHOUR
