
        SUBROUTINE WRSHOUR( FNAME, JDATE, JTIME, NCNTY,
     &                      CNTYCODES, HOURTEMP, HOURQV, HOURBP )
   
C***********************************************************************
C  subroutine WRSHOUR body starts at line 62
C
C  DESCRIPTION:
C      Write by county hourly temperature, relative humidity, and
C      barometric pressure data
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
        INTEGER     , INTENT (IN) :: CNTYCODES( NCNTY )   ! county FIPS codes
        REAL        , INTENT (IN) :: HOURTEMP( NCNTY )    ! hourly temperatures
        REAL        , INTENT (IN) :: HOURQV( NCNTY )      ! hourly mixing ratios
        REAL        , INTENT (IN) :: HOURBP( NCNTY )      ! hourly barometric pressures

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
     &              TRIM( FNAME ) // '".'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
        END IF

C.........  Write one hour of temperatures to file        
        IF( .NOT. WRITE3( FNAME, 'TKCOUNTY', JDATE, JTIME,  
     &                    HOURTEMP( : ) ) ) THEN 
            MESG = 'Could not write hourly temperatures to "' //
     &              TRIM( FNAME ) //  '".'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
        END IF

C.........  Write one hour of mixing ratios to file
        IF( .NOT. WRITE3( FNAME, 'QVCOUNTY', JDATE, JTIME,
     &                    HOURQV( : ) ) ) THEN
            MESG = 'Could not write hourly mixing ratio data ' //
     &             'to "' // TRIM( FNAME ) // '".'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
        END IF
        
C.........  Write one hour of barometric pressure data to file
        IF( .NOT. WRITE3( FNAME, 'BPCOUNTY', JDATE, JTIME,
     &                    HOURBP( : ) ) ) THEN
            MESG = 'Could not write hourly barometric pressure data ' //
     &             'to "' // TRIM( FNAME ) // '".'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

        END SUBROUTINE WRSHOUR
