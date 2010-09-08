
        SUBROUTINE HOURMET( NSRC, AVGTYPE, JDATE, JTIME, DAYBEGT,   
     &              SKIPDATA, LDAYSAV, RH_STRHR, RH_ENDHR )

C***********************************************************************
C  subroutine body starts at line 92
C
C  DESCRIPTION:
C       Creates summed hourly temperatures by county. Checks that temperatures
C       are within requested minimum and maximum. Keeps track of total number
C       of sources for averaging later.
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

C...........   MODULES for public variables

C...........   This module is the derived meteorology data for emission factors
        USE MODMET, ONLY: TASRC, QVSRC, PRESSRC, TKHOUR, RHHOUR, NDAYSRC
        
        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS
        CHARACTER(2) CRLF
        INTEGER      ENVINT
        LOGICAL      ISDSTIME
        REAL         CALCRELHUM

        EXTERNAL     CRLF, ENVINT, ISDSTIME, CALCRELHUM
                
C...........   SUBROUTINE ARGUMENTS
        INTEGER,      INTENT    (IN) :: NSRC                  ! no. sources
        CHARACTER(*), INTENT    (IN) :: AVGTYPE               ! avg type
        INTEGER,      INTENT    (IN) :: JDATE                 ! YYYYDDD
        INTEGER,      INTENT    (IN) :: JTIME                 ! HHMMSS
        INTEGER,      INTENT    (IN) :: DAYBEGT ( NSRC )      ! begin. time for day
        LOGICAL,      INTENT    (IN) :: SKIPDATA              ! skip data for non-day averaging
        LOGICAL,      INTENT    (IN) :: LDAYSAV ( NSRC )      ! true: use daylight time
        INTEGER,      INTENT    (IN) :: RH_STRHR              ! HHMMSS
        INTEGER,      INTENT    (IN) :: RH_ENDHR              ! HHMMSS

C...........   Other local variables
        INTEGER     L, S        ! counters and indices
        INTEGER     IOS         ! I/O status
        INTEGER     TIMESLOT    ! array location
        
        INTEGER, SAVE :: MXWARN ! maximum number of warnings
        INTEGER, SAVE :: NWARN  ! total number of warnings printed

        REAL        TEMPVAL     ! temperature value in Farenheight
        REAL        MIXVAL      ! mixing ratio value
        REAL        PRESVAL     ! pressure value
        REAL        RHVAL       ! RH value

        LOGICAL       :: DAYLIT  = .FALSE.  ! true: date is daylight savings
        LOGICAL, SAVE :: INITIAL = .TRUE.   ! true: first time

        CHARACTER(300)     BUFFER    ! formatted source info for messages
        CHARACTER(300)     MESG      ! message buffer
 
        CHARACTER(16) :: PROGNAME = 'HOURMET' ! program name

C***********************************************************************
C   begin body of subroutine HOURTEMP

C.........  For the first time, initialize all entries to zero
        IF( INITIAL ) THEN
            TKHOUR = 0.  ! array
            RHHOUR = 0.
            
C.............  Get maximum number of warnings
            MXWARN = ENVINT( WARNSET, ' ', 100, IOS )
            NWARN = 0
            
            INITIAL = .FALSE.
        END IF
C.........  Loop through sources
        DO S = 1, NSRC

C.............  Calculate time slot in output array for this time step
C               Appropriate 24 hour time will be day starting time (12 AM in local 
C               time zone ) subtracted from met data time (in GMT)
            TIMESLOT = 1 + ( JTIME - DAYBEGT( S ) ) / 10000 

C.............  Restore daylight saving time if necessay
            DAYLIT = ISDSTIME( JDATE )
            IF( DAYLIT .AND. LDAYSAV( S ) ) THEN
                TIMESLOT = TIMESLOT - 1      ! substract 1hr on DST date
            END IF
                
C.............  If timeslot is less than zero, add 24; if better data comes
C               along, the old data will get overwritten (helps in case of
C               one running one day)
            IF( TIMESLOT <= 0 ) THEN
                TIMESLOT = TIMESLOT + 24
            END IF

            TEMPVAL = TASRC( S )   ! unit in Kelvin
            MIXVAL  = QVSRC( S )
            PRESVAL = PRESSRC( S )

C.............  only estimate RH based on user-defined averaging start/end hour
            IF( TIMESLOT >= RH_STRHR .AND. TIMESLOT <= RH_ENDHR ) THEN
                RHVAL = CALCRELHUM( TEMPVAL, PRESVAL, MIXVAL )
            ELSE
                RHVAL = 0.0
            END IF

C.............  Convert K to F degree for temperature
            TEMPVAL = 1.8 * TEMPVAL - 459.67

            IF( TEMPVAL > AMISS3 )THEN

C.................  Store values in hourly arrays                
                IF( AVGTYPE == 'DAILY' ) THEN
                    TKHOUR ( S,TIMESLOT ) = TEMPVAL
                    RHHOUR ( S,TIMESLOT ) = RHVAL
                    NDAYSRC( S,TIMESLOT ) = 1
                ELSE
                    IF( .NOT. SKIPDATA ) THEN
                        TKHOUR( S,TIMESLOT ) =
     &                                    TKHOUR( S,TIMESLOT ) + TEMPVAL

                        RHHOUR( S,TIMESLOT ) = 
     &                                    RHHOUR( S,TIMESLOT ) + RHVAL

                        NDAYSRC( S,TIMESLOT ) = 
     &                                    NDAYSRC( S,TIMESLOT ) + 1
                    END IF
                END IF
c        print*,s,timeslot,TEMPVAL,RHVAL,NDAYSRC(S,TIMESLOT),'S TimeSlot BH'

            END IF

        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

94020   FORMAT( A, 4( 1X, F8.2, 1X, A ) )
 
        END SUBROUTINE HOURMET
