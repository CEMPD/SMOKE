
        SUBROUTINE HOURTEMP( NSRC, JTIME, DAYBEGT, SRCARRAY, 
     &                       MINTEMP, MAXTEMP, SKIPDATA, NDAYSRC )

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

C...........   This module is used for mobile setup
        USE MODMBSET, ONLY: DAILY
        
C...........   This module is the source inventory arrays
        USE MODSOURC, ONLY: CSOURC

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NCHARS
        
C...........   This module is the derived meteorology data for emission factors
        USE MODMET, ONLY: TASRC, QVSRC, PRESSRC, TKHOUR, QVHOUR, BPHOUR
        
        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS
        CHARACTER(2) CRLF
        INTEGER      ENVINT

        EXTERNAL     CRLF, ENVINT
                
C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT    (IN) :: NSRC                  ! no. sources
        INTEGER, INTENT    (IN) :: JTIME                 ! HHMMSS
        INTEGER, INTENT    (IN) :: DAYBEGT ( NSRC )      ! begin. time for day
        INTEGER, INTENT    (IN) :: SRCARRAY( NSRC,2 )    ! per-source FIPS and averaging type
        REAL   , INTENT    (IN) :: MINTEMP               ! minimum temperature
        REAL   , INTENT    (IN) :: MAXTEMP               ! maximum temperature
        LOGICAL, INTENT    (IN) :: SKIPDATA             ! skip data for non-day averaging
        INTEGER, INTENT(IN OUT) :: NDAYSRC ( NSRC,24 )   ! no. days to average per source

C...........   Other local variables
        INTEGER     L, S        ! counters and indices
        INTEGER     IOS         ! I/O status
        INTEGER     TIMESLOT    ! array location
        
        INTEGER, SAVE :: MXWARN ! maximum number of warnings
        INTEGER, SAVE :: NWARN  ! total number of warnings printed

        REAL        TEMPVAL     ! temperature value
        REAL        MIXVAL      ! mixing ratio value
        REAL        PRESVAL     ! pressure value

        LOGICAL, SAVE :: INITIAL = .TRUE.  ! true: first time

        CHARACTER(300)     BUFFER    ! formatted source info for messages
        CHARACTER(300)     MESG      ! message buffer
        CHARACTER(SRCLEN3) CSRC      ! tmp concat source characteristics
 
       CHARACTER(16) :: PROGNAME = 'HOURTEMP' ! program name

C***********************************************************************
C   begin body of subroutine HOURTEMP

C.........  For the first time, initialize all entries to zero
        IF( INITIAL ) THEN
            TKHOUR = 0.  ! array
            QVHOUR = 0.
            BPHOUR = 0.
            
C.............  Get maximum number of warnings
            MXWARN = ENVINT( WARNSET, ' ', 100, IOS )
            NWARN = 0
            
            INITIAL = .FALSE.
        END IF

C.........  Loop through sources
        DO S = 1, NSRC

C.............  Calculate time slot in output array for this time step
C               Appropriate 24 hour time will be day starting time (6 AM in local 
C               time zone ) subtracted from met data time (in GMT)
            TIMESLOT = 1 + ( JTIME - DAYBEGT( S ) ) / 10000 
                
C.............  If timeslot is less than zero, add 24; if better data comes
C               along, the old data will get overwritten (helps in case of
C               one running one day)
            IF( TIMESLOT <= 0 ) THEN
                TIMESLOT = TIMESLOT + 24
            END IF
                
            CSRC = CSOURC( S )
            
            TEMPVAL = TASRC( S )
            MIXVAL  = QVSRC( S )
            PRESVAL = PRESSRC( S )

            IF( TEMPVAL > AMISS3 )THEN

C.................  Check that temperature is within min and max bounds
                IF( TEMPVAL < MINTEMP ) THEN

                    IF( NWARN <= MXWARN ) THEN
                        CALL FMTCSRC( CSRC, NCHARS, BUFFER, L )
                        WRITE( MESG, 94020 )
     &                     'Increasing hourly temperature from',
     &                     TEMPVAL, 'to', MINTEMP, 'for source' //
     &                     CRLF() // BLANK10 // BUFFER( 1:L ) // '.'
                        CALL M3MESG( MESG )
                        NWARN = NWARN + 1
                    END IF

C.....................  Set value to minimum
                    TEMPVAL = MINTEMP 

                ELSEIF( TEMPVAL > MAXTEMP ) THEN

                    IF( NWARN <= MXWARN ) THEN
                        CALL FMTCSRC( CSRC, NCHARS, BUFFER, L )
                        WRITE( MESG, 94020 )
     &                     'Decreasing hourly temperature from',
     &                     TEMPVAL, ' to', MAXTEMP, 'for source' //
     &                     CRLF() // BLANK10 // BUFFER( 1:L ) // '.'
                        CALL M3MESG( MESG )
                        NWARN = NWARN + 1
                    END IF

C.....................  Set value to maximum
                    TEMPVAL = MAXTEMP
                
                END IF

C.................  Store values in hourly arrays                
                IF( SRCARRAY( S,2 ) == DAILY ) THEN
                    TKHOUR ( S,TIMESLOT ) = TEMPVAL
                    QVHOUR ( S,TIMESLOT ) = MIXVAL
                    BPHOUR ( S,TIMESLOT ) = PRESVAL
                    NDAYSRC( S,TIMESLOT ) = 1
                ELSE
                    IF( .NOT. SKIPDATA ) THEN
                        TKHOUR( S,TIMESLOT ) =
     &                                    TKHOUR( S,TIMESLOT ) + TEMPVAL
                        QVHOUR( S,TIMESLOT ) =
     &                                    QVHOUR( S,TIMESLOT ) + MIXVAL
                        BPHOUR( S,TIMESLOT ) =
     &                                    BPHOUR( S,TIMESLOT ) + PRESVAL
                        NDAYSRC( S,TIMESLOT ) = 
     &                                    NDAYSRC( S,TIMESLOT ) + 1
                    END IF
                END IF

            END IF

        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

94020   FORMAT( A, 4( 1X, F8.2, 1X, A ) )
 
        END SUBROUTINE HOURTEMP
