
        SUBROUTINE HOURTEMP( NSRC, JTIME, DAYBEGT, VALBYSRC, SRCARRAY, 
     &                       MINTEMP, MAXTEMP, SKIPDATA, NDAYSRC, 
     &                       HOUROUT )

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
C****************************************************************************

C...........   MODULES for public variables

C...........   This module is used for mobile setup
        USE MODMBSET, ONLY: DAILY
        
C...........   This module is the source inventory arrays
        USE MODSOURC, ONLY: CSOURC

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NCHARS
        
        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS 
        CHARACTER*2  CRLF

        EXTERNAL     CRLF
                
C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT    (IN) :: NSRC                  ! no. sources
        INTEGER, INTENT    (IN) :: JTIME                 ! HHMMSS
        INTEGER, INTENT    (IN) :: DAYBEGT ( NSRC )      ! begin. time for day
        REAL   , INTENT    (IN) :: VALBYSRC( NSRC )      ! per-source values
        INTEGER, INTENT    (IN) :: SRCARRAY( NSRC,2 )    ! per-source FIPS and averaging type
        REAL   , INTENT    (IN) :: MINTEMP               ! minimum temperature
        REAL   , INTENT    (IN) :: MAXTEMP               ! maximum temperature
        LOGICAL, INTENT    (IN) :: SKIPDATA             ! skip data for non-day averaging
        INTEGER, INTENT(IN OUT) :: NDAYSRC ( NSRC,24 )   ! no. days to average per source
        REAL   , INTENT(IN OUT) :: HOUROUT ( NSRC,24 )   ! hourly temp per source 

C...........   Other local variables
        INTEGER     L, S        ! counters and indices
        INTEGER     TIMESLOT    ! array location

        REAL        VAL         ! tmp value

        LOGICAL, SAVE :: INITIAL = .TRUE.  ! true: first time

        CHARACTER*300 BUFFER        ! formatted source info for messages
        CHARACTER*300 MESG          ! message buffer
        CHARACTER(LEN=SRCLEN3) CSRC ! tmp concat source characteristics
        CHARACTER*16 :: PROGNAME = 'HOURTEMP' ! program name

C***********************************************************************
C   begin body of subroutine HOURTEMP

C.........  For the first time, initialize all entries to zero
        IF( INITIAL ) THEN
            HOUROUT = 0.  ! array
            INITIAL = .FALSE.
        END IF

C.........  Loop through sources
        DO S = 1, NSRC

            CSRC = CSOURC( S )
            VAL = VALBYSRC( S )

            IF( VAL > AMISS3 )THEN

C.................  Check that temperature is within min and max bounds

                IF( VAL < MINTEMP ) THEN

C.....................  Round value up to minimum
                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L )
                    WRITE( MESG, 94020 )
     &                     'Increasing hourly temperature from',
     &                     VAL, 'to', MINTEMP, 'for source' //
     &                     CRLF() // BLANK10 // BUFFER( 1:L ) // '.'
ccs                    CALL M3MESG( MESG )

                    VAL = MINTEMP 

                ELSEIF( VAL > MAXTEMP ) THEN

C.....................  Round value down to maximum
                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L )
                    WRITE( MESG, 94020 )
     &                     'Decreasing hourly temperature from',
     &                     VAL, ' to', MAXTEMP, 'for source' //
     &                     CRLF() // BLANK10 // BUFFER( 1:L ) // '.'
ccs                    CALL M3MESG( MESG )

                    VAL = MAXTEMP
                
                END IF

C.............  Store temperature value in appropriate time slot in output array

C.................  Appropriate 24 hour time will be day starting time (6 AM in local 
C                   time zone ) subtracted from met data time (in GMT)
                TIMESLOT = 1 + ( JTIME - DAYBEGT( S ) ) / 10000 
                
C.................  If timeslot is less than zero, add 24; if better data comes
C                   along, the old data will get overwritten (helps in case of
C                   one running one day)
                IF( TIMESLOT <= 0 ) THEN
                    TIMESLOT = TIMESLOT + 24
                END IF
                
                IF( SRCARRAY( S,2 ) == DAILY ) THEN
                    HOUROUT( S,TIMESLOT ) = VAL
                    NDAYSRC( S,TIMESLOT ) = 1
                ELSE
                    IF( .NOT. SKIPDATA ) THEN
                        HOUROUT( S,TIMESLOT ) =
     &                                    HOUROUT( S,TIMESLOT ) + VAL
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
