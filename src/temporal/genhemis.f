
        SUBROUTINE GENHEMIS( NPLE, JDATE, JTIME, TZONE, DNAME, HNAME,
     &                       NAMIN, NAMOUT, EMAC, EMFAC, EMACV, TMAT, 
     &                       EMIST )  

C***********************************************************************
C  subroutine body starts at line 173
C
C  DESCRIPTION:
C      This subroutine is responsible for applying the temporal profiles
C      to each source using the temporal cross-reference tables.
C      NOTE - Must use EMISV because want to keep EMIS constant for all of the
C      time steps.  This is important to support a different number of
C      day- or hour-specific sources per time step. Many of the IN/OUT
C      arrays are like that because of memory allocation
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by M. Houyoux 1/99
C
C************************************************************************
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

C...........   MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC, ONLY: TZONES, TPFLAG

C...........   This module contains the cross-reference tables
        USE MODXREF, ONLY: MDEX, WDEX, DDEX

C...........   This module contains the temporal profile tables
        USE MODTMPRL, ONLY: NHOLIDAY, HOLJDATE, HOLALTDY, HRLFAC

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: INDXD, INDXH, EMACD, EMACH, NDYSRC, NHRSRC,
     &                      LDSPOA, LHSPOA, LHPROF

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC, CATEGORY, NIPPA, EAREAD, EACNV

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)   CRLF
        LOGICAL         ENVYN
        INTEGER         FIND1
        CHARACTER(10)   HHMMSS
        INTEGER         INDEX1
        INTEGER         WKDAY

        EXTERNAL        CRLF, ENVYN, FIND1, HHMMSS, INDEX1, WKDAY

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN)    :: NPLE      ! Number of pols+emis types
        INTEGER     , INTENT (IN)    :: JDATE  ! Julian date (YYYYDDD) in TZONE
        INTEGER     , INTENT (IN)    :: JTIME     ! Time (HHMMSS) in TZONE
        INTEGER     , INTENT (IN)    :: TZONE     ! Output time zone (typcly 0)
        CHARACTER(*), INTENT (IN)    :: DNAME     ! day-spec file name or NONE
        CHARACTER(*), INTENT (IN)    :: HNAME     ! hour-spec file name or NONE
        CHARACTER(*), INTENT (IN)    :: NAMIN ( NPLE )      ! inv pol names
        CHARACTER(*), INTENT (IN)    :: NAMOUT( NPLE )      ! inv pol names
        REAL        , INTENT (IN)    :: EMAC ( NSRC, NPLE ) ! inv emis or actvty
        REAL        , INTENT (IN)    :: EMFAC( NSRC, NPLE ) ! emission factors
        REAL        , INTENT (OUT)   :: EMACV( NSRC, NPLE ) ! work emis/actvy
        REAL        , INTENT (OUT)   :: TMAT ( NSRC, NPLE, 24 ) ! tmprl matrix
        REAL        , INTENT (OUT)   :: EMIST( NSRC, NPLE )     ! hourly emis

C...........   TMAT update variables

        INTEGER, SAVE :: NHRCALC             ! No. of entries in HRCALC
        INTEGER, SAVE :: HRCALC( 24 )        ! List of GMT hrs for calc'g TMAT
        INTEGER, SAVE :: MONTH ( 24, 0:23 )  ! time zone's month 1 ... 12
        INTEGER, SAVE :: DAYOW ( 24, 0:23 )  ! time zone's day   1 ... 7

C...........   Other local variables

        INTEGER          C, H, I, J, K, L, M, S, V !  indices and counters

        INTEGER, SAVE :: TZMIN   ! minimum time zone in inventory
        INTEGER, SAVE :: TZMAX   ! maximum time zone in inventory

        INTEGER          DAY        ! tmp emissions day of week (1=monday)
        INTEGER          HCORR      ! hour correction factor
        INTEGER          HOUR       ! hour of day (1 ... 24)
        INTEGER          IOS        ! i/o status
        INTEGER, SAVE :: LDATE = -1 ! date used in previous subroutine call
        INTEGER, SAVE :: LTIME = -1 ! time used in previous subroutine call
        INTEGER          MON        ! tmp month number (1=Jan)
        INTEGER          PIDX       ! tmp pollutant/activity index
        INTEGER          TDATE      ! date for computing time zones update arr 
        INTEGER          TTIME      ! time for computing time zones update arr 

        REAL             UFAC       ! tmp units conversion factor

        LOGICAL, SAVE :: DFLAG              ! true: day-specific data
        LOGICAL, SAVE :: EFLAG    = .FALSE. ! true: error found
        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first call to subrtn
        LOGICAL, SAVE :: FIRSTSTP = .TRUE.  ! true: first time step
        LOGICAL, SAVE :: HFLAG              ! true: hour-specific data
        LOGICAL, SAVE :: OUTMSG = .TRUE.    ! true: output message for new day
        LOGICAL       :: RDFLAG = .TRUE.    ! true: read dy data for this iter
        LOGICAL       :: RHFLAG = .TRUE.    ! true: read hr data for this iter
        LOGICAL, SAVE :: TMATCALC           ! true: need to calculate new TMAT
        LOGICAL, SAVE :: UFLAG  = .FALSE.   ! true: use src-spec hr profiles
        LOGICAL, SAVE :: WKEMSG = .FALSE.   ! true: wkend-profile msg written
        LOGICAL, SAVE :: ZONE4WM        !  True: src zone for week/mon temp prof

        CHARACTER(300)     BUFFER    ! source info buffer 
        CHARACTER(300)     MESG      ! message buffer 
        CHARACTER(SRCLEN3) CSRC      ! tmp source chars string
        CHARACTER(IOVLEN3) NAMBUF    ! variable name buffer

        CHARACTER(16) :: PROGNAME = 'GENHEMIS' ! program name

C***********************************************************************
C   begin body of subroutine GENHEMIS

C.........  If current date is less than previous date or if the date is the
C           same but the time is earlier, then we know that the calling
C           program has started over with a new set of pollutants.  So
C           reset the flag for first timestep.
        IF( JDATE .LT. LDATE .OR. FIRSTIME .OR.
     &    ( JDATE .EQ. LDATE .AND. JTIME .LT. LTIME )) FIRSTSTP = .TRUE.

C.........  For the first time the subroutine is called, 
        IF( FIRSTIME ) THEN

C.............  Check source category name
            IF( CATEGORY .NE. 'AREA'   .AND.
     &          CATEGORY .NE. 'MOBILE' .AND. 
     &          CATEGORY .NE. 'POINT'        ) THEN

                MESG = 'INTERNAL ERROR: Source category ' // CATEGORY //
     &                 ' is not recognized by program ' // PROGNAME
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

            END IF  ! End category selection

C.............  Retrieve environment variables
            MESG = 'Assign weekly/monthly profiles using time zones'
            ZONE4WM = ENVYN ( 'ZONE4WM', MESG, .TRUE., IOS )

C.............  Define the minimum and maximum time zones in the inventory
            TZMIN = MINVAL( TZONES )
            TZMAX = MAXVAL( TZONES )

C.............  Adjust TZMIN for possibility of daylight savings
            TZMIN = MAX( TZMIN - 1, 0 )

C.............  Determine hours of output day in GMT for updating TMAT
            NHRCALC = TZMAX - TZMIN + 1

            K = 0        
            DO I = TZMIN, TZMAX
                K = K + 1    
                HRCALC( K ) = MOD( I - TZONE + 25, 24 )
            END DO

C.............  Set flags for daily and hourly data
            DFLAG = ( DNAME .NE. 'NONE' )
            HFLAG = ( HNAME .NE. 'NONE' )

            FIRSTIME = .FALSE.

        END IF  ! End of first time section

C.........  For new date...
        IF( JDATE .NE. LDATE ) THEN

C.............  Store month and day of week for this output date for
C.............  all source's time zones and for the *next* 24 hours in the
C               simulation
C.............  Use same loop for case where this feature is turned off
            DO I = TZMIN, TZMAX

                TDATE = JDATE
                TTIME = 0     ! Set for loop below from 1 to 24

C.................  Adjust time zone I based on output time zone and correct
C                   for negative value.
                K = I - TZONE
c                IF( K .LT. 0 ) K = 24 + K  !     (e.g., K= -1 -> K= 23)

C.................  Convert output time to local time of I time zone, adjusted
C                   by output time zone.
                CALL NEXTIME( TDATE, TTIME, -K * 10000 )

                DO H = 1, 24

C.....................  When using time zones to set monthly and weekly profiles
C.....................  NOTE - this is more correct
                    IF( ZONE4WM ) THEN

                        CALL DAYMON( TDATE, MON, DAY ) ! month & scratch day
                        DAY = WKDAY( TDATE )           ! get day-of-week
c                       DAYADJ = EMWKDAY( TDATE )

C.....................  When not using time zones as in previous emissions
C                       systems
                    ELSE
                        CALL DAYMON( JDATE, MON, DAY )
                        DAY = WKDAY( JDATE )
                        TDATE = JDATE                  ! set for holiday check
c                       DAYADJ = EMWKDAY( TDATE )
                    END IF

C.....................  Check if the date is a holiday.  If so, reset the
C                       day based on the holiday arrays settings.
C.....................  NOTE - this approach will not work for region-specific
C                       setting of holidays.
                    J = FIND1( TDATE, NHOLIDAY, HOLJDATE )
                    IF( J .GT. 0 ) DAY = HOLALTDY( J )

                    MONTH( H,I ) = MON
                    DAYOW( H,I ) = DAY

                    CALL NEXTIME( TDATE, TTIME, 10000 )

                END DO

            END DO

C.............  Set day of week based on output day
            DAY = WKDAY( JDATE )
c           DAYADJ = EMWKDAY( JDATE )
cSTOPPED HERE:
c note: Need to create EMWKDAY and update WRDAYMSG also to use it.

C.............  Write message for day of week and date
            CALL WRDAYMSG( JDATE, MESG )

        END IF   ! End new date check

C.........  Initialize day-specific emissions array with average emissions
C           Only need to do this for a new day because the sources with 
C           day-specific data might change for each day.
        EMACV = EMAC  ! array

C.............  If day-specific emissions, prepare day-corrections.
C.............  Read day-specific data for each hour, because different time
C               zones may be used for different sources and this approach
C               is much more workable.
        IF( DFLAG ) THEN

            RDFLAG = .TRUE.
C.................  Read source index for this day
            IF ( .NOT. READ3( DNAME, 'INDXD', ALLAYS3,
     &                        JDATE, JTIME, INDXD      ) ) THEN

                WRITE( MESG,94010 ) 'WARNING: Could not read "INDXD" '//
     &               'from file "'// DNAME//'", at', JDATE, ':', JTIME
                CALL M3MESG( MESG )

                RDFLAG = .FALSE.
                INDXD = 0   ! array

            END IF      !  if read3() failed on dname

        END IF          ! if using day-specific emissions

C.........  Set integer hour of day for output time 

        HOUR = 1 + MOD( JTIME / 10000 , 24 )

C.........  Determine if this TMAT needs to be updated 
        TMATCALC = ( JDATE .NE. LDATE .OR. FIRSTSTP )

C.........  Construct TMAT -- array of effective composite profile coefficients

        IF( TMATCALC ) THEN
                
C.............  Build temporal allocation matrix
            CALL MKTMAT( NSRC, NPLE, JDATE, TZONE, TZONES, 
     &                   TPFLAG, MDEX, WDEX, DDEX,
     &                   MONTH, DAYOW, TMAT )

        END IF         ! if TMAT is to be calculated

C.........  Write to screen because WRITE3 only writes to LDEV
        WRITE( *, 94030 ) HHMMSS( JTIME )

C.........  If hour-specific emissions, profiles, or activity data...
        IF( HFLAG ) THEN

            RHFLAG = .TRUE.
C.............  Read source index for this hour
            IF( .NOT. READ3( HNAME, 'INDXH', ALLAYS3,
     &                       JDATE, JTIME, INDXH      ) ) THEN

                WRITE( MESG,94010 ) 'WARNING: Could not read "INDXH" '//
     &               'from file "'// HNAME//'", at', JDATE, ':', JTIME
                CALL M3MESG( MESG )

                RHFLAG = .FALSE.
                INDXH = 0   ! array

            END IF      !  if read3() failed on HNAME

        END IF  !  End if hour-specific data

C.........  Precompute hour/zone correction factor
        HCORR = TZONE + 23

C.........  Loop through the emissions and emission types and compute the
C           hourly emissions, depending on if the input data is a pollutant
C           or an activity
        DO V = 1, NPLE

            NAMBUF = NAMIN( V )

C.............  Skip blanks that can occur when NGRP > 1
            IF ( NAMBUF .EQ. ' ' ) CYCLE

C.............  Find pollutant/activity in list of all.  Use EAREAD b/c
C               EANAM has been update to contain emission types.
C.............  NOTE - this is sloppy because NIPPA has a larger dimension
C               than EAREAD for emission types
            PIDX = INDEX1( NAMBUF, NIPPA, EAREAD )

C.............  Set units conversion factor for this pollutant/activity
            UFAC = EACNV( PIDX )

C.............  Read hourly emissions/profile/activity if current pollutant or 
C               emission type has them
            IF( LHSPOA( V ) .AND. RHFLAG ) THEN

                IF( .NOT. READ3( HNAME, NAMBUF, ALLAYS3,
     &                           JDATE, JTIME, EMACH     ) ) THEN

                    MESG = 'Could not read "' // 
     &                     NAMBUF( 1:LEN_TRIM( NAMBUF ) ) // 
     &                     '" from file "' // HNAME //'".'
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                END IF      !  if read3() failed on hourly data

            END IF  ! if this pollutant is hour-specific

C.............  If source-specific profiles are being used, update temporal 
C               matrix for the current hour
            IF( LHPROF( V ) ) THEN

                CALL UPDTMAT( NSRC, NPLE, JDATE, TZONE, V, HOUR, MONTH, 
     &                        DAYOW, TMAT )

            END IF

C.............  For all pollutants and activities...

C.............  Apply hourly factors to all sources for current pollutant or
C               activity. Also apply units conversion.
            DO S = 1, NSRC
                EMIST( S,V ) = UFAC * EMACV( S,V ) * TMAT( S,V,HOUR )
            END DO

C.............  If day-specific data are available for current pollutant
            IF( LDSPOA( V ) .AND. RDFLAG ) THEN

C.................  Read day-specific data
                NAMBUF = NAMIN( V )
                IF( .NOT. READ3( DNAME, NAMBUF, ALLAYS3,
     &                           JDATE, JTIME, EMACD    )) THEN

                    MESG = 'Could not read "' // 
     &                     NAMBUF( 1:LEN_TRIM( NAMBUF ) ) // 
     &                     '" from file "' // DNAME // '".'
                    CALL M3EXIT( PROGNAME,JDATE,JTIME,MESG,2 )

                END IF      !  if read3() failed on day-specific data

C.................  Loop through day-specific sources
                DO I = 1, NDYSRC

                    S = INDXD( I )                      ! Get source index
                    IF( S .EQ. 0 ) CYCLE                ! If no source, skip
                    IF( EMACD( I ) .LE. AMISS3 ) CYCLE  ! No day-specific emis

C.....................  Override annual adjusted emissions with day-specific 
C                       emissions and hourly profile adjustments
                    L   = DDEX( S,V )
                    DAY = DAYOW( HOUR, TZONES( S ) )                        
                    K   = 1 + MOD( HOUR + HCORR - TZONES( S ), 24 )

                    EMIST( S,V ) = UFAC * EMACD( I ) * HRLFAC( K,L,DAY )

                END DO

            END IF  ! if this pollutant is day-specific

C.............  If hourly data are available for current pollutant/activity and 
C               the values are emissions (not profiles), then overwrite with
C               this data 
            IF( RHFLAG .AND. LHSPOA( V ) .AND. .NOT. LHPROF( V ) ) THEN
                DO I = 1, NHRSRC
                    S = INDXH( I )
                    IF( S .EQ. 0 ) CYCLE
                    IF( EMACH( I ) .GT. AMISS3 ) 
     &                  EMIST( S,V ) = EMACH( I )
                END DO

            END IF

C.............  If input data is in an activity (and output an emission type)...
C               E.G. this is for mobile sources
            IF( NAMIN( V ) .NE. NAMOUT( V ) ) THEN

C.................  Loop through sources and apply emission factors to
C                   hourly activity for non-diurnal emissions    
                DO S = 1, NSRC

C.....................  Apply emission factors tp hourly activity data
C.....................  Convert to tons (assuming EFs are in grams)
                    EMIST( S,V ) = EMIST( S,V ) * EMFAC( S,V )

                END DO  ! End loop on sources
 
            END IF              ! End namin != namout

        END DO                  ! End pollutant loop

C.........  Reset previous date and time for future iterations
        LDATE = JDATE
        LTIME = JTIME

C.........  Set controller to turn off first time step setting
        FIRSTSTP = .FALSE.

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
94030   FORMAT( 8X, 'at time ', A8 )

        END SUBROUTINE GENHEMIS

