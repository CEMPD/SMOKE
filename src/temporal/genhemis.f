
        SUBROUTINE GENHEMIS( NPLE, JDATE, JTIME, TZONE, DNAME, HNAME,
     &                       NAMIN, NAMOUT, EMAC, EMACV, TMAT, EMIST )  

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
C***************************************************************************

C...........   MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC

C...........   This module contains the cross-reference tables
        USE MODXREF

C...........   This module contains the temporal profile tables
        USE MODTMPRL

C.........  This module contains emission factor tables and related
        USE MODEMFAC

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR

C...........   This module is the derived meteorology data for emission factors
        USE MODMET

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures
        INCLUDE 'M5CNST3.EXT'   !  MOBILE5a/b constants

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        LOGICAL         ENVYN
        INTEGER         FIND1
        CHARACTER*10    HHMMSS
        INTEGER         INDEX1
        LOGICAL         ISDSTIME
        INTEGER         WKDAY

        EXTERNAL        CRLF, ENVYN, FIND1, HHMMSS, INDEX1, 
     &                  ISDSTIME, WKDAY

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
        REAL        , INTENT (OUT)   :: EMACV( NSRC, NPLE ) ! work emis/actvy
        REAL        , INTENT (OUT)   :: TMAT ( NSRC, NPLE, 24 ) ! tmprl matrix
        REAL        , INTENT (OUT)   :: EMIST( NSRC, NPLE )     ! hourly emis

C...........   TMAT update variables

        INTEGER, SAVE :: NHRCALC             ! No. of entries in HRCALC
        INTEGER, SAVE :: HRCALC( 24 )        ! List of GMT hrs for calc'g TMAT
        INTEGER, SAVE :: MONTH ( 24, 0:23 )  ! time zone's month 1 ... 12
        INTEGER, SAVE :: DAYOW ( 24, 0:23 )  ! time zone's day   1 ... 7

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE, SAVE :: VIDX( : ) ! index for vehicle type

C...........   Other local variables

        INTEGER          C, H, I, J, K, L, M, S, V !  indices and counters

        INTEGER, SAVE :: TZMIN   ! minimum time zone in inventory
        INTEGER, SAVE :: TZMAX   ! maximum time zone in inventory

        INTEGER          DAY        ! tmp emissions day of week (1=monday)
        INTEGER          HCORR      ! hour correction factor
        INTEGER          HOUR       ! hour of day (1 ... 24)
        INTEGER          IOS        ! i/o status
        INTEGER          J1, J2, J3, J4 ! min/max temperature indices
        INTEGER       :: JMAXALL = 0   ! max value of JMAX
        INTEGER          JMIN, JMAX ! min and max value of j1-j4 for one iter
        INTEGER       :: JMINALL = 9999 ! min value of JMIN
        INTEGER, SAVE :: LDATE = -1 ! date used in previous subroutine call
        INTEGER          MON        ! tmp month number (1=Jan)
        INTEGER          NACT       ! activity index
        INTEGER          NETP       ! emission type index
        INTEGER          NIDX       ! tmp emission factor name index
        INTEGER       :: OSRC = 0   ! src count of sources outside the grid
        INTEGER          PIDX       ! tmp pollutant/activity index
        INTEGER          PSI        ! parameter scheme index
        INTEGER          T1, T2     ! temperature indices
        INTEGER          TDATE      ! date for computing time zones update arr 
        INTEGER          TTIME      ! time for computing time zones update arr 

        REAL             FAC( 4 )   ! tmp emission factors
        REAL             PP, QQ     ! ndiu em factor interpolation factors
        REAL             R1, R2     ! diur em factor interpolation factors
        REAL             S1, S2     ! diur em factor interpolation factors
        REAL   , SAVE :: TINCINV    ! tmpr increment inverse
        REAL             TMAX0      ! tmp max temperature for diurnal interp 
        REAL             TMIN0      ! tmp min temperature for diurnal interp
        REAL             UFAC       ! tmp units conversion factor

        LOGICAL, SAVE :: DAYLIT   = .FALSE. ! true: TZONES are in daylight time
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

        CHARACTER*1            EFTYPE    ! tmp em fac type (N=nondiur, D=diur)
        CHARACTER*300          BUFFER    ! source info buffer 
        CHARACTER*300          MESG      ! message buffer 
        CHARACTER(LEN=SRCLEN3) CSRC      ! tmp source chars string
        CHARACTER(LEN=IOVLEN3) NAMBUF    ! variable name buffer

        CHARACTER*16  :: PROGNAME = 'GENHEMIS' ! program name

C***********************************************************************
C   begin body of subroutine GENHEMIS

C.........  If current date is less than previous date, then we know that the
C           calling program has started over with a new set of pollutants.  So
C           reset the flag for first timestep.  This 
        IF( JDATE .LT. LDATE .OR. FIRSTIME ) FIRSTSTP = .TRUE.

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

C.............  Compute the inverse temperature increment
            IF( NTMPR .GT. 0 ) TINCINV = 1./ TMMINVL

C.............  For mobile sources, define the vehicle type index for all 
C               sources so it doesn't need to be looked up each time
            IF( CATEGORY .EQ. 'MOBILE' ) THEN
        	ALLOCATE( VIDX( NSRC ), STAT=IOS )
        	CALL CHECKMEM( IOS, 'VIDX', PROGNAME )

C.................  Convert mobile vehicle type values to their index values
        	DO S = 1, NSRC
                    I = INDEX1( CVTYPE( S ), MXM5VTYP, M5VTYPES )

                    IF( I .LE. 0 ) THEN
                	L = LEN_TRIM( CVTYPE( S ) )
                	MESG = 'WARNING: Inventory vehicle type '// 
     &                     CVTYPE( S )( 1:L ) //
     &                     ' not used in MOBILE model.' // CRLF()//
     &                     BLANK10 // 'No emission factors available '//
     &                     'so emissions will not be computed.'
                	CALL M3MSG2( MESG )                   
                    END IF

                    VIDX( S ) = I

        	END DO

            END IF  ! End of mobile only section

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
c NOTE: NEXTIME Routine has bug that will cause this to not work on January 1st.
c    N: NEXTTIME will not go backwards to the previous year, but will instead change
c    N: the hour on the current year and not go back to the previous year.  The
c    N: following hack works around this problem.

C...............  For the first day of the year, work around I/O API bug
                IF( (TDATE - (TDATE/1000)*1000 ) .EQ. 001 ) THEN  ! integer math
                    J = 0
                    CALL NEXTIME( TDATE, J, -250000 )
                END IF

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

C.............  Adjust sources' time zones to account for daylight time...
C.............  Subtract 1 if date is daylight time and TZONES is not already
C               converted.  Add 1 if date is standard and TZONES has been
C               converted.
C.............  FLTRDAYL is a source-array of 0s and 1s to permit sources
C               to not get daylight time conversion.
            IF( ISDSTIME( JDATE ) .AND. .NOT. DAYLIT ) THEN

                DAYLIT = .TRUE.

                TZONES = TZONES - 1 * FLTRDAYL   ! arrays

            ELSE IF( .NOT. ISDSTIME( JDATE ) .AND. DAYLIT ) THEN

                DAYLIT = .FALSE.

                TZONES = TZONES + 1 * FLTRDAYL   ! arrays

            END IF 
                
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

C.................  Set position of activity in activity list
                NACT = INDEX1( NAMIN( V ), NIACT, ACTVTY )

C.................  Set position of emission type in emission type list
                NETP = INDEX1( NAMOUT( V ), NETYPE, EMTNAM )

C.................  Set emission factor type and units conversion
                EFTYPE = EMTEFT( NETP, NACT )

C.................  If non-diurnal emission type...
                IF( EFTYPE .EQ. 'N' ) THEN

C.....................  Set position in non-diurnal emission factor list
                    NIDX = INDEX1( NAMOUT( V ), NNDI, NDINAM )

C.....................  Loop through sources and apply emission factors to
C                       hourly activity for non-diurnal emissions    
                    DO S = 1, NSRC

C.........................  Get factors for interpolating non-diurnal emis 
C                           factors and get temperature indices
C.........................  Trap temperatures against min and max available 
C                           temperatures and output messages to that affect
                        CALL TMPINTRP( S, TASRC( S ), OSRC, 
     &                                 T1, T2, PP, QQ       )

C.........................  Adjust hour index to local time zone for setting 
C                           index for this sources emission factor for this hour
                        K   = 1 + MOD( HOUR + HCORR - TZONES(S), 24 )

C.........................  Get PSI from unsorted emission factor table
                        M   = EFSIDX( S,NACT )
                        PSI = IPSIA ( M,K )

C.........................  Find index of PSI in emission factor tables using
C                           sorted list of PSIs
                        M = FIND1( PSI, NPSI( NACT ), PSILIST( 1,NACT ))

C.........................  Set vehicle type code
                        C = VIDX( S )

C.........................  Skip computation if no emission factors available 
C                           for inventory type
                        IF( C .LE. 0 ) CYCLE

C.........................  Apply emission factors tp hourly activity data
C.........................  Convert to tons (assuming EFs are in grams)
                        EMIST( S,V ) = EMIST( S,V ) * 
     &                               ( PP * EFNDIALL( T1,C,M,NIDX ) +
     &                                 QQ * EFNDIALL( T2,C,M,NIDX )   )

                    END DO  ! End loop on sources
 
C.................  If diurnal emisison type...
                ELSE IF ( EFTYPE .EQ. 'D' ) THEN

C.....................  Set position in non-diurnal emission factor list
                    NIDX = INDEX1( NAMOUT( V ), NDIU, DIUNAM )

C.....................  Loop through sources and apply emission factors and
C                       to hourly activity for diurnal emissions    
                    DO S = 1, NSRC

C.........................  Adjust hour index to local time zone for setting 
C                           index for this sources emission factor for this hour
                        K   = 1 + MOD( HOUR + HCORR - TZONES(S), 24 )

C.........................  Get PSI from unsorted emission factor table
                        M   = EFSIDX( S,NACT )
                        PSI = IPSIA ( M,K )

C.........................  Find index of PSI in emission factor tables using
C                           sorted list of PSIs
                        M = FIND1( PSI, NPSI( NACT ), PSILIST( 1,NACT ))

C.........................  Set vehicle type code
                        C = VIDX( S )

C.........................  Skip computation if no emission factors available 
C                           for inventory type
                        IF( C .LE. 0 ) CYCLE

C.........................  Get min/max temperatures indices for diurnal  
C                           emission factors
                	J1 = METIDX( S,1 )
                	J2 = METIDX( S,2 )
                	J3 = METIDX( S,3 )
                	J4 = METIDX( S,4 )

C.........................  Make sure indices are not out of bounds. A more 
C                           elaborate approach is not being taken because
C                           of performance concerns
                        JMIN    = MIN( J1, J2, J3, J4 )
                        JMAX    = MAX( J1, J2, J3, J4 )
                        JMINALL = MIN( JMINALL, JMIN )
                        JMAXALL = MAX( JMAXALL, JMAX )
                        IF( JMAX .GT. NVLDTMM ) CYCLE
    
C.........................  Initialize interpolation factor to make emis zero
                	R1 = 0.
                	R2 = 0.
                	S1 = 0.
                	S2 = 0.

C.........................  When the min/max index is non-zero, then diurnal
C                           emission factors will be expected.
                	IF( JMIN .GT. 0 ) THEN

                            FAC(1) = EFDIUALL( J1, C, M, NIDX )
                            FAC(2) = EFDIUALL( J2, C, M, NIDX )
                            FAC(3) = EFDIUALL( J3, C, M, NIDX )
                            FAC(4) = EFDIUALL( J4, C, M, NIDX )

C.............................  Check if diurnal emission factors have been 
C                               calculated for this temperature. If one  
C                               emission type has been, they all have been.
                            IF( FAC(1) .LT. 0. .OR.
     &                          FAC(2) .LT. 0. .OR.
     &                          FAC(3) .LT. 0. .OR.
     &                          FAC(4) .LT. 0.   ) THEN

                                CSRC = CSOURC( S ) 
                                CALL FMTCSRC( CSRC, NCHARS, BUFFER, L ) 
                                
                        	WRITE( MESG,94010 ) 'WARNING: '// 
     &                             'Diurnal EFs are not available ' //
     &                             'for PSI', PSI, CRLF() // BLANK10 //
     &                             'Setting emissions to zero for' //
     &                             CRLF()// BLANK10 // BUFFER( 1:L )
                        	CALL M3MESG( MESG )

C.............................  Get factors for interpolating diurnal em facs
                            ELSE
                		TMIN0 = VLDTMIN( J1 )
                		TMAX0 = VLDTMAX( J1 )

                                R1 = TINCINV*( TKMIN(S)-TMIN0 )
                		R1 = AMOD( R1, 1.0 )
                		R2 = 1.0 - R1
                		S1 = TINCINV*( TKMAX(S)-TMAX0 )
                		S1 = AMOD( S1, 1.0 )
                		S2 = 1.0 - S1

                            END IF

C.........................  Otherwise, there are no diurnal emissions because
C                           the minimum/maximum temperatures were too close
C                           together.
                	ELSE
                            
                	    FAC = 0.  ! array

                	END IF  ! End check that diurnal emissions are zero

                	EMIST( S,V ) = EMIST( S,V ) *
     &                                 ( S2 * ( R2*FAC( 1 ) +
     &                                          R1*FAC( 2 )  ) +
     &                                   S1 * ( R2*FAC( 3 ) +
     &                                          R1*FAC( 4 )  )   )

                    END DO      ! End of loop on sources
                END IF          ! End non-diurnal or diurnal check
            END IF              ! End namin != namout

        END DO                  ! End pollutant loop

C.........  Error found if J indices were out of bounds
        IF ( JMINALL .LT. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'INTERNAL ERROR: Minimum min/max temperature ' //
     &             'index from MINMAXT was < 0'
            CALL M3MSG2( MESG )

        END IF

        IF ( JMAXALL .GT. NVLDTMM ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'INTERNAL ERROR: Maximum min/max ' //
     &             'temperature index from MINMAXT was > ', NVLDTMM
            CALL M3MSG2( MESG )

        END IF

        IF ( EFLAG ) THEN
            MESG = 'Problem with meteorology indices.'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
        END IF

C.........  Reset previous date for future iterations
        LDATE = JDATE

C.........  Set controller to turn off first time step setting
        FIRSTSTP = .FALSE.

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
94030   FORMAT( 8X, 'at time ', A8 )

        END SUBROUTINE GENHEMIS

