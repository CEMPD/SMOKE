
        SUBROUTINE GENHEMIS( CATEGORY, NSRC, NPOL, NDAYSRC, NHRSRC, 
     &                       JDATE, JTIME, TZONE, DNAME, HNAME, LDSPOL, 
     &                       LHSPOL, PNAM, ZONES, TPF, INDXD, EMISD, 
     &                       INDXH, EMISH, EMIS, EMISV, TMAT, EMIST )                       

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine is responsible for applying the temporal profiles
C      to each source using the temporal cross-reference tables.
C NOTE: Must use EMISV because want to keep EMIS constant for all of the
C       time steps.  This is important to support a different number of
C       day- or hour-specific sources per time step. NOTE: many of the IN/OUT
C       arrays are like that because of memory allocation
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by M. Houyoux 1/99
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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

C...........   MODULES for public variables
C...........   This module contains the cross-reference tables
        USE MODXREF

C...........   This module contains the temporal profile tables
        USE MODTPRO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL         ENVYN
        INTEGER         FIND1
        CHARACTER*10    HHMMSS
        LOGICAL         ISDSTIME
        CHARACTER*14    MMDDYY
        INTEGER         WKDAY

        EXTERNAL        ENVYN, FIND1, HHMMSS, ISDSTIME, MMDDYY, WKDAY

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN)     :: CATEGORY  ! name of source category
        INTEGER     , INTENT (IN)     :: NSRC      ! Number of sources
        INTEGER     , INTENT (IN)     :: NPOL      ! Number of pollutants
        INTEGER     , INTENT (IN)     :: NDAYSRC   ! Number of day-spec sources
        INTEGER     , INTENT (IN)     :: NHRSRC    ! Number of hour-spec sources
        INTEGER     , INTENT (IN)     :: JDATE     ! Julian date (YYYYDDD)
        INTEGER     , INTENT (IN)     :: JTIME     ! Time (HHMMSS)
        INTEGER     , INTENT (IN)     :: TZONE     ! Output time zone (typcly 0)
        CHARACTER(*), INTENT (IN)     :: DNAME     ! day-spec file name or NONE
        CHARACTER(*), INTENT (IN)     :: HNAME     ! hour-spec file name or NONE
        LOGICAL     , INTENT (IN)     :: LDSPOL( NPOL ) ! indicates dy-spec pol
        LOGICAL     , INTENT (IN)     :: LHSPOL( NPOL ) ! indicates hr-spec pol
        CHARACTER(*), INTENT (IN)     :: PNAM  ( NPOL ) ! inventory pol names
        INTEGER     , INTENT (IN OUT) :: ZONES( NSRC    ) ! time zones
        INTEGER     , INTENT (IN OUT) :: TPF  ( NSRC    ) ! temporal type
        INTEGER     , INTENT (IN OUT) :: INDXD( NDAYSRC ) ! src indx for dy-spec
        REAL        , INTENT (IN OUT) :: EMISD( NDAYSRC ) ! dy-spec emis
        INTEGER     , INTENT (IN OUT) :: INDXH( NHRSRC  ) ! src indx for hr-spec
        REAL        , INTENT (IN OUT) :: EMISH( NDAYSRC    ) ! hr-spec emis
        REAL        , INTENT (IN)     :: EMIS ( NSRC, NPOL ) ! inventory emis
        REAL        , INTENT (IN OUT) :: EMISV( NSRC, NPOL ) ! work emissions
        REAL        , INTENT (IN OUT) :: TMAT ( NSRC, NPOL, 24 ) ! tmprl matrix
        REAL        , INTENT (IN OUT) :: EMIST( NSRC, NPOL ) ! hourly emissions

C...........   TMAT update variables

        INTEGER, SAVE :: NHRCALC             ! No. of entries in HRCALC
        INTEGER, SAVE :: HRCALC( 24 )        ! List of hrs for calc'g TMAT
        INTEGER          MONTH ( 24, 0:23 )  ! time zone's month 1 ... 12
        INTEGER          DAYOW ( 24, 0:23 )  ! time zone's day   1 ... 7

C...........   Other local variables

        INTEGER          I, J, S, T, V !  indices and counters.

        INTEGER, SAVE :: TZMIN   ! minimum time zone in inventory
        INTEGER, SAVE :: TZMAX   ! maximum time zone in inventory

        INTEGER          DAY        ! tmp emissions day of week (1=monday)
        INTEGER          HOUR       ! hour of day (1 ... 24)
        INTEGER          IOS        ! i/o status
        INTEGER, SAVE :: LDATE = -1 ! date used in previous subroutine call
        INTEGER          NS         ! tmp no. of srcs w/day or hour-spec data
        INTEGER          MON        ! tmp month number (1=Jan)
        INTEGER          TDATE      ! date for computing time zones update arr 
        INTEGER          TTIME      ! time for computing time zones update arr 

        LOGICAL, SAVE :: DAYLIT   = .FALSE. ! true: ZONES are in daylight time
        LOGICAL, SAVE :: DFLAG              ! true: day-specific data
        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first call to subrtn
        LOGICAL, SAVE :: FIRSTSTP = .TRUE.  ! true: first time step
        LOGICAL, SAVE :: HFLAG              ! true: hour-specific data
        LOGICAL, SAVE :: OUTMSG             ! true: outputs date/time msg needed
        LOGICAL, SAVE :: TMATCALC           ! true: need to calculate new TMAT
        LOGICAL, SAVE :: ZONE4WM        !  True: src zone for week/mon temp prof

        CHARACTER*300          MESG      !  message buffer 
        CHARACTER(LEN=IOVLEN3) NAMBUF    !  variable name buffer

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

            ENDIF  ! End category selection

C.............  Retrieve environment variables
            MESG = 'Assign weekly/monthly profiles using time zones'
            ZONE4WM = ENVYN ( 'ZONE4WM', MESG, .TRUE., IOS )

C.............  Define the minimum and maximum time zones in the inventory
            TZMIN = MINVAL( ZONES )
            TZMAX = MAXVAL( ZONES )

C.............  Adjust TZMIN for possibility of daylight savings
            TZMIN = TZMIN - 1

C.............  Determine hours of output day for updating TMAT
            NHRCALC = TZMAX - TZMIN + 1
        
            DO I = TZMIN, TZMAX
                J = J + 1    
                HRCALC( J ) = MOD( I - TZONE + 25, 24 )
            ENDDO

C.............  Set flags for daily and hourly data
            DFLAG = ( DNAME .NE. 'NONE' )
            HFLAG = ( HNAME .NE. 'NONE' )

            FIRSTIME = .FALSE.

        ENDIF

C.........  For new date...
        IF( JDATE .NE. LDATE ) THEN

C.............  Store month and day of week for this output date for
C.............  all source's time zones and hours of day
C.............  Use same loop for case where this feature is turned off
            DO I = TZMIN, TZMAX

                TDATE = JDATE
                TTIME = 0     ! Set for loop below from 1 to 24
                CALL NEXTIME( TDATE, TTIME, ( TZONE - I ) * 10000 )

                DO T = 1, 24

                    IF( ZONE4WM ) THEN
                        CALL DAYMON( TDATE, MON, DAY ) ! month & scratch date
                        DAY = WKDAY( TDATE )           ! get day-of-week
                    ELSE
                        CALL DAYMON( JDATE, MON, DAY )
                        DAY = WKDAY( JDATE )
                    ENDIF

                    MONTH( T,I ) = MON
                    DAYOW( T,I ) = DAY

                    CALL NEXTIME( TDATE, TTIME, 10000 )

                ENDDO

            ENDDO

C.............  Set day of week based on output day
            DAY = WKDAY( JDATE )

C.............  Turn on message for day of week and date
            OUTMSG = .TRUE.

C.............  Initialize day-specific emissions array with average emissions
C               Only need to do this for a new day because the sources with 
C               day-specific data might change for each day.
            EMISV = EMIS  ! array

C.............  If day-specific emissions, prepare day-corrections:
                
            IF( DFLAG ) THEN

C.................  Read source index for this day
                IF ( .NOT. READ3( DNAME, 'INDXD', ALLAYS3,
     &                            JDATE, JTIME, INDXD      ) ) THEN

                    MESG= 'Could not read "INDXD" from file "' // 
     &                    DNAME //'".'
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                END IF      !  if read3() failed on dname

C.................  Find location of last non-zero entry in index
                NS = FIND1( 0, NDAYSRC, INDXD ) - 1

C.................  Read emissions for day-specific pollutants and store
C                   data in EMISV using index INDXD
                DO V = 1, NPOL

                    IF( LDSPOL( V ) ) THEN

                        NAMBUF = PNAM( V )
                        IF( .NOT. READ3( DNAME, NAMBUF, ALLAYS3,
     &                                   JDATE, JTIME, EMISD     )) THEN

                            MESG = 'Could not read "' // 
     &                             NAMBUF( 1:LEN_TRIM( NAMBUF ) ) // 
     &                             '" from file "' // DNAME // '".'
                            CALL M3EXIT( PROGNAME,JDATE,JTIME,MESG,2 )

                        ENDIF      !  if read3() failed on day-specific data

                        DO I = 1, NS
                            S = INDXD( I )
                            EMISV( S,V ) = EMISD( I )
                        ENDDO

                    END IF  ! if this pollutant is day-specific
                ENDDO       ! end loop on day-specific pollutants

            END IF          ! if using day-specific emissions

C.............  Reset previous date for future iterations
            LDATE = JDATE

        ENDIF   ! End new date check

C.........  Set integer hour of day for output time 

        HOUR = 1 + MOD( JTIME / 10000 , 24 )

C.........  Determine if this hour is one that requires a new TMAT (because one
C           or more sources are starting a new day)
C.........  Allow for case where ZONE4WM is false, then only need to update TMAT
C           once per day
        IF( ZONE4WM ) THEN
            I        =  FIND1( HOUR, NHRCALC , HRCALC )
            TMATCALC = ( I .GT. 0 .OR. FIRSTSTP )
 
        ELSE
            TMATCALC = ( JDATE .NE. LDATE  )
 
        ENDIF

C.........  Construct TMAT -- array of effective composite profile coefficients

        IF( TMATCALC ) THEN

C.............  Adjust sources' time zones to account for daylight time...
C.............  Subtract 1 if date is daylight and ZONES is not and
C               add 1 if date is not daylight and ZONES is daylight
            IF( ISDSTIME( JDATE ) .AND. .NOT. DAYLIT ) THEN

                DAYLIT = .TRUE.

                ZONES = ZONES - 1  ! array

            ELSEIF( .NOT. ISDSTIME( JDATE ) .AND. DAYLIT ) THEN

                DAYLIT = .FALSE.

                ZONES = ZONES + 1  ! array

            ENDIF 
                
C.............  If there are no weekend packets...
            IF ( DAY .LT. 6 .OR. NEND .EQ. 0 ) THEN  !  no weekend packet

                IF( OUTMSG ) THEN
                    CALL M3MSG2( 'Processing ' //
     &                           DAYS( DAY ) // MMDDYY( JDATE ) )
                    OUTMSG = .FALSE.
                ENDIF

                IF( CATEGORY .NE. 'MOBILE' ) THEN

                    CALL MKTMAT( NSRC, NPOL, JDATE, TZONE, NDIU, 
     &                           NDAYSRC, INDXD, ZONES, TPF,
     &                           MDEX, WDEX, DDEX, MONTH, DAYOW,
     &                           LDSPOL, DIUFAC, TMAT )

                ELSE

C                   CALL MKVMAT( )

                END IF ! End category selection

C.............  If there are weekend packets, and it's a weekend...
            ELSE IF( DAY .GE. 6 .AND. NEND .GT. 0 ) THEN

                IF( OUTMSG ) THEN
                    CALL M3MSG2( 'Processing weekend day ' //
     &                           DAYS( DAY ) // MMDDYY( JDATE ) )
                    OUTMSG = .FALSE.
                ENDIF

                IF( CATEGORY .NE. 'MOBILE' ) THEN

                    CALL MKTMAT( NSRC, NPOL, JDATE, TZONE, NEND, 
     &                           NDAYSRC, INDXD, ZONES, TPF,
     &                           MDEX, WDEX, EDEX, MONTH, DAYOW,
     &                           LDSPOL, ENDFAC, TMAT )
                ELSE

C                   CALL MKVMAT( )

                END IF ! End category selection

            ENDIF      ! End day selection

        END IF      ! if jdate not ldate

C.........  Write to screen because WRITE3 only writes to LDEV
        WRITE( *, 94030 ) HHMMSS( JTIME )

C.........  Apply diurnal factors to all sources and pollutants. Note that
C           for mobile emissions, NSRC = NMSRC * NVTYPE
        IF( CATEGORY .NE. 'MOBILE' ) THEN
            DO V = 1, NPOL
                DO S = 1, NSRC

                    EMIST( S,V ) = EMISV( S,V ) * TMAT( S,V,HOUR )

                END DO
            END DO

        ELSE

C           Insert later for mobile, it is different becuse also must multiply
C           in interpolated emission factors, unless I've already done that. 
C           May be possible to not need an IF statement here if plan carefully.

        ENDIF

C.........  If hour-specific emissions, overwrite hourly emissions with 
C           hour-specific data for pollutants in the hour-specific input file
        IF( HFLAG ) THEN

C.............  Read source index for this hour
            IF( .NOT. READ3( HNAME, 'INDXH', ALLAYS3,
     &                        JDATE, JTIME, INDXH      ) ) THEN

                MESG= 'Could not read "INDXH" from file "' // 
     &                HNAME //'".'
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

            END IF      !  if read3() failed on dname

C.............  Find location of last non-zero entry in index
            NS = FIND1( 0, NHRSRC, INDXH ) - 1

C.............  Read all variables from hour-specific file because
C.............  file structure is IDDATA3

            DO V = 1, NPOL

                IF( LHSPOL( V ) ) THEN

                    NAMBUF = PNAM( V )
                    IF( .NOT. READ3( HNAME, NAMBUF, ALLAYS3,
     &                               JDATE, JTIME, EMISH     ) ) THEN

                        MESG = 'Could not read "' // 
     &                         NAMBUF( 1:LEN_TRIM( NAMBUF ) ) // 
     &                         '" from file "' // HNAME //'".'
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                    ENDIF      !  if read3() failed on hourly data

                    DO I = 1, NS
                        S = INDXH( I )
                        EMIST( S,V ) = EMISH( I )
                    END DO

                END IF  ! if this pollutant is hour-specific
            END DO      ! end loop on hour-specific pollutants
                   
        ENDIF           ! if using hour-specific emissions 

C.........  Set controller to turn off first time step setting
        FIRSTSTP = .FALSE.

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
94030   FORMAT( 8X, 'at time ', A8 )

        END SUBROUTINE GENHEMIS

