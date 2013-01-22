
        SUBROUTINE HOURMET( NSRC, CNTYSRC, TA, QV, PRES, JDATE, JTIME, 
     &                      DAYBEGT, LDAYSAV, PDTEMP, HFLAG )

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
        USE MODMBSET, ONLY: NREFC, MCREFSORT, MCREFIDX,
     &                      NREFF, FMREFSORT, NFUELC, FMREFLIST

C...........   This module is the derived meteorology data for emission factors
        USE MODMET, ONLY: TKHOUR, RHHOUR, RHDAY,
     &                    NTKHOUR, NRHHOUR, RHTBIN, NRHTBIN, NFUEL,
     &                    FUELIDX, MINTSRC, MAXTSRC, MAXTDAY, MINTDAY

C...........   This module contains the gridding surrogates tables
        USE MODSURG, ONLY: NSRGFIPS, SRGFIPS, NCELLS, FIPCELL
     
        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS
        CHARACTER(2) CRLF
        INTEGER      ENVINT
        LOGICAL      ISDSTIME
        REAL         CALCRELHUM
        INTEGER      FIND1
        INTEGER      FIND1FIRST

        EXTERNAL     CRLF, ENVINT, ISDSTIME, CALCRELHUM, FIND1, FIND1FIRST
                
C...........   SUBROUTINE ARGUMENTS
        INTEGER,      INTENT    (IN) :: NSRC                  ! no. sources
        INTEGER,      INTENT    (IN) :: CNTYSRC( NSRC )       ! no. counties
        REAL   ,      INTENT    (IN) :: TA( * )               ! gridded temp data
        REAL   ,      INTENT    (IN) :: QV( * )               ! gridded mixing ratio data
        REAL   ,      INTENT    (IN) :: PRES( * )             ! gridded pressure data
        INTEGER,      INTENT    (IN) :: JDATE                 ! YYYYDDD
        INTEGER,      INTENT    (IN) :: JTIME                 ! HHMMSS
        INTEGER,      INTENT    (IN) :: DAYBEGT ( NSRC )      ! begin. time for day
        LOGICAL,      INTENT    (IN) :: LDAYSAV ( NSRC )      ! true: use daylight time
        INTEGER,      INTENT    (IN) :: PDTEMP                ! RPP temperature increment
        LOGICAL,      INTENT    (IN) :: HFLAG                 ! true: use specific humidity (no RH)

C...........   Other local variables
        INTEGER     C, K, L, LL, N, S, I, J, NR, NF, NT, T       ! counters and indices
        INTEGER     IOS         ! I/O status
        INTEGER     MONTH,DAY   ! processing month and date
        INTEGER     REFCOUNTY   ! ref. county FIPS code
        INTEGER     TIMESLOT    ! array location
        INTEGER     CURMONTH, NMON
        
        INTEGER, SAVE :: MXWARN ! maximum number of warnings
        INTEGER, SAVE :: NWARN  ! total number of warnings printed

        REAL        MINTMP      ! min temperature value in Farenheight
        REAL        MAXTMP      ! max temperature value in Farenheight
        REAL        TEMPVAL     ! temperature value in Farenheight
        REAL        TEMPTMP     ! tmp temperature value in Farenheight
        REAL        RHVAL       ! RH value
        REAL        RHTMP       ! tmp RH value

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

C.........  If last day of month, process monthly averages
        CALL DAYMON( JDATE, MONTH, DAY )

C.........  Loop through sources
        DO S = 1, NSRC

C.........  Apply ungridding matrix from a (possible) subgrid to data on base 
C           grid.  If no subgrid, then XOFF and YOFF will be 1 and no problem.
            LL = FIND1( CNTYSRC( S ), NSRGFIPS, SRGFIPS )

            IF( LL < 1 ) CYCLE

            IF( NCELLS( LL ) > 0 ) THEN
                RHVAL = 0.0
                TEMPVAL = 0.0
            ELSE
		RHVAL = BADVAL3
                TEMPVAL = BADVAL3
            END IF

            N = 0
            DO I = 1, NCELLS( LL )

C.................  Count no of cell used in county-level averaging
                N = N + 1

C.................  Get column and row from subgrid
                C = FIPCELL( I,LL )

C.................  Convert K to F
                TEMPTMP = 1.8 * TA( C ) - 459.67

                MAXTSRC( S ) = MAX( TEMPTMP, MAXTSRC( S ) )
                MINTSRC( S ) = MIN( TEMPTMP, MINTSRC( S ) )

                MAXTDAY( S ) = MAX( TEMPTMP, MAXTDAY( S ) )
                MINTDAY( S ) = MIN( TEMPTMP, MINTDAY( S ) )

C.................  Calculate RH using Temp. Pressure and mixing ratio values
                IF( HFLAG ) THEN
                    RHTMP = QV( C ) / ( QV( C ) + 1 )     ! Specifici Humidity
                ELSE
                    RHTMP = CALCRELHUM( TA( C ), PRES( C ), QV( C ) ) ! relative humidity
                END IF

C.................  Store RH into temperature bins
                NT = 0
                DO T = -150, 150-PDTEMP, PDTEMP
                    NT = NT + 1
                    MINTMP = REAL( T ) - ( REAL( PDTEMP ) / 2.0 )
                    MAXTMP = MINTMP + REAL( PDTEMP )

                    IF ( MINTMP < TEMPTMP .AND. TEMPTMP <= MAXTMP ) THEN

                        REFCOUNTY = MCREFSORT( S,2 )
                        NR = FIND1( REFCOUNTY,NREFC, MCREFIDX( :,1 ) )

                        L = FIND1FIRST( REFCOUNTY, NREFF, FMREFSORT( :,1 ) )
                        K = FIND1FIRST( REFCOUNTY, NFUELC,FMREFLIST( :,1 ) )
                        NMON = FMREFLIST( K, 2 )   ! no month of ref county

C.........................  Loop over months per ref. county
                        DO J = L, L + NMON - 1
                           CURMONTH  = FMREFSORT( J,3 )    ! processing  current month per ref. cou
                           IF( CURMONTH == MONTH ) NF = FMREFSORT( J,2 )  ! processing fuelmonth/co
                        END DO

                        RHTBIN ( NR,NF,NT ) = RHTBIN ( NR,NF,NT )  + RHTMP
                        NRHTBIN( NR,NF,NT ) = NRHTBIN( NR,NF,NT ) + 1

                    END IF

                END DO

                RHVAL = RHVAL + RHTMP
                TEMPVAL = TEMPVAL + TEMPTMP

            END DO

            TEMPVAL = TEMPVAL / N    ! averaged Temp by source
            RHVAL   = RHVAL   / N    ! averaged RH by srouce

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

            IF( TEMPVAL > AMISS3 .AND. RHVAL > 0.0 )THEN

C.................  Store daily RH values for daily SMOKE-ready output
                RHDAY( S,TIMESLOT ) = RHDAY( S,TIMESLOT ) + RHVAL

C.................  Store values in hourly arrays for mothly SMOKE-ready output               
                TKHOUR( S,TIMESLOT )  = TKHOUR( S,TIMESLOT ) + TEMPVAL
                NTKHOUR( S,TIMESLOT ) = NTKHOUR( S,TIMESLOT ) + 1

                RHHOUR( S,TIMESLOT )  = RHHOUR( S,TIMESLOT ) + RHVAL
                NRHHOUR( S,TIMESLOT ) = NRHHOUR( S,TIMESLOT ) + 1

            END IF

        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

94020   FORMAT( A, 4( 1X, F8.2, 1X, A ) )
 
        END SUBROUTINE HOURMET
