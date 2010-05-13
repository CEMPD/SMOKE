
        SUBROUTINE WRAVGMET( NSRC, AVGTYPE, ODEV, SDATE ) 

C***********************************************************************
C  subroutine body starts at line 78
C
C  DESCRIPTION:
C       Averages hourly meteorology data based on number of sources
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     2/10: Created by B.H. Baek
C
C***********************************************************************
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
C***********************************************************************

C...........   MODULES for public variables
C.........  This module is used for MOBILE6 setup information                
        USE MODMBSET, ONLY: NINVC, NREFC, MCREFSORT, MCREFIDX,
     &                      NREFF, FMREFSORT, NFUELC, FMREFLIST

 
C...........   This module is the derived meteorology data for emission factors
        USE MODMET, ONLY: TKHOUR, RHHOUR, MAXTSRC, MINTSRC, MINTEMP,
     &                    MAXTEMP, RHFUEL, MAXTFUEL, MINTFUEL, NFUEL,
     &                    FUELIDX, NDAYSRC

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS 
        CHARACTER(2)    CRLF
        INTEGER         FIND1
        INTEGER         FIND1FIRST
        CHARACTER(14)   MMDDYY

        EXTERNAL     CRLF, FIND1, FIND1FIRST, MMDDYY

C...........   SUBROUTINE ARGUMENTS
        INTEGER,     INTENT   (IN) :: NSRC      ! no. sources
        CHARACTER(*),INTENT   (IN) :: AVGTYPE   ! temporal resolution
        INTEGER,     INTENT   (IN) :: ODEV      ! SMOKE ready output file
c        INTEGER,     INTENT   (IN) :: ODEV2     ! MOVES ready output file 
        INTEGER,     INTENT   (IN) :: SDATE     ! outpout date  
c        INTEGER,     INTENT   (IN) :: PPTEMP    ! rateperprofile temp increment
c        LOGICAL,     INTENT   (IN) :: COMPLETE  ! true: fuelmonth is already processed for refcounty

C.........  Local array
        REAL   , ALLOCATABLE :: TKPRO( : ) !  24hr temp. profile
        REAL   , ALLOCATABLE :: TMPRO( : ) !  tmp 24hr temp. profile
        REAL   , ALLOCATABLE :: TKREFHR( : ) !  24hr temp. profile

C...........   Other local variables
        INTEGER I, IT, J, K, L, LT, N, NR, MX, MN, R, S, T           ! counters and indices                     
        
        INTEGER IOS                       ! I/O status
        INTEGER IC                        ! number of inven cnty per ref cnty
        INTEGER REFCOUNTY                 ! ref. county FIPS code
        INTEGER INVCOUNTY                 ! inv. county FIPS code
        INTEGER PRCOUNTY                  ! previous ref. county
        INTEGER IMAXT, IMINT              ! current integer min/max temp.
        INTEGER MAXT , MINT               ! tmp remainder of max/min
        INTEGER MONTH, DAY                ! processing month and day
        INTEGER NMON                      ! no of month per ref. county
        INTEGER FUELMONTH                 ! current fuelmonth for ref. county
        INTEGER CURMONTH                  ! current month
        INTEGER NF, PRNF                  ! current and prvious county fuelmonths

        REAL    RHSUM                     ! sum of RH
        REAL    RHAVG                     ! avg of RH
        REAL    RHREFSUM                  ! sum of RH for ref county
        REAL    RHREFAVG                  ! avg of RH for ref county
        REAL    MAXTREF, MINTREF          ! ref county min/max temperatures
        REAL    TKMAX                     ! tmp max temperatures
        REAL    TKMIN                     ! tmp min temperatures
        REAL    TKDIF, DT                 ! tmp DIFF of min/max temperatures
        REAL    TKMED                     ! tmp median temperatures
        
        CHARACTER(1)  :: HD = 'M'           ! temporal resolution header
        CHARACTER(32)    TPROID             ! temporal resolution header
        CHARACTER(300)   MESG               ! message buffer
        CHARACTER(16) :: PROGNAME = 'WRAVGMET' ! program name

C***********************************************************************
C   begin body of subroutine WRAVGMET

C.........  Allocate local arrays
        ALLOCATE( TKPRO( 24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TKPRO', PROGNAME )
        ALLOCATE( TMPRO( 24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPRO', PROGNAME )
        ALLOCATE( TKREFHR( 24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TKREFHR', PROGNAME )
        TKPRO = 0.0
        TMPRO = 0.0
        TKREFHR = 0.0

C.........  Define temporal resolution header for temperatureProfileID
        IF( AVGTYPE == 'DAILY' ) HD = 'D'
        IF( AVGTYPE == 'MONTHLY' ) HD = 'M'
        IF( AVGTYPE == 'EPISODE' ) HD = 'E'
        
        CALL DAYMON( SDATE, MONTH, DAY )
        
C.........  Loop through all counties
        IC = 0
        PRCOUNTY = 0
        DO S = 1, NSRC

            INVCOUNTY = MCREFSORT( S,1 )
            REFCOUNTY = MCREFSORT( S,2 )

            L = 0
            N = 0
            RHSUM = 0.0
            DO T = 1, 24
                IF( RHHOUR( S,T ) > 0 ) N = N + 1
                RHSUM = RHSUM + RHHOUR( S,T )
c           print*,S,T,RHHOUR(S,T),'S,T,RHHOUR'
            END DO

C.............  Calculation hour profiles AND avg RH based on daily total temp.
c               Convert K to F
            RHAVG  = RHSUM / N
            MAXTEMP= MAXTSRC( S )
            MINTEMP= MINTSRC( S )
            
            L = FIND1FIRST( REFCOUNTY, NREFF, FMREFSORT( :,1 ) )
            K = FIND1FIRST( REFCOUNTY, NFUELC,FMREFLIST( :,1 ) )
            NMON = FMREFLIST( K, 2 )   ! no month of ref county

C.............  Loop over months per ref. county
            FUELMONTH = 0
            DO J = L, L + NMON - 1

                CURMONTH  = FMREFSORT( J,3 )    ! processing  current month per ref. county

C.................  Skip other months
                IF( CURMONTH  == MONTH ) THEN
                    FUELMONTH = FMREFSORT( J,2 )    ! processing fuelmonth/county
                END IF

            END DO

c         print*,INVCOUNTY,month,fuelmonth,RHAVG,MINTEMP,MAXTEMP

C.............  write inventory county min/max and avg RH 
            WRITE( ODEV,94040 ) INVCOUNTY, FUELMONTH, MONTH, RHAVG,
     &                           MINTEMP, MAXTEMP
     
        END DO

C.........  Re-initialize arrays
        NDAYSRC = 0
        TKHOUR = 0.0
        RHHOUR = 0.0
        MAXTSRC = BADVAL3
        MINTSRC = -1*BADVAL3
       
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I6, :, 1X ) )

94020   FORMAT( A, 4( 1X, F8.2, 1X, A ) )

94030   FORMAT( A, I7.7, I3.3 )   

94040   FORMAT( I6.6, I5, 3X, I5, 3F10.2 )   
 
 
        END SUBROUTINE WRAVGMET
