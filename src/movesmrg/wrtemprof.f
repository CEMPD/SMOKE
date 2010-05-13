
        SUBROUTINE WRTEMPROF( ODEV, MDATE, HDR, COUNTY, PMONTH, 
     &                        PPTEMP, THOUR ) 

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
     &                    FUELIDX

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
C.............  Subroutine arguments
        INTEGER  , INTENT(IN)  ::  ODEV
        INTEGER  , INTENT(IN)  ::  MDATE
        CHARACTER, INTENT(IN)  ::  HDR
        INTEGER  , INTENT(IN)  ::  COUNTY
        INTEGER  , INTENT(IN)  ::  PMONTH
        INTEGER  , INTENT(IN)  ::  PPTEMP
        REAL     , INTENT(IN)  ::  THOUR( 24 )

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
        INTEGER NF                        ! current county fuelmonths

        REAL    RHSUM                     ! sum of RH
        REAL    RHAVG                     ! avg of RH
        REAL    RHREFSUM                  ! sum of RH for ref county
        REAL    RHREFAVG                  ! avg of RH for ref county
        REAL    MAXTREF, MINTREF          ! ref county min/max temperatures
        REAL    TKMAX                     ! tmp max temperatures
        REAL    TKMIN                     ! tmp min temperatures
        REAL    TKDIF, DT                 ! tmp DIFF of min/max temperatures
        REAL    TKMED                     ! tmp median temperatures

        LOGICAL, SAVE :: FIRSTIME = .TRUE.

        CHARACTER(32)    TPROID             ! temporal resolution header
        CHARACTER(300)   MESG               ! message buffer
        CHARACTER(16) :: PROGNAME = 'WRTEMPROF' ! program name

C***********************************************************************
C   begin body of subroutine WRTEMPROF

c        print*, ODEV, MDATE, HDR, COUNTY, PMONTH, PPTEMP, 'BH11' 
C.........  Allocate local arrays
        IF( FIRSTIME ) THEN
            ALLOCATE( TKPRO( 24 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TKPRO', PROGNAME )
            ALLOCATE( TMPRO( 24 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TMPRO', PROGNAME )
            ALLOCATE( TKREFHR( 24 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TKREFHR', PROGNAME )
            FIRSTIME = .FALSE.
        END IF
        TKPRO = 0.0
        TMPRO = 0.0
        TKREFHR = 0.0

C.........  Replace ref county min/max and RH with fuelmonth min/max and RH
C.........  Choose month-specific fulemonth county
        L = FIND1FIRST( COUNTY, NREFF, FMREFSORT( :,1 ) )
        K = FIND1FIRST( COUNTY, NFUELC,FMREFLIST( :,1 ) )

        NMON = FMREFLIST( K, 2 )   ! no month of ref county
        FUELMONTH = 0
C.........  Loop over months per ref. county
        DO J = L, L + NMON - 1

            CURMONTH  = FMREFSORT( J,3 )    ! processing  current month per ref. county

C.............  Skip other months
            IF( CURMONTH  == PMONTH ) THEN
                FUELMONTH = FMREFSORT( J,2 )    ! processing fuelmonth/county

C.................  Define fuelmonth min/max and RH per ref. county
                NR = FIND1( COUNTY,  NREFC,MCREFIDX (:,1) )
                NF = FIND1( CURMONTH,NFUEL,FUELIDX (NR,:) )

                RHREFAVG = RHFUEL  ( NR,NF )
                MAXTREF  = MAXTFUEL( NR,NF )
                MINTREF  = MINTFUEL( NR,NF )
c          print*,COUNTY,NR,NF,RHREFAVG,MAXTREF,MINTREF,FUELMONTH,'REF,,,'
            END IF

        END DO

C.........  Write out last ref. county min/max temp and avg RH
        WRITE( ODEV,94040 ) COUNTY, FUELMONTH, 'min_max',
     &                      RHREFAVG, MINTREF, MAXTREF

        IMAXT = INT( MAXTREF )
        IMINT = INT( MINTREF )

        TKMAX  = MAXVAL( THOUR )
        TKMIN  = MINVAL( THOUR )
        TKDIF  = ABS( TKMAX - TKMIN )
        TKMED  = ( TKMAX + TKMIN ) / 2
        
        TMPRO = 0.0
        DO T = 1, 24
            TMPRO( T ) = ( THOUR( T ) - TKMED ) / TKDIF
        END DO

C.........  Calculate temp profiles based on a combination of temp bins
        MAXT = IMAXT + ( PPTEMP - MOD( IMAXT,PPTEMP ) )
        MINT = IMINT - MOD( IMINT,PPTEMP )

        IT = 0
C.........  Determine temperature bins based on PPTEMP
        DO MX = MAXT,MINT,-PPTEMP

            DO MN = MINT,MX,PPTEMP

                IT = IT + 1
                DT = MX - MN
                WRITE( TPROID,94030 ) HDR, MDATE, IT

                TKPRO = 0.0
C.................  Create 24 hr temp profile per temp bin
                DO T = 1,24
                    TKPRO( T ) = TMPRO( T ) * DT + (MX+MN)/2
                END DO
                    
C.................  Output for SMOKE ready input file
                WRITE( ODEV,94050 ) COUNTY, FUELMONTH, 
     &                 TRIM(TPROID), RHREFAVG, (TKPRO( T ), T=1,24)

            END DO

        END DO
        
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I6, :, 1X ) )

94020   FORMAT( A, 4( 1X, F8.2, 1X, A ) )

94030   FORMAT( A, I7.7, I3.3 )   

94040   FORMAT( I6.6, I5, 3X, A, 3F10.2 )   

94050   FORMAT( I6.6, I5, 3X, A, 25F10.2 )   

94060   FORMAT( I6.6, I5, 3X, A, F10.2, 2I5 )   
 
        END SUBROUTINE WRTEMPROF
