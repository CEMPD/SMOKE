
        SUBROUTINE WRTEMPROF( ODEV, MYEAR, HDR, COUNTY, PMONTH,
     &                        PPTEMP, RHAVG, THOUR, MAXTEMP, MINTEMP ) 

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

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS 
        CHARACTER(2)    CRLF
        CHARACTER(14)   MMDDYY

        EXTERNAL     CRLF, MMDDYY

C...........   SUBROUTINE ARGUMENTS
C.............  Subroutine arguments
        INTEGER  , INTENT(IN)  ::  ODEV
        INTEGER  , INTENT(IN)  ::  MYEAR
        CHARACTER, INTENT(IN)  ::  HDR
        INTEGER  , INTENT(IN)  ::  COUNTY
        INTEGER  , INTENT(IN)  ::  PMONTH
        INTEGER  , INTENT(IN)  ::  PPTEMP
        REAL     , INTENT(IN)  ::  RHAVG
        REAL     , INTENT(IN)  ::  THOUR( 24 )
        REAL     , INTENT(IN)  ::  MAXTEMP 
        REAL     , INTENT(IN)  ::  MINTEMP

C.........  Local array
        REAL   :: TKPRO( 24 ) = 0.0 !  24hr temp. profile
        REAL   :: TMPRO( 24 ) = 0.0 !  tmp 24hr temp. profile

C...........   Other local variables
        INTEGER I, IT, J, K, L, LT, N, NR, MX, MN, R, S, T           ! counters and indices                     
        
        INTEGER IOS                       ! I/O status
        INTEGER IMAXT, IMINT              ! current integer min/max temp.
        INTEGER MAXT , MINT               ! tmp remainder of max/min

        REAL    TKMAX                     ! tmp max temperatures
        REAL    TKMIN                     ! tmp min temperatures
        REAL    TKDIF, DT                 ! tmp DIFF of min/max temperatures
        REAL    TKMED                     ! tmp median temperatures

        CHARACTER(32)    TPROID             ! temporal resolution header
        CHARACTER(300)   MESG               ! message buffer
        CHARACTER(16) :: PROGNAME = 'WRTEMPROF' ! program name

C***********************************************************************
C   begin body of subroutine WRTEMPROF

C.........  Write out last ref. county min/max temp and avg RH
        WRITE( ODEV,94040 ) COUNTY, PMONTH, 'min_max',
     &                      RHAVG, MINTEMP, MAXTEMP

        IMAXT = INT( MAXTEMP )
        IMINT = INT( MINTEMP )

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
        IF( MINT <= 0 .AND. MINTEMP < 0  ) MINT = MINT - PPTEMP
        
        IT = 0
C.........  Determine temperature bins based on PPTEMP
        DO MX = MAXT,MINT,-PPTEMP

            DO MN = MINT,MX,PPTEMP

                IT = IT + 1
                DT = MX - MN
                WRITE( TPROID,94030 ) HDR, MYEAR, PMONTH, IT

                TKPRO = 0.0
C.................  Create 24 hr temp profile per temp bin
                DO T = 1,24
                    TKPRO( T ) = TMPRO( T ) * DT + (MX+MN)/2
                END DO
                    
C.................  Output for SMOKE ready input file
                WRITE( ODEV,94050 ) COUNTY, PMONTH, 
     &                 TRIM(TPROID), RHAVG, (TKPRO( T ), T=1,24)

            END DO

        END DO
        
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I6, :, 1X ) )

94020   FORMAT( A, 4( 1X, F8.2, 1X, A ) )

94030   FORMAT( A, I4.4, I2.2, I4.4 )   

94040   FORMAT( I6.6, I5, 3X, A, 3F10.2 )   

94050   FORMAT( I6.6, I5, 3X, A, 25F10.2 )   

94060   FORMAT( I6.6, I5, 3X, A, F10.2, 2I5 )   
 
        END SUBROUTINE WRTEMPROF
