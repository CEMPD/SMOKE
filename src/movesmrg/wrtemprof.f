
        SUBROUTINE WRTEMPROF( ODEV1, ODEV2, MYEAR, HDR, COUNTY, PMONTH,
     &                        PDTEMP, PPTEMP, RHAVG, THOUR, MAXTEMP,
     &                        MINTEMP, TEMPBIN )

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
C...........   This module is the derived meteorology data for emission factors
        USE MODMET, ONLY:   RHTBIN, NRHTBIN
        
        USE MODMBSET, ONLY: NREFC, MCREFIDX

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS 
        CHARACTER(2)    CRLF
        CHARACTER(14)   MMDDYY
        INTEGER         FIND1

        EXTERNAL     CRLF, MMDDYY, FIND1

C...........   SUBROUTINE ARGUMENTS
C.............  Subroutine arguments
        INTEGER  , INTENT(IN)  ::  ODEV1           ! MOVES RPP output file
        INTEGER  , INTENT(IN)  ::  ODEV2           ! MOVES RPD/RPV output file
        INTEGER  , INTENT(IN)  ::  MYEAR           ! modeling year
        CHARACTER, INTENT(IN)  ::  HDR             ! Averaging type for RPP output file
        INTEGER  , INTENT(IN)  ::  COUNTY          ! refCounty
        INTEGER  , INTENT(IN)  ::  PMONTH          ! fuel month
        INTEGER  , INTENT(IN)  ::  PDTEMP          ! RPP/RPV temperature increment
        INTEGER  , INTENT(IN)  ::  PPTEMP          ! RPP temperature increment
        REAL     , INTENT(IN)  ::  RHAVG           ! avg RH for refcounty/fuelmonth
        REAL     , INTENT(IN)  ::  THOUR( 24 )     ! avg temp 24-hr profiles for refcounty/fuelmonth
        REAL     , INTENT(IN)  ::  MAXTEMP         ! max temp for refcounty/fuelmonth
        REAL     , INTENT(IN)  ::  MINTEMP         ! min temp for refcounty/fuelmonth
        REAL     , INTENT(IN)  ::  TEMPBIN         ! temp buffer

C.........  Local array
        REAL   :: TKPRO( 24 ) = 0.0 !  24hr temp. profile
        REAL   :: TMPRO( 24 ) = 0.0 !  tmp 24hr temp. profile

C...........   Other local variables
        INTEGER I, IT, J, K, L, LT, N, NF, NR, NT, NTB, MX, MN, R, S, T          ! counters and indices                     
        
        INTEGER IOS                       ! I/O status
        INTEGER IMAXT, IMINT              ! current integer min/max temp.
        INTEGER MAXT , MINT               ! tmp remainder of max/min

        REAL    TKMAX                     ! tmp max temperatures
        REAL    TKMIN                     ! tmp min temperatures
        REAL    TKDIF, DT, ST             ! tmp DIFF of min/max temperatures
        REAL    TKMED                     ! tmp median temperatures
        REAL    RHVAL                     ! tmp relative humidity

        CHARACTER(32)    TPROID             ! temporal resolution header
        CHARACTER(300)   MESG               ! message buffer
        CHARACTER(16) :: PROGNAME = 'WRTEMPROF' ! program name

C***********************************************************************
C   begin body of subroutine WRTEMPROF

C.........  Write out last ref. county min/max temp and avg RH
        WRITE( ODEV1,94040 ) COUNTY, PMONTH, 'min_max',
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
        MAXT = IMAXT + ( PPTEMP - ABS( MOD( IMAXT,PPTEMP ) ) )
        MINT = IMINT - ABS( MOD( IMINT,PPTEMP ) )

        IF( MINT <= 0 .AND. MINTEMP < 0  ) THEN
            MINT = IMINT - ( PPTEMP - ABS( MOD( IMINT,PPTEMP ) ) )
        END IF

        IT = 0
C.........  Determine temperature bins based on PPTEMP
        DO MX = MAXT,MINT,-PPTEMP

            DO MN = MINT,MX,PPTEMP

                IT = IT + 1
                DT = MX - MN
                ST = MX + MN
                WRITE( TPROID,94030 ) HDR, MYEAR, PMONTH, IT

                TKPRO = 0.0
C.................  Create 24 hr temp profile per temp bin
                DO T = 1,24
                    TKPRO( T ) = TMPRO( T ) * DT + ST/2.0
                END DO

C.................  Output for SMOKE ready input file
                WRITE( ODEV1,94050 ) COUNTY, PMONTH, 
     &                 TRIM(TPROID), RHAVG, (TKPRO( T ), T=1,24)

            END DO

        END DO

C.........  Output avg RH by temperature bins for RPD/RPV modes

C.........  Calculate max/min temp bins based on RPD_TEMP_INCREMENT 
        MAXT = IMAXT + ( PDTEMP - ABS( MOD( IMAXT,PDTEMP ) ) )
        MINT = IMINT - ABS( MOD( IMINT,PDTEMP ) )

        IF( MINT <= 0 .AND. MINTEMP < 0  ) THEN
            MINT = IMINT - ( PDTEMP - ABS( MOD( IMINT,PDTEMP ) ) )
        END IF

        NR = FIND1( COUNTY,NREFC, MCREFIDX( :,1 ) )
        NF = PMONTH
        NT = 0
        DO T = -150, 150, PDTEMP
             NT = NT + 1
             IF( MINT <= T .AND. T <= MAXT ) THEN

                  IF( NRHTBIN( NR,NF,NT ) < 1 ) THEN 
                      NTB = NT + INT( TEMPBIN/PDTEMP )
                      IMINT = MINT + INT( TEMPBIN )
                      IF( T <= IMINT ) THEN
                          IF( NRHTBIN( NR,NF,NTB ) < 1 ) NTB = NTB + 1 
                          RHTBIN( NR,NF,1:NTB-1 ) = RHTBIN( NR,NF,NTB ) 
                          NRHTBIN( NR,NF,1:NTB-1 ) = NRHTBIN( NR,NF,NTB ) 
                      ELSE
                          RHTBIN( NR,NF,NT: ) = RHTBIN( NR,NF,NT-1 ) 
                          NRHTBIN( NR,NF,NT: ) = NRHTBIN( NR,NF,NT-1 ) 
                      END IF
                  END IF

                  RHVAL = RHTBIN(NR,NF,NT) / NRHTBIN(NR,NF,NT)

C...................  Write out refcounty min/max temp and avg RH by temperature bin for RPD/RPV
                  WRITE( ODEV2,94060 ) COUNTY, PMONTH, RHVAL, MINTEMP,
     &                                 MAXTEMP, T

             END IF
        END DO
        
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I6, :, 1X ) )

94020   FORMAT( A, 4( 1X, F8.2, 1X, A ) )

94030   FORMAT( A, I4.4, I2.2, I4.4 )   

94040   FORMAT( I6.6, I5, 3X, A, 3F10.2 )   

94050   FORMAT( I6.6, I5, 3X, A, 25F10.2 )   

94060   FORMAT( I6.6,',',I5,',',F12.6,',',F10.2,',',F10.2,',',I5 )   
 
        END SUBROUTINE WRTEMPROF
