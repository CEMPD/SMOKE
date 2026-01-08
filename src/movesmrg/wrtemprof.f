
        SUBROUTINE WRTEMPROF( ODEV1, ODEV2, MYEAR, COUNTY, PMONTH,
     &                        PDTEMP, PPTEMP, THOUR, MAXTEMP,
     &                        MINTEMP, TEMPBIN, MINNORH )

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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***********************************************************************

C...........   MODULES for public variables
C...........   This module is the derived meteorology data for emission factors
        USE M3UTILIO

        USE MODMET, ONLY:   RHTBIN, NRHTBIN
        
        USE MODMBSET, ONLY: NREFC, MCREFIDX

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS 
C       CHARACTER(2)    CRLF
C       CHARACTER(14)   MMDDYY
C       INTEGER         FINDC

C        EXTERNAL     CRLF, MMDDYY, FINDC

C...........   SUBROUTINE ARGUMENTS
C.............  Subroutine arguments
        INTEGER     , INTENT(IN)  ::  ODEV1           ! MOVES RPP output file
        INTEGER     , INTENT(IN)  ::  ODEV2           ! MOVES RPD/RPV output file
        INTEGER     , INTENT(IN)  ::  MYEAR           ! modeling year
        CHARACTER(*), INTENT(IN)  ::  COUNTY          ! refCounty
        INTEGER     , INTENT(IN)  ::  PMONTH          ! fuel month
        INTEGER     , INTENT(IN)  ::  PDTEMP          ! RPP/RPV temperature increment
        INTEGER     , INTENT(IN)  ::  PPTEMP          ! RPP temperature increment
        REAL        , INTENT(IN)  ::  THOUR( 24 )     ! avg temp 24-hr profiles for refcounty/fuelmonth
        REAL        , INTENT(IN)  ::  MAXTEMP         ! max temp for refcounty/fuelmonth
        REAL        , INTENT(IN)  ::  MINTEMP         ! min temp for refcounty/fuelmonth
        REAL        , INTENT(IN)  ::  TEMPBIN         ! temp buffer
        INTEGER     , INTENT(IN)  ::  MINNORH         ! min no of data points for avg RH by tempbin

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
        REAL    RHAVG                     ! tmp relative humidity

        LOGICAL FIRSTIME 

        CHARACTER(32)    TPROID             ! temporal resolution header
        CHARACTER(300)   MESG               ! message buffer
        CHARACTER(16) :: PROGNAME = 'WRTEMPROF' ! program name

C***********************************************************************
C   begin body of subroutine WRTEMPROF

C.........  Write out last ref. county min/max temp and avg RH
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
                WRITE( TPROID,94040 ) 'M', MYEAR, PMONTH, IT

                TKPRO = 0.0
C.................  Create 24 hr temp profile per temp bin
                DO T = 1,24
                    TKPRO( T ) = TMPRO( T ) * DT + ST/2.0
                END DO

C.................  Output for SMOKE ready input file
                WRITE( ODEV1,94050 ) COUNTY, PMONTH, 
     &              TRIM(TPROID), MINTEMP, MAXTEMP, (TKPRO( T ), T=1,24)

            END DO

        END DO

C.........  Output avg RH by temperature bins for RPD/RPV modes
C.........  Fill empty temperature bins for RPD/RPV modes
        NR = FINDC( COUNTY,NREFC, MCREFIDX( :,1 ) )
        NF = PMONTH
        NT  = 0
        NTB = 0
        DO T = -150, 200-PDTEMP, PDTEMP
            NT = NT + 1
            IF( NRHTBIN( NR,NF,NT ) < MINNORH ) THEN
                RHTBIN ( NR,NF,NT ) = 0.0
                NRHTBIN( NR,NF,NT ) = 0
            END IF
        END DO
        
        NT = 0
        FIRSTIME = .TRUE.
        DO T = -150, 200-PDTEMP, PDTEMP
            NT = NT + 1
            IF( FIRSTIME .AND. NRHTBIN( NR,NF,NT ) == 0 ) THEN
                NTB = NT
                FIRSTIME = .FALSE.
                CYCLE
            ELSE IF( RHTBIN( NR,NF,NT ) > 0.0 ) THEN
                RHTBIN ( NR,NF,NTB:NT ) = RHTBIN ( NR,NF,NT )
                NRHTBIN( NR,NF,NTB:NT ) = NRHTBIN( NR,NF,NT )
                NTB = NT + 1
                FIRSTIME = .TRUE.
                CYCLE
            ELSE IF( T == 200-PDTEMP ) THEN
                RHTBIN ( NR,NF,NTB:NT ) = RHTBIN ( NR,NF,NTB-1 )
                NRHTBIN( NR,NF,NTB:NT ) = NRHTBIN( NR,NF,NTB-1 )
            END IF
        END DO

C.........  Write out refcounty min/max temp and avg RH by temperature bin for RPD/RPV
C.........  Calculate max/min temp bins based on RPD_TEMP_INCREMENT 
        MAXT = IMAXT + ( PDTEMP - ABS( MOD( IMAXT,PDTEMP ) ) )
        MINT = IMINT - ABS( MOD( IMINT,PDTEMP ) )

        IF( MINT <= 0 .AND. MINTEMP < 0  ) THEN
            MINT = IMINT - ( PDTEMP - ABS( MOD( IMINT,PDTEMP ) ) )
        END IF

        NT = 0
        DO T = -150, 200-PDTEMP, PDTEMP
            NT = NT + 1
            IF( MINT <= T .AND. T <= MAXT ) THEN
            
                RHAVG = RHTBIN(NR,NF,NT) / NRHTBIN(NR,NF,NT)

                WRITE( ODEV2,94060 ) COUNTY, PMONTH, RHAVG, MINTEMP,
     &                 MAXTEMP, T

           END IF
        END DO
        
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94040   FORMAT( A, I4.4, I2.2, I4.4 )   

94050   FORMAT( A,',', I5,',', 3X, A, 26(',',F10.2) )   

94060   FORMAT( A,',',I5,',',F12.6,',',F10.2,',',F10.2,',',I5 )   
 
        END SUBROUTINE WRTEMPROF
