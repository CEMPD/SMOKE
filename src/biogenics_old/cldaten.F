
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE  CLDATEN( JDATE, JTIME, LAT, LON, CLDATN )

C***********************************************************************
C  subroutine body starts at line
C
C  FUNCTION:
C       Computes gridded cloud attenuation factor CLDATN for SMOKE BEIS2
C
C       Reads records from the SURMET-formatted data file
C       then call a BARNES analysis scheme to compute the
C       gridded cloud attenuation factors.
C
C       in-lines computation of cos( zenith angle ) originally found
C       in ROM and UAM routine SOLELV()
C
C  PRECONDITIONS REQUIRED:
C	SURMET-formatted surface meteorology file.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       m3io
C       CBARNES1, READSMET
C
C  REVISION  HISTORY:
C       Prototype 12/95 by CJC:  dumb implementation using EDSS clouds
C
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

      INCLUDE 'CONST3.EXT'      ! physical and mathematical constants
      INCLUDE 'PARMS3.EXT'      ! I/O API constants
      INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
      INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
      INCLUDE 'GRDIMS3.EXT'     ! grid-related constants
      INCLUDE 'BIODIMS3.EXT'    ! biogenic-related constants


C...........   ARGUMENTS and their descriptions:

        INTEGER         JDATE   	!  current simulation date (YYYYDDD)
        INTEGER         JTIME   	!  current simulation time (HHMMSS)
        REAL            LAT   ( NCOLS, NROWS )	!  latitude,  deg
        REAL            LON   ( NCOLS, NROWS )	!  longitude, deg
        REAL            CLDATN( NCOLS, NROWS )	!  cloud attenuation factor

C...........   PARAMETERS and their descriptions:

        INTEGER		MAXMET
        REAL            SIGA, D60, D15, D24, D100
        REAL            ROTDAY , ROTHR, ROTLON, SDEC
        REAL            WL50

        PARAMETER ( 
     &          MAXMET = 1000 ,
     &            SIGA = 279.9348  ,
     &            SDEC =   0.39784984 , ! SIN (23^26'37.8") the declination angle
     &          WL50   = 273.83678 ,
     &            D100 = 1.0 /100.0,
     &             D60 = 1.0 / 60.0,
     &             D15 = 1.0 / 15.0 ,
     &             D24 = 1.0 / 24.0 ,
     &          ROTDAY = 360.0 / 365.242,   ! rotation per day
     &          ROTHR  = D24 * ROTDAY,      ! rotation per hour
     &          ROTLON = PI180 * ROTHR )    ! rot (radians) per degree lon


C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         PROMPTFFILE
        LOGICAL         READSMET

        EXTERNAL        PROMPTFFILE,  READSMET


C...........   SCRATCH LOCAL VARIABLES and their descriptions:
C.......   Surface meteorology data from READSMET:
C.......   has the following structure
C              RMET (  1, * )  --  latitude
C              RMET (  2, * )  --  longitude
C              RMET (  3, * )  --  opaque sky cover
C              RMET (  4, * )  --  total sky cover descriptor
C              RMET (  5, * )  --  lowest cloud cover
C              RMET (  6, * )  --  lowest cloud height
C              RMET (  7, * )  --  2nd lowest cloud cover
C              RMET (  8, * )  --  2nd lowest cloud height
C              RMET (  9, * )  --  3rd lowest cloud cover
C              RMET ( 10, * )  --  3rd lowest cloud height
C              RMET ( 11, * )  --  sea level pressure
C              RMET ( 12, * )  --  wind direction
C              RMET ( 13, * )  --  wind speed
C              RMET ( 14, * )  --  temperature
C              RMET ( 15, * )  --  dew point
C              RMET ( 16, * )  --  station pressure
C.......   Values for the total sky cover descriptor are as follows:
C              0 = clear
C              1 = partial obscuration
C              2 = thin scattered
C              3 = thin broken
C              4 = thin overcast
C              5 = scattered
C              6 = broken
C              7 = overcast
C              8 = obscured
C           -999 = indeterminate
C.........................................................................


        INTEGER		NMET		     !  actual number of stations found
        INTEGER		IMET( MAXMET )	     !  array of station ID's
        REAL		RMET( 16 , MAXMET )  !  station data
        REAL            SLON( MAXMET )       !  longitudes for CGRID1
        REAL            SLAT( MAXMET )       !  latitudes  for CGRID1
        REAL            ATEN( MAXMET )       !  cloud attenuation

C.......   loop counters, etc:

        REAL            ALAT, ALON
        INTEGER         I, M, N, NT, R, C

C.......   Variables for zenith angle calculation:

        REAL            TIMLOC, T, DADM, DF
        REAL            DESIN, DECOS, DESIN2, DECOS2 
        REAL            SIG, DECSIN, DECCOS, EQT, TST, HRANGL, CZ

C.......   Variables for cloud attenuation calculation:

        INTEGER		DESCRI
        REAL            COVER, TTYPE
        REAL		PCT( 3 )
        REAL            HGT( 3 )
        
C...........   SAVED LOCAL VARIABLES and their descriptions:
C...........   NOTE:  the ANSI standard requires the use of SAVE statements
C...........   for variables which must retain their values from call to call.

        INTEGER         SDEV		!  unit number for surface met file
        LOGICAL         FIRSTIME
        REAL            NBORD, SBORD, EBORD, WBORD
        REAL            AA( 5 ), BB( 5 )

        DATA            FIRSTIME / .TRUE. /
        DATA            AA / 0.165, 0.18,  0.365,  0.47,  0.89 /
        DATA            BB / 0.005, 0.02, -0.015, -0.01, -0.04 /

        SAVE   SDEV, FIRSTIME, NBORD, SBORD, EBORD, WBORD, AA, BB


C***********************************************************************
C   begin body of subroutine  CLDATEN

        IF ( FIRSTIME ) THEN

            FIRSTIME = .FALSE.

            CALL M3MESG( 'Using OBSERVATIONAL MET version of CLDATEN' )
            SDEV = PROMPTFFILE( 
     &              'Enter name for SURFACE MET input file',
     &              .TRUE., .TRUE., 'SURMET', 'TMPBIO:CLDATEN' )

            NBORD = MAX( LAT( 1, NROWS ), LAT( NCOLS, NROWS ) ) + 10.0
            SBORD = MIN( LAT( 1,     1 ), LAT( NCOLS,     1 ) ) - 10.0
            EBORD = MAX( LON( NCOLS, 1 ), LON( NCOLS, NROWS ) ) + 10.0
            WBORD = MIN( LON(     1, 1 ), LON(     1, NROWS ) ) - 10.0

        END IF          	!  if firstime

        IF ( .NOT. READSMET( SDEV, JDATE, JTIME, 
     &                       NBORD, SBORD, EBORD, WBORD,
     &                       MAXMET, NMET  , IMET  , RMET ) ) THEN

            CALL M3EXIT( 'CLDATEN', JDATE, JTIME,
     &                   'Error reading SURFACE MET file', 2 )

        END IF			!  if readsmet() failed

C.......   Compute each station's cloud attenuation and locations,
C.......   ordered for CBARNES1:

        TIMLOC = ROTHR * ( FLOAT( JTIME / 10000 ) +
     &                   D60 * ( FLOAT( MOD( JTIME / 100 , 100 ) ) +
     &                           D60 * FLOAT( MOD( JTIME , 100 ) ) ) )
        DADM   = ROTDAY * FLOAT( MOD( JDATE , 1000 ) ) + TIMLOC
        NT     = 0

        DO  33  N = 1, NMET

            ALAT = RMET(  1,N )
            ALON = RMET(  2,N )

            SLAT( N ) = ALAT
            SLON( N ) = ALON

C...........   Compute CZ = cos( zenith angle ) (from SOLELV() )

            T      = D15 * ALON
            DF     = ROTLON * T +  DADM
            DESIN  =  SIN( DF )             !  SINE   of this angle
            DECOS  =  COS( DF )             !  COSINE of this angle
            DESIN2 =  2.0 * DESIN * DECOS            !  SINE   of double-angle
            DECOS2 =  DECOS * DECOS - DESIN * DESIN  !  COSINE of double-angle
            SIG    =  SIGA + DF
     &                     + 1.914827 * DESIN  - 0.079525 * DECOS
     &                     + 0.019938 * DESIN2 - 0.00162  * DECOS2

C...........   The equation of time adjustment, uses sine of the declination

            EQT  =  0.123470 * DESIN  - 0.004289 * DECOS
     &            + 0.153809 * DESIN2 + 0.060783 * DECOS2

C...........   The true solar time

            TST    = T + TIMLOC - EQT
            HRANGL = 15.0 * PI180 * ( TST - 12.0 )

C...........   Compute the sine of the solar elevation (= cos of zenith angle)

            DECSIN = SDEC * SIN( PI180 * SIG )
            DECCOS = SQRT( 1.0 - DECSIN )

            CZ = DECSIN * SIN( ALAT ) + 
     &           DECCOS * COS( ALAT ) * COS( HRANGL )


C...........   Now compute this stations cloud attenuation, in terms of
C...........   its sky descriptor:

            DESCRI = NINT( RMET( 4,N ) )

            IF ( DESCRI .LT. 0 )  THEN

                GO TO  33

            ELSE IF ( DESCRI .EQ. 0 ) THEN

                NT = NT + 1
                SLAT( NT ) = ALAT
                SLON( NT ) = ALON
                ATEN( NT ) = 1.0

            ELSE IF ( DESCRI .EQ. 1 ) THEN

C...........   if the sky descriptor is 1 or 8, then a ground level obscuration
C...........   is indicated.   There is not a differentiation made between 
C...........   "thin" and otherwise, but if the sky is not totally (100%)
C...........   obscured, the partial obscuration (=1) is given.  In this
C...........   case there will be a value in the opaque sky cover column
C...........   between 10 and 90.

                NT = NT + 1
                TTYPE = AA(1) + BB(1) / MAX( CZ, 0.2 )
                COVER = RMET( 3,N )
                ATEN( NT ) = 1.0 - D100 * COVER * ( 1.0 - TTYPE )
                SLAT( NT ) = ALAT
                SLON( NT ) = ALON

            ELSE IF ( DESCRI .LE. 7 ) THEN
                
C...............   if clouds are reported and the sky is not obscured, then
C...............   compute the transmittance for all layers.
C...............   Take cloud layer information one layer at a time.
C...............   Ccheck the height of the cloud layer to decide which
C...............   equation should be used.  compute the transmittanc
C...............   for the type of cloud cover given.  Then compute the
C...............   layer transmittance including the fractional coverage
C...............   of the layer.  
C...............   If the sky descriptor indicated 'thin' sky cover (=2,3,4)
C...............   adjust the total transmittance half way back to 1.0

                NT = NT + 1
                HGT( 1 ) = RMET( 6,N )
                IF ( HGT( 1 ) .LT. 0.0 ) GO TO  33
                
                CZ    = 1.0 / MAX( CZ, 0.2 )
                COVER = RMET( 3,N )
                PCT( 1 ) = RMET(  5,N )
                HGT( 2 ) = RMET(  7,N )
                PCT( 2 ) = RMET(  8,N )
                HGT( 3 ) = RMET(  9,N )
                PCT( 3 ) = RMET( 10,N )

                IF ( HGT( 2 ) .LT. 0.0 ) THEN
                    M = 1
                ELSE IF ( HGT( 3 ) .LT. 0.0 ) THEN
                    M = 2
                ELSE
                    M = 3
                END IF

                T = 1.0
                DO  11  I = 1, M
                    IF      ( HGT(I) .LT. 1000.0 ) THEN
                        TTYPE = AA(2) + BB(2) * CZ
                    ELSE IF ( HGT(I) .LT. 1000.0 ) THEN
                        TTYPE = AA(3) + BB(3) * CZ
                    ELSE IF ( HGT(I) .LT. 1000.0 ) THEN
                        TTYPE = AA(4) + BB(4) * CZ
                    ELSE
                        TTYPE = AA(5) + BB(5) * CZ
                    END IF
                    T = T * ( 1.0 - D100 * ( 1.0 - TTYPE ) *  PCT( I ) )
11              CONTINUE

                IF ( DESCRI .LE. 4 ) THEN
                    ATEN( NT ) = 0.5 + 0.5 * T
                    SLAT( NT ) = ALAT
                    SLON( NT ) = ALON
                ELSE
                    ATEN( NT ) = T
                    SLAT( NT ) = ALAT
                    SLON( NT ) = ALON
                END IF

            ELSE		!  else sky descriptor is 8

C...............   If the sky descriptor is 8, then a ground level obscuration
C...............   is indicated.   Descriptors of 8 are always accompanied by
C...............   an opaque sky cover value of 100.

                NT = NT + 1
                ATEN( NT ) = AA(1) + BB(1) / MAX( CZ, 0.2 )
                SLAT( NT ) = ALAT
                SLON( NT ) = ALON

            END IF	!  if descri zero, <= 7, or is 8

33      CONTINUE		!  end loop on stations N


C...........   Grid the sky coverage fraction:

        CALL CBARNES1( NCOLS * NROWS, LAT, LON, 
     &                 NT, SLAT, SLON, ATEN, WL50, CLDATN )


C...........   Trap the cloud attenuation factor between 0.05 and 1.0:
        
        DO  55  R = 1, NROWS
        DO  44  C = 1, NCOLS
           T = CLDATN( C,R )
           IF ( T .LT. 0.05 ) THEN
               CLDATN( C,R ) = 0.05
           ELSE IF ( T .GT.  1.0 ) THEN
               CLDATN( C,R ) = 1.0
           END IF
44      CONTINUE
55      CONTINUE

        RETURN

        END

