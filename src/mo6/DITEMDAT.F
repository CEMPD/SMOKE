      MODULE DITEMDAT 
C    This block has the data for temperature profile, vapor pressure 
C    parameters, and temperature related parameters and results used in 
C    the Hourly Diurnal Fraction subroutine, DIRNLFRC.FOR and the 
C    diurnal temperature subroutine, DIRNLTEM.FOR 
C 
C 
C  Changes:  (Last change first) 
C 
C 23 Jun 00 @EPA-bag bug67  Renamed NEW_TEMP to DIRNL_TEMP 
C  07 Jan 00 @Dyntel-bg 1-005 updated diurnal code(new) 
C 
C 
C  MODULE dictionary: 
C 
C     Name           Type              Description 
C    ------          ----  --------------------------------------------- 
C  LESS_THN10_FCTR    R    ratio for diurnal emissions which have a daily temperature 
C                          range of less than 10 degrees fahrenheit 
C  TEMP               R    standard temperature distribution 
C  MIN_TEMPST         R    minimum temperature from standard temperature distribution 
C  MAX_TEMPST         R    maximum temperature from standard temperature distribution 
C  LOWTEMP            R    diurnal emission low temperature limit, 40 degrees fahrenheit 
C  MIN_TEMP_TIME      I    time when minimum temperature occurs 
C  DIRNL_TEMP           R    diurnal temperature distribution derived from the user supplied temperatures 
C  CTEMPCHANGE        R    temperature change, current hour 
C  TTEMPCHANGE        R    temperature change, from beginning of diurnal emissions 
C  PTEMPCHANGE        R    temperature change, previous hour 
C  VP                 R    vapor pressure coefficients 
C  VPLO               R    vapor pressure at minimum temperature 
C  VPHI               R    vapor pressure at maximum temperature 
C  VPMI               R    average of VPHI and VPLO 
C  VPPR               R    product of VPMI and (VPHI-VPLO) 
C 
C  Local array subscripts: 
C 
C TEMP(24)            -      TEMP(IH) 
C CTEMPCHANGE(24)     -      CTEMPCHANGE(IH) 
C PTEMPCHANGE(24)     -      PTEMPCHANGE(IH) 
C TTEMPCHANGE(24,20)  -      TTEMPCHANGE(IH,ISOAK) 
C VP(5)               -      5 parameters for vapor pressure equation 
C 
C 
      IMPLICIT NONE 
      INTEGER I,II 
C 
C For standard temperature profile 
C 
C IH         :     1            2            3           4          ...      24 
C clock time : 12a.m.-1a.m.  1a.m.-2a.m.  2a.m.-3a.m.  3a.m.-4a.m.  ...   11p.m.-12a.m. 
C 
      REAL, DIMENSION(24), PARAMETER :: TEMP = (/ 
C IH =   1     2     3     4     5     6     7     8     9    10 
     *  75.8, 74.7, 73.9, 73.3, 72.6, 72.0, 72.5, 75.5, 80.3, 85.2, 
C IH =  11    12    13     14   15    16    17    18    19    20 
     *  89.4, 93.1, 95.1, 95.8, 96.0, 95.5, 94.1, 91.7, 88.6, 85.5, 
C IH =  21    22    23    24 
     *  82.8, 80.9, 79.0, 77.2/) 
C 
      REAL, PARAMETER :: MIN_TEMPST       = 72.0 !standard minimum temperature from 
                                                     !standard temperature distribution 
C 
      REAL, PARAMETER :: MAX_TEMPST       = 96.0 !standard maximum temperature from 
                                                     !standard temperature distribution 
C 
      REAL, PARAMETER :: LOWTEMP          = 40.0 !diurnal emission low temperature limit 
C 
      INTEGER, PARAMETER :: MIN_TEMP_TIME = 6    !the standard minimum temperature occurs between 
                                                 !the hours of 5a.m. and 6a.m. 
C 
      REAL, DIMENSION(24) :: DIRNL_TEMP  = 0.0 
C 
      REAL, DIMENSION(24) :: CTEMPCHANGE = -1.0 
      REAL, DIMENSION(24) :: PTEMPCHANGE = -1.0 
      REAL, DIMENSION(24,20) :: TTEMPCHANGE = -1.0 
C 
C The coefficient for VP_midpoint 
C 
      REAL, DIMENSION(5), PARAMETER :: VP = 
     *      (/658728.1,-0.2267976,126.92048,6417.4874,459.67/) 
C 
      REAL :: VPLO 
      REAL :: VPHI 
      REAL :: VPMI 
      REAL :: VPPR 
      REAL :: LESS_THN10_FCTR 
C 
      END MODULE DITEMDAT 
 
