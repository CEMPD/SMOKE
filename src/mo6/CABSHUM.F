      REAL FUNCTION CABSHUM(IH)
C
C  This function computes the absolute humidity (in grains/pound)
C  given the temperature (F), relative humidity (%), and barometric
C  pressure (inches Hg).
C
C  Changes: (Last change first)
C
C  31 Jul 02 BUG394 EPA -bag Equation to convert from relative humidity
C            to absolute humidity was coded incorrectly
C  02 Jul 02 AIR Task 36: New function.
C
C Input on call:
C
C   argument list : 
C
      IMPLICIT NONE
C
      INCLUDE 'ALUHIN.I'   !RH_HRLY,GBARO
      INCLUDE 'TEMPS.I'    !TEMHRLY

      REAL, DIMENSION(6), SAVE :: C=(/  ! constant terms, see Task 36 write-up/spreadsheet
     * 4.34780E+03,  ! C(1)
     * 1.53195E-04,  ! C(2)
     * 3.24370E+00,  ! C(3)
     * 5.87826E-03,  ! C(4)
     * 1.17020E-08,  ! C(5)
     * 2.18785E-03/) ! C(6)
C
      INTEGER IH      ! Hour index
      REAL BETA       ! Beta term of equation, see Task 36 write-up/spreadsheet
      REAL BP         ! Barometric pressure (default or user-supplied)
      REAL PV         ! PV term of equation, see Task 36 write-up/spreadsheet
      REAL T0         ! T0 term of equation, see Task 36 write-up/spreadsheet
      REAL TK         ! TK term of equation, temp in degrees Kelvin, see Task 36 write-up/spreadsheet
C
C  Set default or user supplied barometric pressure
C
      IF(BARO_FLAG) THEN
        BP=BARO  ! user-supplied
      ELSE
        BP=GBARO ! default (29.92 in Hg)
      ENDIF
C
C  Convert input temperature to Kelvin and compute T0 term
C 
      TK=5.*(TEMHRLY(IH)-32.)/9.+273.15
      T0=647.27-TK
C
C  Compute Beta and PV terms
C
      BETA=(T0/TK)*(C(3)+C(4)*T0+C(5)*T0**3)/(1.0+C(6)*T0)
      PV=BP*C(2)*(10.**BETA)*(100./RH_HRLY(IH))
C
C  Compute absolute humidity. (Relative humidity (%) is converted to a fraction.)
C
      CABSHUM=C(1)/(PV-1.)
C
      RETURN
      END
