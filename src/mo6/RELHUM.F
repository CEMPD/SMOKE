      REAL FUNCTION RELHUM(X) 
C 
C  This function calculates the relative humidity RELHUM (%) given 
C  the absolute humidity AHUM (grain/lb) and temperature X (F). The 
C  relative humidity is calculated using formula RELHUM=AHUM/S*100 
C  where S is the saturation humidity (grain/lb) and is a function 
C  of X (F) given by the arrays XA and SA. 
C 
C 
C Called by HEAT. 
C 
C Calls 
C 
C Changes: 
C 
C  04 Aug 00 @EPA-bg  bug203 removed blanket use of SAVE; applied it only to 
C            appropriate local variables; removed absolute humidity from 
C            parameter list; used common ALUHIN.I to pass ABSHUM 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C  25 Apr 00 AIR Task 02 Added SAVE statement to keep local variables/arrays 
C  15 Mar 00 @DynTel-ddj 1-013, Put a limit on relative humidity, i.e., 100%  
C  21 Jan 00 @DynTel-km 1-013,  New Function 
C 
C 
C Input on call: 
C 
C   argument list : X 
C 
C   common blocks : 
C   /ALUHIN/ ABSHUM 
C 
C Output on return: 
C 
C   function      : RELHUM 
C 
C Local variable / array dictionary : 
C 
C Name            Type         Description 
C -----          ------    ----------------------------------------- 
C X                R        Temperature in degrees Fahrenheit (F) 
C S                R        Saturation humidity (grain/lb) 
C RELHUM           R        Relative humidity (%) 
C XA               R        Temperature array in degrees Fahrenheit (F) 
C SA               R        Saturation humidity array (grain/lb)  
C I                I        Index variable 
C N                I        XA and SA arrays size 
C 
C Notes : 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ALUHIN.I' 
C 
C     Declare parameter list 
C 
      REAL, INTENT(IN)    :: X 
C 
C     Declare local variables 
C 
      REAL S 
      REAL, SAVE, DIMENSION(23) :: XA=( 
     *                      /10.,15.,20.,25.,30.,35.,40.,45., 
     *                       50.,55.,60.,65.,70.,75.,80.,85., 
     *                       90.,95.,100.,105.,110.,115.,120./) 
      REAL, SAVE, DIMENSION(23) :: SA=( 
     *                      /10.,12.,16.,20.,25.,30.,37.,44., 
     *                       53.,65.,79.,93.,111.,132.,157.,184., 
     *                       218.,257.,290.,338.,392.,456.,528./) 
C 
      INTEGER, SAVE :: N=23 
      INTEGER I 
C 
C 
C 
      IF (X.GE.XA(1) .AND. X.LT.XA(N)) THEN 
         I=1 
         DO WHILE (X.GE.XA(I)) ! Find the nearest neighbors of X 
            I=I+1              ! in the XA array 
         ENDDO 
         S=SA(I-1)+(X-XA(I-1))*(SA(I)-SA(I-1))/(XA(I)-XA(I-1)) 
      ELSE                ! X is outside of the values of the XA array 
         IF (X.LT.XA(1)) S=SA(1)-(XA(1)-X)*(SA(2)-SA(1))/(XA(2)-XA(1)) 
         IF (X.GE.XA(N)) S=SA(N)+(X-XA(N))*(SA(N)-SA(N-1))/(XA(N) 
     &      -XA(N-1)) 
      ENDIF 
C 
      RELHUM=ABSHUM/S*100 
      RELHUM = MIN(RELHUM, 100.0) 
C 
      RETURN 
      END 
