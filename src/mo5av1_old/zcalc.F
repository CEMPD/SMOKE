C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION ZCALC(X, Y, Z)
C
C  ZCALC perform linear interpolation given two values (X and Y) and
C  the intermediate ratio Z.
C
C  Called by ADJUST.
C
C  Calls none.
C
C  Input on call:
C
C    parameter list: X,Y,Z
C
C  Output on return:
C
C    function: ZCALC
C
C  Notes:
C
C  ZCALC was added for MOBILE5v1.
C
C
C  Assign base year value. Then if that value is non-zero, interpolate.
C  This prevents a model year from suddenly appearing before it exists.
C
      ZCALC = X
      IF (X .GT. 0.) ZCALC = X + (Y - X) * Z
C
      RETURN
      END
