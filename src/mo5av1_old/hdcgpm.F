C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE HDCGPM
C
C  HDCGPM converts user entered new emission rates from grams per brake-
C  horsepower units to grams per mile units.
C
C  Called by GETERQ.
C
C  Calls HDCPIC.
C
C  Input on call:
C
C    common block: /BASEQ3/ NEWVEH,NEWMYF,NEWMYL,ZMLNEW,DRNEW
C
C  Output on return:
C
C    common block: /BASEQ3/ ZMLNEW,DRNEW
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  CONV     R    conversion factor to be used for a given MY and IV
C  DRSUM    R    sum of the converted deterioration rates for the model years
C                in the range specified by the user
C  RANGE    R    number of model years in the my range specified by the user
C                for which the new emission rate is to be used.
C  ZMLSUM   R    sum of the converted zero mile levels for the model years
C                in the range specified by the user
C
C  Notes:
C
C  None.
C
C
      COMMON /BASEQ3/ NEWVEH,NEWMYF,NEWMYL,ZMLNEW,DRNEW,KINKS
C
C  Apply each model year's conversion factor to the user's emission rate
C  and add the result to the running sums.
C
      ZMLSUM=0.0
      DRSUM=0.0
      DO 10 MY=NEWMYF,NEWMYL
      CONV=HDCPIC(MY,NEWVEH)
      ZMLSUM=ZMLSUM+(CONV*ZMLNEW)
      DRSUM=DRSUM+(CONV*DRNEW)
  10  CONTINUE
C
C  Divide the sums by the number of my's in the range to get the average
C  emission rate.
C
      RANGE=NEWMYL-NEWMYF+1
      ZMLNEW=ZMLSUM/RANGE
      DRNEW=DRSUM/RANGE
C
C  Alternative, less "user friendly", but more efficient method:
C  Compute the average conversion rate for the interval and multiply
C  times the entered rates.
C
      RETURN
      END
