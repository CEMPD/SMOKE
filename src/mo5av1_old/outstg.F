C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTSTG
C
C  OUTSTG prints the Stage II program parameters.
C
C  Called by OUTHD3 and OUTHD4.
C
C  Input on call:
C
C    common blocks:
C    /IOUCOM/ IOUREP
C    /RLCOM2/ IS2SY,NPHASE,S2EFF(4)
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  S2EFF1   R    for output of percent efficiency for LDGV and LDGT
C  S2EFF4   R    for output of percent efficiency for HDGV
C
C  Notes:
C
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /RLCOM2/ IOBMY,IVOB(4),IS2SY,NPHASE,S2EFF(4)
C
C  Calculate percent efficency
C
      S2EFF1=S2EFF(1)*100.0
      S2EFF4=S2EFF(4)*100.0
C
      WRITE(IOUREP,200)
  200 FORMAT('0','Stage II program selected:'/)
      WRITE(IOUREP,250) IS2SY,NPHASE,S2EFF1,S2EFF4
  250 FORMAT('0','    Start year (January 1):             ',I4/
     *       ' ','    Phase-in period (yrs.):                ',I1/
     *       ' ','    Percent Efficiency for LDGV & LDGT: ',F3.0,'%'/
     *       ' ','    Percent Efficiency for HDGV:        ',F3.0,'%'/)
C
      RETURN
      END
