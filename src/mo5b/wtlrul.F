C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE WTLRUL(INERR,*)
C
C  WTLRUL reads in the 6 user-supplied %VMT trip length values for
C  either each scenario (SPDFLG=3) OR ALL scenarios (SPDFLG=4)
C  WTLRUL previously read in 3 weighting values (never used) for TLVMT
C
C  Called by ONESEC and PARSEC.
C
C  Calls QUITER.
C
C  Input on call:
C
C    parameter list: INERR
C    common blocks:
C    /FLAGS1/ PROMPT,SPDFLG
C    /IOUCOM/ IOUGEN,IOUASK
C    /RUNLS1/ MAXTL,TLVMTU
C
C  Output on return:
C
C    common blocks:
C    /RUNLS1/ TLVMTU
C
C  Local array subscripts
C
C  TLVMTU(6)  -  TLVMTU( ITL )
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  TLSUM    R    sum of the trip length VMT fractions
C
C  Notes:
C
C  WTLRUL was activated beginning with Version D02 of MOBILE4.1
C
C  WTLRUL is not used in Release or In-House MOBILE4.
C
C
      INCLUDE 'FLAGS1.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'RUNLS1.I'
C
C  Read in alternate %VMT trip length fractions.
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)
  200 FORMAT(' ','Enter the trip length VMT fractions ( 6(1X,F4.1)):')
      READ(IOUGEN,100,END=99,ERR=99) TLVMTU
  100 FORMAT(6(1X,F4.1))
C
C  Verify %VMT trip length fraction values and that they sum to 1.0.
C
      TLSUM=0.0
      DO 20 ITL=1,MAXTL
      IF(TLVMTU(ITL).LT.0.0.OR.TLVMTU(ITL).GT.100.0)
     *    CALL QUITER(TLVMTU(ITL),0,93,INERR)
      IF(TLVMTU(ITL).LT.0.0.OR.TLVMTU(ITL).GT.100.0)
     *    TLVMTU(ITL)=0.0
      TLVMTU(ITL)=TLVMTU(ITL)/100.0
      TLSUM=TLSUM+TLVMTU(ITL)
   20 CONTINUE
      IF(TLSUM+1.E-6.LT.1..OR.TLSUM-1.E-6.GT.1.)
     *     CALL QUITER(TLSUM,0,94,INERR)
      IF(TLSUM.LE.0.0) GO TO 98
C
C  Normalize TLVMTU
C
      DO 30 ITL=1,MAXTL
      TLVMTU(ITL)=TLVMTU(ITL)/TLSUM
   30 CONTINUE
   98 RETURN
   99 RETURN 1
      END
