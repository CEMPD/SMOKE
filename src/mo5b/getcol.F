C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETCOL(INERR,*)
C
C  GETCOL obtains the cold CO standard parameters.
C
C  Called by PARSEC.
C
C  Calls QUITER.
C
C  Input on call:
C
C    common blocks:
C    /FLAGS1/ PROMPT
C    /IOUCOM/ IOUGEN,IOUASK
C    /YEARS4/ IYEND
C
C  Output on return:
C
C    parameter list: INERR
C    common block: /TEMPC5/ NCC,MYCOL,CCSTD
C    common block: /TEMPC6/ PIVPCT
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  IPIV     I    Phase-In model year pointer MYCOL=1, MYCOL+1=2, etc.
C  LOWCOL   I    lowest (earliest) allowed cold CO standard 1st model year
C  STDMAX   R    maximum allowed cold CO standard
C  STDMIN   R    minimum allowed cold CO standard
C                (i.e., a normalized vmt distribution has been entered)
C
C  Notes:
C
C  GETCOL was added to MOBILE4.1, but is disabled.
C  The CAA COLD CO values are set in BD17.
C
C
      INCLUDE 'FLAGS1.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'TEMPC5.I'
      INCLUDE 'TEMPC6.I'
      INCLUDE 'YEARS4.I'
C
C  Change STDMAX back to 15.0 after validating the version.  This will prevent
C  the user from entering a standard which would produce worse results than no
C  cold CO standard at all (the break even point is 16.6222...).  The code has
C  a range check in it to prevent CSUB from being negative and CMUL from
C  exceeding 1.0, but STDMAX = 15.0 eliminates the need for those checks.
C
      DATA LOWCOL/1980/,STDMIN/2.0/,STDMAX/17.0/
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)
  200 FORMAT(/'&','Enter the number of cold CO standards (1 or 2):')
      READ(IOUGEN,100,ERR=99,END=99) NCC
  100 FORMAT(I1)
      IF(NCC.LT.1.OR.NCC.GT.2) CALL QUITER(0.,NCC,12,INERR)
C
C  Exit if request other than 1 or 2 standards.
C
      IF(NCC.LT.1.OR.NCC.GT.2) GOTO 99
C
      DO 20 ICOLD=1,NCC
      IF(PROMPT.EQ.2) WRITE(IOUASK,210) ICOLD
  210 FORMAT(/' ','Enter cold CO standard parameter record #',I1,':')
      READ(IOUGEN,110,ERR=99,END=99)
     *    MYCOL(ICOLD),(CCSTD(IVLDG,ICOLD),IVLDG=1,3),
     *    ((PIVPCT(IPIV,IVLDG,ICOLD),IPIV=1,3),IVLDG=1,3)
  110 FORMAT(I2,3F5.1,9F6.1)
C
      CALL YRTEST(MYCOL(ICOLD),13,LOWCOL,IYEND,INERR)
      DO 10 IVLDG=1,3
      IF(CCSTD(IVLDG,ICOLD).LT.STDMIN.OR.CCSTD(IVLDG,ICOLD).GT.STDMAX)
     *    CALL QUITER(CCSTD(IVLDG,ICOLD),ICOLD,14,INERR)
   10 CONTINUE
      DO 15 IVLDG=1,3
      DO 15 IPIV=1,3
      IF(PIVPCT(IPIV,IVLDG,ICOLD).LT.0.0 .OR.
     *   PIVPCT(IPIV,IVLDG,ICOLD).GT.100.0)
     *    CALL QUITER(PIVPCT(IPIV,IVLDG,ICOLD),ICOLD,16,INERR)
   15 CONTINUE
   20 CONTINUE
C
      IF(NCC.EQ.1) RETURN
      IF(MYCOL(1).GE.MYCOL(2)) CALL QUITER(0.,0,15,INERR)
C
      RETURN
C
C  Read EOF or illegal NCC => take alternate return 1 => run aborts.
C
   99 CALL QUITER(0.,0,136,INERR)
C
      RETURN1
      END
