C Version @(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD38
C
C  BLOCK DATA Subprogram 38: /EVAL/,/X1/,/X2/,/X3/,/X4/,/X5/
C
C  BLOCK DATA 38 initializes the storage arrays used for by-month calculations.
C
C  Common block array subscripts:
C
C  MONTH(12)     -  MONTH ( MEVAL )
C
C  Common block dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------
C  /EVAL/:
C
C  MEVAL    I    month of evaluation for by-month output, 1=Jan, 2=Feb, etc.
C  MONTH    C*4  4 character abbreviations of the months
C
C  Notes:
C
C  BD 38 was added for MOBILE5v1.
C
C
      CHARACTER*4 MONTH
C
      COMMON /EVAL/ MEVAL,MONTH(12)
      COMMON /X1/ XEVP(4,2),XGREVP(4,9,2),XVGREV(4,2),XGRPD(3,9,2)
      COMMON /X2/ XFTP(3,9,2),XEXH(9,2),XEVAP(9,2),XREFL(9,2),XRUNL(9,2)
      COMMON /X3/ XRLG(9,2),XRSTG(9,2),XRSTL(9,2),XIDLE(3,9,2)
      COMMON /X4/ XVFTP(3,2),XVEXH(2),XVEVAP(2),XVLOSS(2),XVRUNL(2),
     *            XVRSTL(2),XVIDLE(3,2)
      COMMON /X5/ XREG(8,2),XVMT(8,2)
C
C  /EVAL/:
C
      DATA MEVAL/1/
      DATA MONTH/'Jan.','Feb.','Mar.','Apr.','May.','June',
     *           'July','Aug.','Sep.','Oct.','Nov.','Dec.'/
C
C  /X1/: temporary storage for /EVAPGR/ variables
C
      DATA XEVP,XGREVP,XVGREV,XGRPD/142*0.0/
C
C  /X2/: temporary storage for /RESUL1/ variables
C
      DATA XFTP,XEXH,XEVAP,XREFL,XRUNL/126*0.0/
C
C  /X3/: temporary storage for /RESUL2/ variables
C
      DATA XRLG,XRSTG,XRSTL,XIDLE/108*0.0/
C
C  /X4/: temporary storage for /RESUL3/ variables
C
      DATA XVFTP,XVEXH,XVEVAP,XVLOSS,XVRUNL,XVRSTL,XVIDLE/22*0.0/
C
C  /X5/: temporary storage for /VMXCOM/ variables
C
      DATA XREG,XVMT/32*0.0/
C
      END
