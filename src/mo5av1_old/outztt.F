C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTZTT(ICY)
C
C  OUTZTT Calculates Zero Emission Vehicle (ZEV)'s percent of
C  VMTAGE for LDGV's LDGT1's and LDGT2's
C
C  Called by OUTPUT
C
C  Calls PCTLEV
C
C    Input on call:
C    parameter list: ICY
C
C    common blocks:
C    /MAXIMA/ MAXYRS
C    /MYRCAL/ TF
C
C    Output on call:
C    common blocks:
C    /BASE10/ ZEVTF, ZEVTOT
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  IDX      I    see top of program
C  IV       I    see top of program
C  MY       I    see top of program
C
C
C  Notes:
C
C   OUTZTT was added for MOBILE5v3
C   OUTZTT was modified for MOBILE5v13 for calling PCTLEV
C
C
      COMMON /BASE10/ ZEVTF(25,3),ZEVTOT(3)
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /MYRCAL/ XMYM(25,8),JANMYR(25,8),TF(25,8),TFMYM(25,8)
C
C     IF(ICY.LE.1999) IDX = ICY - MAXYRS
C
      DO 20 IV = 1, 2
C
        ZEVTOT(IV) = 0.0
C
        DO 10 IDX = 1,MAXYRS
C
C          get the Model Year MY
C
           MY=(ICY+IDX)-MAXYRS
C
C          add up the ZEV's TF
C
           ZEVTF(IDX,IV) = 0.0
C
           IF(MY.GE.1994) THEN
             ZEVTF(IDX,IV) = TF(IDX,IV) * PCTLEV(MY,IV,8)
             ZEVTOT(IV) = ZEVTOT(IV) + ZEVTF(IDX,IV)
           ENDIF
  10    CONTINUE
C
C       after totaling convert the total TF to a percentage
C
        ZEVTOT(IV) = ZEVTOT(IV) * 100
  20  CONTINUE
C
      RETURN
      END
