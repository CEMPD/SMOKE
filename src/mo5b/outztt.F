C Version "@(#)$Id$ $Source$ $Date$ 

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
C    /LEVBLK/ APPLEV
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
C  (22-April-1994) @ CSC ljn Request 2-420 added LDGT2 (IV=3) to the LEV's
C  October-17-1994 @CSC-pme request 2-472. Checks for APPLEV flag were added
C  to identify the type of messages that need to be printed whether or not
C  the LEV standards are applied to LDGT2. Also it is checked to identify
C  if the ZEV TF percentages should be calculated.
C
      INCLUDE 'BASE10.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'LEVBLK.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'MYRCAL.I'
C
C     (22-April-1994) @ CSC ljn Request 2-420
      DO 20 IV = 1, APPLEV + 1
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
