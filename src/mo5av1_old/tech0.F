C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE TECH0(*)
C
C  TECH0 determines whether or not to compute the effects groups for the
C  passed in TGSUSE/IQG/IPG values.  If the effects loop is to be completed,
C  TECH0 decides what subset of BTR is to be used.
C
C  Called by EFFGRP.
C
C  Input on call:
C
C    common blocks:
C    /LOOKUP/ IQG,IPG
C    /SIZCAL/ TGSUSE
C    /TAMEQ4/ BTR
C
C  Output on return:
C
C    common blocks:
C    /MYCODE/ IAY
C    /SIZCAL/ ATR
C
C  Notes:
C
C  MOBILE4 excludes evaporative HC ( ID(1-7) is now IDBAG(1-5) ).
C
C  None.
C
C
      COMMON /LOOKUP/ IVTAM,IQG,IPG,JPGD,IHG,IGCSF
      COMMON /MYCODE/ MY,IDX,JDX,LDXSY,LMYRVT,IAY,IMDXSY,IMKINK
      COMMON /SIZCAL/ ATR(5,2),TGSUSE,CSIZE(12,2),CSAE(12)
      COMMON /TAMEQ4/ BTR(9,2)
C
C  If technology size is 0, branch to the end of the IQG loop.  If no part
C  of the fleet is so equipped, then no tampering effects based on that
C  equipment type can occur.
C
      IF(TGSUSE.EQ.0) RETURN1
C
C  If the technology type does not affect the pollutant group, branch to end
C  of IQG loop.  Tampering vehicles so equipped will not add to the level of
C  the pollutant(s) being emitted.
C
      IF(IQG.GE.4.AND.IPG.EQ.1) RETURN1
      IF(IQG.LE.3.AND.IPG.EQ.2) RETURN1
C
C  Technology size is positive and tampering it affects the pollutant group.
C  The effects groups will be calculated.  First assign the base rates.
C
      DO 20 IAY=1,2
      DO 10 IDBAG=1,5
      ATR(IDBAG,IAY)=BTR(IDBAG,IAY)
   10 CONTINUE
C
C  Zero out rate(s) for case(s) where the type tampering cannot occur, given
C  the technology type.
C
      IF(IQG.EQ.3) ATR(1,IAY)=0.0
      IF(IQG.EQ.6) ATR(5,IAY)=0.0
      IF(IQG.NE.1.AND.IQG.NE.4) GOTO 20
      ATR(2,IAY)=0.0
      ATR(3,IAY)=0.0
      ATR(4,IAY)=0.0
   20 CONTINUE
C
      RETURN
      END
