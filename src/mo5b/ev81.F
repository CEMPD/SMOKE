C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION EV81(KEYEVP,IFDS,IVGAS,RIVAL)
C
C  EV81 calculates a hot soak or diurnal rate for the 1981+ LDGV and
C  LDGT1/2 emission standards, given the fuel delivery system and the
C  fuel RVP or Wade Index, respectively.
C
C  Called by EVMAIN, EVSTD and EV85.
C
C  Input on call:
C
C    parameter list: KEYEVP,IFDS,IVGAS,RIVAL
C
C    common blocks:
C    /EVADU1/ DUEQ
C    /EVAHS1/ HSEQ
C    /EVAPHD/ IVVMAP
C
C  Output on return:
C
C    function: EV81
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  IVV      I    vehicle (IVV=1) or truck (IVV=2) classification switch
C
C
C  Notes:
C
C  Modified for MOBILE4.1 to handle pressure/purge fail effects.
C  Modified for MOBILE5v13 to reflect changes in diurnal.
C     The changes account for the difference
C     between a 1 hour measurement for diurnal and an 8 hour
C     measurement of diurnal emissions.
C  Sep-20-1994 @ CSC-ked (for bsg) Request 438 Fix 1990 and earlier calendar
C  year removed IY and ICY from the parameter list
C
C
      INCLUDE 'EVADU1.I'
      INCLUDE 'EVAHS1.I'
      INCLUDE 'EVAPHD.I'
C
C  Map for car or truck.
C
      IVV=IVVMAP(IVGAS)
C
      IF(KEYEVP.EQ.2) GOTO 10
C
C  Hot soak
C
      IF(RIVAL.LT.9.0) EV81=HSEQ(IFDS,IVV,1)
     *                     +HSEQ(IFDS,IVV,2)*RIVAL
      IF(RIVAL.GE.9.0) EV81=HSEQ(IFDS,IVV,3)
     *                     +HSEQ(IFDS,IVV,4)*RIVAL
     *                     +HSEQ(IFDS,IVV,5)*RIVAL**2
      EV81=AMAX1(EV81,0.0)
      RETURN
C
C  Diurnal
C  Use LDGV parameters for low fuel injected RIVAL values.
C
   10 IF(IFDS.GE.2.AND.RIVAL.LT.0.9) IVV=1
C
      IF(IFDS.EQ.1.AND.RIVAL.LT.1.0) THEN
      IF(IVV.EQ.1) EV81=1.9488*RIVAL-0.2488
      IF(IVV.GT.1) EV81=1.9488*RIVAL+0.1412
      GOTO 20
      ENDIF
C
      EV81=DUEQ(IFDS,IVV,1)
     *    +DUEQ(IFDS,IVV,2)*RIVAL
     *    +DUEQ(IFDS,IVV,3)*RIVAL**2
C
   20 EV81=AMAX1(EV81,0.0)
      RETURN
      END
