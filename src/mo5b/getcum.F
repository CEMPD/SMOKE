C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETCUM
C
C  GETCUM computes average January 1 fleet cumulative mileages.  GETCUM is
C  called on the first scenario of every run.  If the user enters AMAR, that
C  is the only call.  If the default AMAR is used and therefore ICY dependent
C  HDDV AMAR is calculated for every scenario whose ICY differs from the
C  prior scenario, then GETCUM is correspondingly called to recompute CUMMIL.
C  (Actually, only HDDV CUMMIL need be recalculated - a MOBILE5 fix)
C
C  Called by EFCALX and ONESEC.
C
C  Input on call:
C
C    common blocks:
C    /MAXIMA/ MAXVEH,MAXYRS
C    /MYRCAL/ XMYM
C
C  Output on return:
C
C    common block: /CUMCOM/ CUMMIL
C
C  Notes:
C
C  GETCUM header description was modified for MOBILE4 In-House Version 02 to
C  explain the circumstances under which EFCALX calls GETCUM.
C
C
      INCLUDE 'CUMCOM.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'MYRCAL.I'
C
      DO 10 IV=1,MAXVEH
      CUMMIL(1,IV)=XMYM(1,IV)*100000.
C
      DO 10 JDX=2,MAXYRS
      CUMMIL(JDX,IV)=XMYM(JDX,IV)*100000.+CUMMIL(JDX-1,IV)
   10 CONTINUE
C
      RETURN
      END
