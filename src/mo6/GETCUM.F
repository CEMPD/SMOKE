      SUBROUTINE GETCUM 
C 
C  GETCUM computes average January 1 fleet cumulative mileages.  GETCUM is 
C  called on the first scenario of every run.  If the user enters AMAR, that 
C  is the only call.  If the default AMAR is used and therefore ICY dependent 
C  HDDV AMAR is calculated for every scenario whose ICY differs from the 
C  prior scenario, then GETCUM is correspondingly called to recompute CUMMIL. 
C  (Actually, only HDDV CUMMIL need be recalculated - a MOBILE5 fix) 
C 
C  Called by CHKINP and EFCALX. 
C 
C  Changes: (Last change first) 
C 
C  08 Jun 00 AIR Task 03: Removed nonblock DO constructs. 
C  07 Oct 98 @DynTel-ZK 2-000 Explicit Typing 
C  Input on call: 
C 
C    Constants:  
C    MAXVEH from include file IVTYPE.I 
C 
C    common blocks: 
C    /MAXIMA/ MAXYRS 
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
C  11 Jun 97 DynTel@MLA 2-622  Moved constant MAXVEH form include file 
C            MAXIMA to include file IVTYPE. Added declaration for IVTYPE. 
      IMPLICIT NONE  
      INCLUDE 'CUMCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'MYRCAL.I' 
      INTEGER IV, JDX 
       
      DO IV=1,MAXVEH 
        CUMMIL(1,IV)=XMYM(1,IV)*100000. 
C 
        DO JDX=2,MAXYRS 
          CUMMIL(JDX,IV)=XMYM(JDX,IV)*100000.+CUMMIL(JDX-1,IV) 
        END DO 
C 
      END DO 
C 
      RETURN 
      END 
