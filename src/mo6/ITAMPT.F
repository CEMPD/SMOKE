      FUNCTION ITAMPT(MY,ICASE) 
C 
C  ITAMPT returns the pointer to the technology groups cell location holding 
C  the desired data point, given MY.  ICASE identifies the calling CSECT 
C  location and hence the array to be indexed and thus the my groups array 
C  to search for the pointer.   
C 
C  Called by DISCAL, EFFGRP, EMIRAT, EVPEME, FAIL and P3WPOX. 
C 
C  Input on call: 
C 
C  Changes: (Last change first) 
C 
C  04 Dec 00 @EPA-bag Task X6, removed MY, JDX, and IDX from MYCODE.I; removed 
C            MYCODE.I; model year now passed through parameter lists 
C  23 Nov98 @DynTel-ZK 2-000 Explicit Typing 
C 
C    argument list: ICASE 
C    common blocks: 
C    /LOOKUP/ IVTAM,IQG,IHG 
C    /MYCODE/ MY 
C    /MYUB1/  MYTGS,MYEGR,MYVTS 
C    /MYUB2/  MYCCEI 
C    /TAMEQ2/ MYGTAM 
C 
C  Output on return: 
C 
C    function: ITAMPT 
C 
C  Notes: 
C 
C  Differences from the other pointer functions: 
C    (1) 5 different my group arrays are handled by this 1 subroutine, instead 
C        of having a separate function for each case. 
C    (2) ITAMPT assumes no default value, as the others do.  I.e. the typical 
C        code: 
C                C 
C                C..Set default 
C                C 
C                      IG=1 
C                      CALL QUITER(0.,6,61,MY) 
C        is not inserted after each search loop.  All the my arrays have 2051, 
C        the highest ICY allowed by subroutine input, as the last non-dummy 
C        upper bound.  Hence, an in-range match has to occur. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'LOOKUP.I' 
      INCLUDE 'MYCODE.I' 
      INCLUDE 'MYUB1.I' 
      INCLUDE 'MYUB2.I' 
      INCLUDE 'TAMEQ2.I' 
C 
      INTEGER, INTENT(IN)  :: MY 
      INTEGER ICASE, IG, ITAMPT 
C 
C  Branch to array pointer loop being asked for by ICASE. 
C 
      GOTO(10,20,30,40,50),ICASE 
C 
C  Case 1: find MY's cell in TGS. 
C 
  10  DO 15 IG=1,20 
      IF(MY.LE.MYTGS(IG,IQG,IVTAM)) GOTO 90 
  15  CONTINUE 
C 
C  Case 2: find MY's cell in EMEGR. 
C 
  20  DO 25 IG=1,4 
      IF(MY.LE.MYEGR(IG,IVTAM)) GOTO 90 
  25  CONTINUE 
C 
C  Case 3: find MY's cell in EVPTGS. 
C 
  30  DO 35 IG=1,3 
      IF(MY.LE.MYVTS(IG,IHG,IVTAM)) GOTO 90 
  35  CONTINUE 
C 
C  Case 4: find MY's cell in CCEMI. 
C 
  40  DO 45 IG=1,14 
      IF(MY.LE.MYCCEI(IG,IVTAM)) GOTO 90 
  45  CONTINUE 
C 
C  Case 6: find MY's cell in TAMZML/DR/A50 or CSF1ST 
C 
   50 DO 55 IG=1,3 
      IF(MY.LE.MYGTAM(IG,IVTAM))  GOTO 90 
   55 CONTINUE 
C 
C  Assign cell number. 
C 
   90 ITAMPT=IG 
C 
      RETURN 
      END 
