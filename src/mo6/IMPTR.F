      FUNCTION IMPTR(MY,IV) 
C 
C  IMPTR selects the emission control technology group pointer into the I/M 
C  credits arrays for LDGV/T.  In the case of TECH 4+, the returned value will 
C  have to be decreased by 3 to properly index the TECH 4+ arrays. 
C 
C  Called by PCLEFT. 
C 
C  Changes: (Last change first) 
C 
C  20 Nov 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: MY,IV 
C    common blocks: 
C    /IMPAR5/ MYGIM2,LBIM4P 
C    /VVTYPS/ VTLDGV  (IVTYPE.I) 
C 
C  Output on return: 
C 
C    function: IMPTR 
C 
C 
C  Notes: 
C 
C  June-15-1993 @ ARC-bk Subtask 244 the common block IMPAR5 
C  was modified, the HDGV I/M credit was expanded to include 
C  an initial value for NOx. 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. 
C 
      IMPLICIT NONE  
C 
      INCLUDE 'IMPAR5.I' 
      INCLUDE 'IVTYPE.I' 
C 
      INTEGER IG, MY, IV, IMPTR 
C 
C  First check if Technology 1-2 case. 
C 
      DO 10 IG=1,3 
      IF(MY.LE.MYGIM2(IG,IV)) GOTO 20 
   10 CONTINUE 
C 
C  Technology 4+ case. 
C  Add 3 to get relative slot (first Tech 4+ slot is 1 + 3 = 4th overall). 
C  When used to index a Tech 4+ array, the 3 will be subtracted first. 
C  Note that the MY to credit mapping differs between LDGV and LDGT, both 
C  indexing the same array, but at different locations for a given MY. 
C 
      IG=IG-1 
      IF(MY.LT.LBIM4P) IG=MY-MYGIM2(IG,IV)+3 
      IF(MY.GE.LBIM4P) IG=LBIM4P-MYGIM2(IG,IV)+3 
      IF(IV.NE.VTLDGV.AND.MY.GE.1988) IG=IG+3 
C 
   20 IMPTR=IG 
      RETURN 
      END 
