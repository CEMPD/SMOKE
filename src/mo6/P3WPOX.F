      SUBROUTINE P3WPOX(MYP3W,IV,P3W,POX) 
C 
C  Computes the number of catalyst and non-catalyst vehicles, 
C 
C  Called by BIGRVP, EMIRAT and FUEL 
C 
C  Calls ITAMPT 
C 
C  Changes: (Last change first) 
C 
C  04 Dec 00 @EPA-bag Task X6, removed MY, JDX, and IDX from MYCODE.I; removed 
C            MYCODE.I; model year now passed through parameter lists 
C  20 Jun 00 @EPA-bg  2-000 Added INTENT qualifier to 
C                           parameter variable declarations 
C  16 Dec 98 @ DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C  argument list: IV 
C    common blocks: 
C    /EGSCAL/ TGS 
C 
C  Output on call: 
C 
C  argument list: P3W, POX 
C    common blocks: 
C    /LOOKUP/ IVTAM,IQG 
C    /VVNAMS/ VVGAS  (IVTYPE.I) 
C    /VVTYPS/ VTGASBUS  (IVTYPE.I) 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  IG1      I    Return value of ITAMPT 
C  IQG      I    Pointer in TGS array 
C  IVTAM    I    See parameter dictionary 
C 
C 
C  Notes:  
C 
C  Sep-22-1994 @ CSC-bsg Request 426  P3WPOX was created to modularize 
C  BIGRVP, EMIRAT and FUEL.  BIGRVP is needed to calculate the NOx RFG 
C  effects based on percentage of vehicles equipped with 3-way catalysts. 
C 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. 
C 
      IMPLICIT NONE 
      INCLUDE 'EGSCAL.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'LOOKUP.I' 
      INTEGER, INTENT(IN)  :: MYP3W 
      INTEGER, INTENT(IN)  :: IV 
      REAL,    INTENT(OUT) :: P3W 
      REAL,    INTENT(OUT) :: POX 
      INTEGER IG1, ITAMPT 
C 
      P3W=0.0 
      POX=0.0 
C 
      IF(VVGAS(IV).EQ.1) THEN 
        IVTAM=IV 
        IF (IV.EQ.VTGASBUS) THEN 
          IVTAM = 14 
        ENDIF 
        DO 10 IQG=5,6 
          IG1=ITAMPT(MYP3W,1) 
          P3W=P3W+TGS(IG1,IQG,IVTAM) 
   10   CONTINUE 
C 
        DO 20 IQG=2,3 
          IG1=ITAMPT(MYP3W,1) 
          POX=POX+TGS(IG1,IQG,IVTAM) 
   20   CONTINUE 
        POX=POX-P3W 
      ENDIF 
C 
      RETURN 
      END 
