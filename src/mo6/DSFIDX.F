      FUNCTION DSFIDX(ICY,MY,IVTEMP) 
C 
C  DSFIDX returns the sales fraction of gas vehicles based on the 
C  given MY for the given IV and IREJN, using user-input diesel sales 
C  fractions, if appropriate (if DSFLAG=2) 
C 
C  Called by EFCALX and REGMOD. 
C 
C  Calls IGSFPT. 
C 
C  Changes: (Last change first) 
C 
C 29 Oct 01 @EPA-djb    Changed all code to upper case letters. 
C  4 Nov 1998 @DynTel-ZK 2-000 Explicit Typing  
C  Input on call: 
C 
C    argument list: ICY,MY,IVTEMP 
C    common blocks: 
C    /FLAGS1/ DSFLAG 
C    /DSFCOM/ MAXGSF,DSFRAC 
C    /REGISF/ USRGSF,DSFLAG 
C    /VVNAMS/ VVHEAVY, VVLDGAS  (IVTYPE.I) 
C    /VVTYPS/ VTLDDV, VTLDGV, VTMC  (IVTYPE.I) 
C 
C  Output on return: 
C 
C    function: DSFIDX 
C 
C  Notes: 
C 
C  DSFIDX was added in MOBILE4.1 to support user-entered gas/diesel 
C  sales fractions. 
C 
C 29 Oct 01 @EPA-djb    Changed all code to upper case letters. 
C  25 Jan 2000 @EPA-djb Changed from gas fraction to diesel fraction. 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. 
C 
      IMPLICIT NONE  
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'DSFCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'REGISF.I' 
C 
      REAL  DSFIDX 
      INTEGER IVTEMP, JDX, ICY, MY, IGSFPT 
C 
C  Initialize DSFIDX value to indicate no diesel sales 
C 
      DSFIDX = 0.0000 
C 
C  Motorcycles do not have diesels. 
C  Urban/Transit busses are assumed to be all diesel. 
C 
      IF(IVTEMP.GT.14) RETURN 
C 
C  Mapping of vehicle classes to array groupings: 
C 
C     1 : LDV 
C     2 : LDT1 
C     3 : LDT2 
C     4 : LDT3 
C     5 : LDT4 
C     6 : HDV2B 
C     7 : HDV3 
C     8 : HDV4 
C     9 : HDV5 
C    10 : HDV6 
C    11 : HDV7 
C    12 : HDV8A 
C    13 : HDV8B 
C    14 : HDBS (school bus) 
C 
C  Switch to user-entered gas/diesel sales fractions 
C  if user supplied fractions are available. 
C 
      IF(DSFLAG.NE.2) GO TO 10 
        JDX=ICY-MY+1 
        DSFIDX=USRGSF(JDX,IVTEMP) 
      RETURN 
C 
C  Default "built-in" gas/diesel sales fraction case 
C 
   10 DSFIDX=DSFRAC(IGSFPT(MY),IVTEMP) 
      RETURN 
C 
      END 
