      SUBROUTINE EFFGRP(MY) 
C 
C  EFFGRP builds the tampering effects group sizes array.  The tampering 
C  rates from DISATP are modified  to account for technology vs tampering 
C  types, overlap factors, anti-tampering program effectiveness factors and 
C  technology group sizes.  The output rates are divided into effects groups 
C  corresponding to the significant emission impact tampering groups. 
C 
C  Called by TAMPER. 
C 
C  Calls ADJ12, ATPEFF, A8TO11, C1TO7, C8TO12, ITAMPT, SETEGS and TECH0. 
C 
C  Changes: (Last change first) 
C 
C  04 Dec 00 @EPA-bag Task X6, removed MY, JDX, and IDX from MYCODE.I; added 
C            model year, MYEFFGRP, to parameter list; ITAMPT also has model 
C            year passed in the parameter lsit 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C  18 Sep 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    common blocks: 
C    /EGSCAL/ TGS 
C    /LOOKUP/ IVTAM 
C    /MYCODE/ LDXSY,LMYRVT 
C    /SIZCAL/ CSIZE 
C 
C  Output on return: 
C 
C    common blocks: 
C    /EGSCAL/ EGS 
C    /LOOKUP/ IQG,IPG,JPGD,IGCSF 
C    /MYCODE/ IAY 
C    /SIZCAL/ TGSUSE,CSAE 
C 
C  Local array subscripts: 
C 
C  CSF1ST(7,2,3) -  CSF1ST ( IC7, IPG, IGCSF ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------ 
C  CSF1ST   R    category size factor for 1st calculation of 
C                CSIZE, expressed as <overlap case's CSF> = 
C                <o.c. freq> / <base case freq>, where 
C                <base case> = AorE for IC7=1-3, CATS for 
C                IC7=4-7 
C 
C 
C  Notes: 
C 
C  For MOBILE4.1 EFFGRP was modified to expand CSF1ST to 3 MYGs & 
C  updates its rates. 
C  March 9, 1993 IY was added to the argument list to guarantee 
C  the 1990 inventory emission factors remain the same in newer 
C  versions of MOBILE as in MOBILE5.  And ICY is needed. 
C  Sep-20-1994 @ CSC-ked (for bsg) Request 438 Fix 1990 and earlier calendar 
C  year removed IY from the argument list, and the calculation of ICY is no 
C  longer needed 
C 
      IMPLICIT NONE  
C 
      INCLUDE 'EGSCAL.I' 
      INCLUDE 'LOOKUP.I' 
      INCLUDE 'MYCODE.I' 
      INCLUDE 'SIZCAL.I' 
C 
      INTEGER, INTENT(IN) :: MY 
      INTEGER IG1, ITAMPT, IC, IEG  
      REAL, DIMENSION(7,2,3) :: CSF1ST=RESHAPE((/ 
C 
C  Continuation code is IPG = H (HC & CO) & N (NOx) 
C  MYG 1 applies to all MY for HDGV; the indicated ranges are for LDGV/T. 
C                            CSIZE 
C         1      2      3      4      5      6      7 
C  MYG 1 (pre-1981) 
     H  .0708, .1420, .1538, .5529, .0000, .1513, .2500, 
     N  .0948, .1209, .0033, .3059, .0039, .3843, .0118, 
C  MYG 2 (1981-83) 
     H  .1421, .0519, .0364, .2840, .0406, .1917, .1469, 
     N  .0698, .0233, .0116, .2727, .0260, .2338, .0130, 
C  MYG 3 (1984+) 
     H  .1739, .0750, .0938, .1011, .0000, .0470, .0556, 
     N  .1212, .0000, .0000, .0000, .0000, .0588, .0000/),(/7,2,3/)) 
C 
C  Zero out tampering effects group sizes (EGS) array. 
C 
      DO 20 IPG=1,2 
      DO 10 IEG=1,7 
      EGS(IEG,IPG)=0.0 
   10 CONTINUE 
   20 CONTINUE 
C 
C  Add in each equipment (technology) group's before and after anti-tampering 
C  program (ATP) start year contribution to each EGS category for each 
C  pollutant group. 
C 
      DO 80 IPG=1,2 
C 
C  Set ATR pointer for 1st 3 sizes: HC/CO => use AIRS, NOX => use EGRS. 
C 
      IF(IPG.EQ.1) JPGD=1 
      IF(IPG.EQ.2) JPGD=5 
C 
      DO 70 IQG=1,6 
C 
C  Lookup technology group size to be used on this pass. 
C 
      IG1=ITAMPT(MY,1) 
      TGSUSE=TGS(IG1,IQG,IVTAM) 
C 
C  Check how technology type & size affect tampering rates. 
C 
      CALL TECH0(*70) 
C 
C  HC & CO effects category sizes are identical => only 2 pollutant passes. 
C 
      DO 65 IAY=1,2 
C 
      IGCSF=ITAMPT(MY,5) 
C 
      CALL C1TO7(CSF1ST(1,IPG,IGCSF)) 
      CALL C8TO12 
C 
C  The first calculation pass can produce negative sizes for the no 
C  overlap (8-11) and the no effects (12) categories.  Correct these 
C  by calling A8TO11 and ADJ12, respectively. 
C 
      DO 30 IC=8,11 
      IF(CSIZE(IC,IAY).LT.0) CALL A8TO11(IC) 
   30 CONTINUE 
      IF(CSIZE(12,IAY).LT.0) CALL ADJ12(CSF1ST) 
C 
C  "Subsequent" (to ATP start year) pass + have split mileage 
C  case => "subsequent" category size = overall size - "previous" 
C  size i.e., in the split mileage case, the program computes the 
C  previous to ATP start year and overall values up to this 
C  point.  After this next loop assigns the subsequent values as 
C  residuals, the program works directly with the previous and 
C  Check if ATP covers this MY.  Yes => apply effectiveness factors to the 
C  category sizes. 
C 
      IF(LMYRVT.EQ.2) CALL ATPEFF(*60) 
C 
C  Otherwise just shift IAY's sizes to CSAE without modification. 
C 
      DO 55 IC=1,12 
      CSAE(IC)=CSIZE(IC,IAY) 
   55 CONTINUE 
C 
C  Set the effects group sizes. 
C 
   60 CALL SETEGS 
C 
   65 CONTINUE 
   70 CONTINUE 
   80 CONTINUE 
C 
      RETURN 
      END 
