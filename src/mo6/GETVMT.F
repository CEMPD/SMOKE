      SUBROUTINE GETVMT(NSCEN,INERR,RC) 
C 
C  GETVMT gets the fleet vehicle miles travelled (vmt) mix by vehicle type. 
C 
C  Called by DSPTCH. 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C 
C  23 Aug oo @EPA-djb bug231 Changed fixed format to free format. 
C  05 Jul 00 @EPA-ddj bug#169 Added bounds checks for vmt's 
C  24 Apr 00 AIR Task 02: Removed spurious INPREC from parameter list. 
C  11 Apr 00 AIR Task 01: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C            so that Lahey Standard for keyboard/monitor I/O is followed. 
C  07 Apr 00 @EPA-djb Renamed arrays. 
C  29 Oct 98 @DynTel-ZK 2-000 Explicit Typing 
C  Input on call: 
C 
C    Constants:  
C    MAXVEH from include file IVTYPE.I 
C 
C    ARGUMENT LIST: INERR 
C    common blocks: 
C    /IOUCOM/ IOUGEN,IOUOUT 
C 
C  Output on return: 
C 
C    argument list: INERR,NSCEN,RC 
C    common block: /FLAGS1/ GLVMFLG,VMFLAG 
C                  /VMXCOM/ GLVMTMIX,VMTUSR 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  IVC      I    composite (gas/diesel) vehicle classes (16) 
C  VMTCHK   R    vmt mix check: sum of the vmt fractions = 1.0 +/- 1.E-6 
C                (i.e., a normalized vmt distribution has been entered) 
C 
C  Notes: 
C 
C 11 Jun 97 DynTel@MLA 2-622  Altered code to read Mobile6 transition 
C           input formats. 28 VMT fractions are required in Mobile6. 
C           Added IVTYPE.I for Mobile6. 
C 
C  Oct-30-97 @ DynTel-gjr 2-645 Eliminating option for interactive input. 
C 
                     
      IMPLICIT NONE 
C 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'VMXCOM.I' 
C 
      INTEGER INERR, NSCEN, RC, IVC 
      REAL VMTCHK 
C 
C  Initialize local variables 
C 
      RC = 0 
      VMFLAG = 2 
      IF(NSCEN.EQ.0) THEN 
        GLVMFLG = 2 
      ENDIF 
C 
C  User supplied vehicle mile traveled (VMT) by vehicle class. 
C  VMT is read by composite (gas/diesel) vehicle class (IVC). 
C  Diesel fractions are used to split these composite values. 
C 
C  Combined diesel/gas sales vehicle class groupings are: 
C 
C         IVC 
C        ----- 
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
C    15 : HDBT (transit bus) 
C    16 : MC 
C 
      READ(IOUGEN,*,END=97,ERR=98) (VMTUSR(IVC),IVC=1,16) 
C 
C  Check user input (sum must be 1.0) 
C 
      VMTCHK=0.0 
      DO IVC=1,16 
        VMTCHK=VMTCHK+VMTUSR(IVC) 
        IF(VMTUSR(IVC).LT.0.0.OR.VMTUSR(IVC).GT.1.) 
     *      CALL QUITER(VMTUSR(IVC),0,17,INERR) 
        IF(NSCEN.EQ.0) THEN 
          GLVMTMIX(IVC) = VMTUSR(IVC) 
        ENDIF 
      END DO !IVC 
C 
      IF(VMTCHK+1.E-6.LT.1.0 .OR. VMTCHK-1.E-6.GT.1.0) THEN 
        CALL QUITER(VMTCHK,0,18,INERR) 
      ELSE 
        CALL QUITER(0.0,0,615,INERR) 
      END IF 
C 
      RETURN 
C 
C  Error or EOF on any attempted read. 
C 
   97 RC = -2 
   98 RC = RC + 1 
C 
      CALL QUITER(0.,0,127,INERR) 
      RETURN 
      END 
