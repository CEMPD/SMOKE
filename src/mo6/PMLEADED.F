      SUBROUTINE PMLEADED(MY,ICY,IV,PLNOCT,PUNOCT,PLYSCT,PUYSCT,INERR) 
C 
C  LEADED calculates Lead Particulate Emission Factors 
C 
C Changes: 
C 
C 15 Nov 01 AIR Task 27: Added INERR support due to PMFLECON update. 
C 12 Feb 01 AIR Task 11: Modified to 28 vehicle classes and 
C                        now uses MOBILE6 fuel economy values. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Called by PMLEAD. 
C 
C  Calls IPBFUN, PMFLECON 
C 
C  Input on call: 
C    parameter list: MY, ICY, IV 
C 
C    common blocks: See Includes below 
C 
C  Output on return: 
C     PLNOCT,PUNOCT,PLYSCT,PUYSCT 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------------- 
C  FECONM   R    Vehicle fuel economy, miles per gallon, MOBILE6 based 
C  ICY      I    Calendar Year 
C  IPBPTR   I    Index for fuel lead content by calendar year 
C  IV       I    Vehicle type 
C  JDX      I    Model year window index (age) 
C  MY       I    Model year 
C  PBL      R    Lead content of leaded fuel 
C  PBNL     R    Lead content of unleaded fuel 
C  PLNOCT   R    Particulate emission factor for non-catalyst vehicles using leaded fuel 
C  PLYSCT   R    Particulate emission factor for catalyst equipped vehicles using leaded fuel 
C  PUNOCT   R    Particulate emission factor for non-catalyst vehicles using unleaded fuel 
C  PUYSCT   R    Particulate emission factor for catalyst equipped vehicles using unleaded fuel 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'DFSBLK.I'  ! ABURN 
      INCLUDE 'IVTYPE.I'  ! VTMC 
      INCLUDE 'LEAD.I'    ! IPBPTR 
      INCLUDE 'PSIZE.I'   ! PSL, PSNL, PSNLCT 
C 
C  Declare external functions. 
C 
      INTEGER, EXTERNAL :: IPBFUN 
      REAL,    EXTERNAL :: PMFLECON 
C 
C  Declare parameter list. 
C 
      INTEGER ICY 
      INTEGER INERR 
      INTEGER IV 
      INTEGER MY 
      REAL    PLNOCT 
      REAL    PLYSCT 
      REAL    PUNOCT 
      REAL    PUYSCT 
C 
C  Declare local variables/arrays. 
C 
      REAL    FECONM 
      REAL    PBL 
      REAL    PBNL 
C 
C  Initialize 
C 
      PLNOCT=0.0 
      PUNOCT=0.0 
      PLYSCT=0.0 
      PUYSCT=0.0 
C 
C     function to find the pointer index for Lead content of Gas(g/m) 
C     by calendar year. 
C 
      IPBPTR = IPBFUN(ICY) 
C 
C     get the lead content for leaded and unleaded fuel: 
C 
      PBL    = PB(IPBPTR,1) 
      PBNL   = PB(IPBPTR,2) 
C 
C  Jun-24-1994 @ CSC-bsg Request 430 changed MC calculation 
C 
      IF(IV.EQ.VTMC.AND.MY.LE.1978) THEN 
        PLNOCT=(0.49*0.100+0.51*0.014)*PSL*PBL*0.5 
        RETURN 
      ELSEIF(IV.EQ.VTMC.AND.MY.GT.1978) THEN 
        PLNOCT=0.014*PSL*PBL*0.5 
        RETURN 
      ENDIF 
C 
C     for each Vehicle Type get the corresponding fuel economy 
C     (Use MOBILE6 values instead of PART5 values.) 
C 
      FECONM=PMFLECON(MY,IV,INERR) 
C 
C     Partical weight ratios: 
C     1.557  = PbClBr / Pb 
C 
C     determine the lead particulate emission factor 
C     for non-catalyst vehicles using leaded fuel 
C 
      PLNOCT = ((ABURN(1)*1.557)/FECONM) * PBL * PSL 
C 
C     determine the lead particulate emission factor 
C     for non-catalyst vehicles using unleaded fuel 
C 
      PUNOCT = ((ABURN(1)*1.557)/FECONM) * PBNL * PSNL 
C 
C  Jul-22-1994 @ CSC-bsg Fixed mistake in derivation of catalyst lead EF 
C     determine the lead particulate emission factor 
C     for catalyst equipped vehicles using leaded fuel 
C 
      PLYSCT = ((ABURN(1)*1.557)/FECONM) * PBL * PSL 
C 
C  Jul-22-1994 @ CSC-bsg Fixed mistake in derivation of catalyst lead EF 
C     determine the lead particulate emission factor 
C     for catalyst equipped vehicles using unleaded fuel 
C 
      IF(MY .LT. 1981) THEN 
        PUYSCT = ((ABURN(2)*1.557)/FECONM) * PBNL * PSNLCT 
      ELSE 
        PUYSCT = ((ABURN(3)*1.557)/FECONM) * PBNL * PSNLCT 
      ENDIF 
C 
      RETURN 
      END 
