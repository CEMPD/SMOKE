      SUBROUTINE RESTOR 
C 
C  RESTOR restores the BLOCK DATA "default" values to the common block 
C  variables and arrays that can be modified through user run level data 
C  option selections (previously called the one-time data section). In 
C  subroutine MOBILE4 these must be restored, since the user may elect 
C  different options on successive calls.  Hence, RESTOR is called at 
C  the outset of each call to MOBILE after the first. The subroutine 
C  BDSAVE should be used to store the default values in variables 
C  typically named VAxxx. (Restoration of the default values for 
C  variables associated with scenario level commands should be 
C  done in the subroutine INISCN.) 
C 
C 
C  Called by MOBILE. 
C 
C  Calls none. 
C 
C  Changes: (Last change first) 
C 
C  30 Apr 01 AIR Task 12: Added VA022 for T2PM_FRAC 
C  21 Nov 00 @EPA-bag bug300 remove soak dist.; proper scenario level 
C            initializations are done in INISCN 
C  26 Jul 00 AIR Task 07: Fixed Bug 201: Failure to reinitialize the 
C            Soak Distribution array 
C  23 May 00 @EPA-djb Add VA019 & VA020 for STARTSPERDAY and STARTDIST. 
C  17 Apr 00 @EPA-djb  Changed VMTMIX save to VMTUSR. 
C  14 Jan 00 @Dyntel - jwr 1-017 Added VA016A for T2_EVAP 
C  24 Jan 00 @DynTel-MLA 1-000  Corrected typos in comments. 
C  30 Dec 99 @DynTel-ddj 1-021, Added VA016 and VA017 for T2_BIN_FRAC 
C            and TIER2_50K respectively. 
C  01 Oct 99 @Dyntel JWR 2-705 Removed subroutine GETIMC which had been 
C            incorrectly redefining cutpoints.   
C             
C  22 Jun 99 @DynTel-MLA 2-694  Removed P94REF and P94SEF. 
C  25 Mar 99 @DynTel-HQ 2-698 Speed related coefficient 
C   3 Mar 99 @DynTel-MLA 2-659  Moved declaration of VA015 from the output 
C            section of the header to the input section. Added declaration 
C            for MAXLYR. 
C  11 Feb 99 @DynTel-MLA 2-663  Added code to initialize cutpoints. 
C  31 Dec 98 @DynTel-ZK 2-000 Explicit Typing 
C   7 Oct 98  @DynTel-HQ 2-678 Comment statements 
C  28 Sep 98  @DynTel-HQ 2-678 Conversion Factor 
C 
C  Input on call: 
C 
C    Constants:  
C    MAXVEH from include file IVTYPE.I 
C 
C     common blocks: 
C     /BASEQ9/ MAXLYR,VA015 
C     /MAXIMA/ MAXPOL,MAXREG,MAXYRS 
C     /SAVE04/ VA007,VA008 
C     /SAVE05/ VA009,VA010 
C     /SAVE06/ VA011,VA012,IVA013,VA014 
C     /SAVE07/ VA015, VA016, VA017,VA016A,VA022 
C     /SAVE08/ VA018 
C     /SAVE09/ STVMT, SHVMT, SFVMT, SSVMT, SLOCAL, SRAMP, 
C              SOCECO, SOCY, SOCX, SLEVEL, SEMSFC 
C     /SAVE10/ SSOAKDIST 
C 
C  Output on return: 
C 
C     common blocks: 
C     /BASEQ9/ LEVIMP 
C     /ETPPAR/ ETP,ETPPIV,NETPYR 
C     /HDCCOM/ HDCFAC 
C     /IMPAR4/ CUTHC,CUTCO,CUTNO 
C     /MYRSAV/ AMAR,JULMYR 
C     /TAMEQ1/ TAMZML,TAMDR,TAMA50 
C     /VMXCOM/ VMTMIX 
C     /SPEED9/ TVMT, HVMT, FVMT, SVMT, LOCAL, RAMP, 
C                OCECO, OCY, OCX, LEVEL, EMSFC 
C     /T2_BIN_FRAC/ 
C     /TIER2-50K/ 
C     /SKDIST/ SOAKDIST 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  I<d>     I    array subscripts for the save areas 
C 
C  Notes: 
C 
C  RESTOR was added for Version 01's conversion of MOBILE4 to a subroutine. 
C 
C  In Version 03, "IF(I2.LE.2)" removed before TAMZML line (dim now 3). 
C  TAMA50 is left alone, since for the release version, users will not be 
C  given the option to enter their own TAMA50. 
C 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6. 
C 
C  20 Apr 98 DynTel@ID 2-659 Added restore procedure for LEVIMP and LEVSTD 
C  arrays. 
C 
      IMPLICIT NONE 
      INCLUDE 'SAVE01.I' 
      INCLUDE 'SAVE04.I' 
      INCLUDE 'SAVE05.I' 
      INCLUDE 'SAVE06.I' 
      INCLUDE 'SAVE07.I' 
      INCLUDE 'SAVE08.I' 
      INCLUDE 'SAVE09.I' 
C 
      INCLUDE 'ETPPAR.I' 
      INCLUDE 'IMPAR4.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'MYRSAV.I' 
      INCLUDE 'TAMEQ1.I' 
      INCLUDE 'VMXCOM.I' 
      INCLUDE 'BASEQ9.I' 
      INCLUDE 'HDCCOM.I' 
      INCLUDE 'SPEED9.I' 

      INCLUDE 'STDIST.I' 
      INCLUDE 'TIER2.I' 
      INCLUDE 'TIER2PM.I' 
      INCLUDE 'SKDIST.I' 
C 
C  Subroutine CHKMRN 
C 
      AMAR   = VA007 
      JULMYR = VA008 
C 
C  Subroutine GETTAM 
C 
      TAMZML = VA009 
      TAMDR  = VA010 
C 
C  Subroutine GETVMT 
C 
      VMTUSR = VA011 
C 
C  Subroutine GETETP 
C 
      ETP = VA012 
C 
      IETPYR = IVA013 
      ETPPIV = VA014 
C 
C 
C  Subroutine GET94P - 94P IMP data 
C 
      LEVIMP = VA015 
 
C 
C  Subroutine T2_BIN_FRAC - Tier2 Fractions 
C 
      T2_BIN_FRAC = VA016 
C 
C 
      T2_EVAP = VA016A 
C 
C Subroutine TIER2_50K - Tier2 50K Certification Standards 
C 
      TIER2_50K = VA017 
C 
C  Subroutine GETCONV 
C 
      HDCFAC = VA018 
C 
C  SUBROUTINE GETFVMT.FOR 
C 
      TVMT = STVMT 
C 
C  SUBROUTINE GETHVMT.FOR 
C 
      HVMT = SHVMT 
      FVMT = SFVMT 
C 
C  SUBROUTINE GETSVMT.FOR 
C 
      SVMT = SSVMT 
C 
C  SUBROUTINE GETSPCO.FOR 
C 
      LOCAL = SLOCAL 
      RAMP  = SRAMP 
      OCY   = SOCY 
      OCX   = SOCX 
      LEVEL = SLEVEL 
C 
      OCECO = SOCECO 
      EMSCF = SEMSCF 
C 
C  SUBROUTINE GETSPD 
C 
      STARTSPERDAY = VA019 
C 
C  SUBROUTINE GETSDIST 
C 
      STARTDIST = VA020 
      SDIST = 0 
C 
C  Subroutine T2PM_FRAC - Tier2 PM Fractions 
C 
      T2PM_FRAC = VA022 
C 
      RETURN 
      END 
