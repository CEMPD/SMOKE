      SUBROUTINE BDSAVE 
C 
C  BDSAVE saves the common block variables and arrays that have to have their 
C  default values restored at the outset of each subsequent call to MOBILE. 
C 
C  Called by MOBILE. 
C 
C  Calls none. 
C 
C  Changes: (Last change first) 
C 
C  30 Apr 01 AIR Task 12: Added VA022 for T2PM_FRAC 
C  14-DEC-01 @EPA-bag        REMOVED SAVE02.I, SAVE03.1 and associated variables 
C  21 Nov 00 @EPA-bag bug300 rename soak dist. save array; now named va021; 
C            removed save10.i 
C  26 Jul 00 AIR Task 07: Fixed Bug 201: Failure to reinitialize the 
C            Soak Distribution array 
C  23 May 00 @EPA-djb Add VA019 & VA020 for STARTSPERDAY and STARTDIST. 
C  07 Apr 00 @EPA-djb Changed VMTMIX restore to VMTUSR. 
C  24 Jan 00 @DynTel-MLA 1-000  Corrected typos in comments. 
C  14 Jan 00 @Dyntel - jwr 1-017 Added VA016A for T2_EVAP 
C  30 Dec 99 @DynTel-ddj 1-021, Added VA016 and VA017 for T2_BIN_FRAC 
C            and TIER2_50K respectively. 
C  22 Jun 99 @DynTel-MLA 2-694  Removed P94REF and P94SEF. 
C  25 Mar 99 @DynTel-HQ 2-698 Speed related coefficient 
C   7 Oct 98 @DynTel-HQ 2-678 Comment statements 
C  29 Sep 98 @DynTel-HQ 2-678 Conversion Factor 
C  30 Sep 96 DynTel@MLA 2-226  Added IVTYPE.I for Mobile6. Changed the I3 
C            loop index for TAMZML and TAMDR from 4 to 14 because the arrays 
C            have been expanded to accomodate new vehicle types in Mobile6. 
C 
C  20 Apr 98 DynTel@ID 2-659 Added the part to save data of the hardcoded 
C  LEVIMP and LEVSTD arrays. 
C 
C  Input on call: 
C 
C    Constants:  
C    MAXVEH from include file IVTYPE.I 
C 
C     common blocks: 
C     /ETPPAR/ ETP,ETPPIV,NETPYR 
C     /IVTYPE/ MAXVEH 
C     /MAXIMA/ MAXPOL,MAXREG,MAXYRS 
C     /MYRSAV/ AMAR,JULMYR 
C     /TAMEQ1/ TAMZML,TAMDR,TAMA50 
C     /VMXCOM/ VMTMIX 
C     /BASEQ9/ LEVIMP 
C     /HDCCOM/ HDCFAC 
C     /SPEED9/ TVMT, HVMT, FVMT, SVMT, LOCAL, RAMP, 
C                OCECO, OCY, OCX, LEVEL, EMSFC 
C     /SKDIST/ SOAKDIST 
C 
C  Output on return: 
C 
C     common blocks: 
C     /SAVE04/ VA007,VA008 
C     /SAVE05/ VA009,VA010 
C     /SAVE06/ VA011,VA012,IVA013,VA014 
C     /SAVE07/ VA015,VA016,VA016A,VA022 
C     /SAVE08/ VA018, VA019, VA020, VA021 
C     /SAVE09/ STVMT, SHVMT, SFVMT, SSVMT, SLOCAL, SRAMP, 
C              SOCECO, SOCY, SOCX, SLEVEL, SEMSFC, 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  I<d>     I    array subscripts for the save areas 
C 
C  Notes: 
C 
C  BDSAVE was added for Version 01's conversion of MOBILE4 to a subroutine. 
C 
C  In Version 03, "IF(I2.LE.2)" removed before TAMZML line (dim now 3). 
C  TAMA50 is left alone, since for the release version, users will not be 
C  given the option to enter their own TAMA50. 
C 
      IMPLICIT NONE 
      INCLUDE 'ETPPAR.I' 
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
      INCLUDE 'SKDIST.I' 
      INCLUDE 'TIER2PM.I' 
C 
      INCLUDE 'SAVE01.I' 
      INCLUDE 'SAVE04.I' 
      INCLUDE 'SAVE05.I' 
      INCLUDE 'SAVE06.I' 
      INCLUDE 'SAVE07.I' 
      INCLUDE 'SAVE08.I' 
      INCLUDE 'SAVE09.I' 
C 
C  Subroutine CHKMRN 
C 
      VA007 = AMAR 
      VA008 = JULMYR 
C 
C  Subroutine GETTAM 
C 
      VA009 = TAMZML 
      VA010 = TAMDR 
C 
C  Subroutine GETVMT 
C 
      VA011 = VMTUSR 
C 
C  Subroutine GETETP 
C 
      VA012 = ETP 
C 
      IVA013 = IETPYR 
      VA014  = ETPPIV 
C 
C  Subroutine GET94P - '94P' IMP data 
C 
      VA015 = LEVIMP 
C 
C  Subroutine GET_TIER2 - TIER2 FRACTIONS 
C 
      VA016 = T2_BIN_FRAC 
C 
      VA016A = T2_EVAP 
C 
C  Subroutine GET_T2_CERT - TIER2 50K CERTIFICATION STANDARDS 
C 
      VA017 = TIER2_50K 
C 
C  Subroutine GETCONV 
C 
      VA018 = HDCFAC 
C 
C  SUBROUTINE GETFVMT.FOR 
C 
      STVMT = TVMT 
C 
C  SUBROUTINE GETHVMT.FOR 
C 
      SHVMT = HVMT 
      SFVMT = FVMT 
C 
C  SUBROUTINE GETSVMT.FOR 
C 
      SSVMT = SVMT 
C 
C  SUBROUTINE GETSPCO.FOR 
C 
      SLOCAL = LOCAL 
      SRAMP  = RAMP 
      SOCY   = OCY            
      SOCX   = OCX 
      SLEVEL = LEVEL 
C 
      SOCECO = OCECO 
      SEMSCF = EMSCF 
C 
C  SUBROUTINE GETSPD 
C 
      VA019 = STARTSPERDAY 
C 
C  SUBROUTINE GETSDIST 
C 
      VA020 = STARTDIST 
C 
C  SUBROUTINE GETSOAK 
C 
      VA021 = SOAKDIST 
C 
C  Subroutine GETT2PMFRC - TIER2 PM FRACTIONS 
C 
      VA022 = T2PM_FRAC 
C 
      RETURN 
      END 
