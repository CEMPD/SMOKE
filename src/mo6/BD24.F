      BLOCK DATA BD24 
C 
C  BLOCK DATA Subprogram 24: /NAMES1/,/STRING/,/PRINTC/,/FLAGS5/ 
C 
C  BLOCK DATA 24 initializes character strings used in prompts & output formats 
C  and zeroes pointers into character string arrays. 
C 
C  Changes : 
C 
C  19 Jan 01 @EPA-bag Added HDR2007 for no 2007 rule for heavy duty 
C            vehicle emissions reductions 
C  13 Oct 00 AIR Task 02: Bug Fix 273: Removed SOAKFLAG code to allow for 
C            multi-scenario SOAK DISTRIBUTION commands. 
C  06 Oct 00 @EPA-bag BUG155&137 Remove FTP_FLAG, RN_FLAG 
C  27 Jun 00 AIR Task 04: Added initialization of EVAP_FLAG 
C  14 Jun 00 AIR Task 03: Added WXRVP_FLAG. Removed CO_FLAG 
C  20 Apr 00 @EPA-ddj E-7 Added FUEL P label for fuel program flag to 
C            NAMFLG array. 
C 28 Apr 00 @DynTel-ddj 1-020, Added SFTP_AC_FLAG and SFTP_SP_FLAG. 
c  24 Apr 00 AIR Task 02: Fixed erroneous use of CHAR in DATA statement. 
C  13 Jan 00 @DynTel-ddj 1-013, Added AC_FLAG. 
C  30 Dec 99 @DynTel-ddj 1-021, Added TIER2_FLAG. 
C   7 Sep 99 @DynTel-ddj 2-710 Added FTP_INPUT_FLG. 
C  13 Jul 99 @DynTel-ddj 2-699 Changed implementation logic of FTP_FLAG 
C            back so that FTP_FLAG is equal to 1 when invoked and 0 
C            otherwise. 
C   1 Jul 99 @DynTel-HXQ 2-705 Add HDCGPM_FLAG. 
C  28 May 99 @DynTel-ddj 2-700  Changed implementation logic of FTP_FLAG 
C            to be consistent with other flags. 
C  21 May 99 @DynTel-ddj, Added FTP_FLAG. 
C  12 May 1999 @DynTel-HQ Added SPD_FLAG, RN_FLAG 
C  12 Apr 1999 @DynTel-ddj 2-698 Added flag SPD_FLAG used to bypass CO offset ... 
C  12 Apr 1999 @DynTel-ddj 2-000 Added flags CO_FLAG, RPV_FLAG, 
C              SOAK_FLAG, and TCF_FLAG used to bypass CO offset 
C              correction, RPV correction, Soak Time correction, 
C              Temperature Correction Factor calculation and 
C              speed correction respectively in MAINTCF.FOR. 
C  11 Mar 1999 @DynTel-MLA 2-701 Added common block PRINTC. Deleted USNAME, 
C              IUSD, NUSD, and the obsolete include file USDATA. 
C  20 Nov 1998 @DynTel-MLA 2-663 Removed RDTTCF (Not used anywhere in code) 
C  14 Sep 1998 DynTel@bag  2-000 Removed TI1FLG 
C 
C  Common block array subscripts: 
C 
C  NAMFLG(10)     -  NAMFLG ( IFG ) 
C  NAMMMR(2,3,8)  -  NAMMMR ( ICH, IG, IV ) 
C  NAMRVP(4)      -  NAMRVP ( IRVP ) 
C  NAMTEM(2)      -  NAMTEM ( ITEM ) 
C  NOYES(2)       -  NOYES  ( <ILDT,INO,IMFLAG,ATPFLG> ) 
C 
C  Common block dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------- 
C  /NAMES1/: 
C 
C  NAMFLG   C*8  control flag names for prompts and error messages 
C  NAMMMR   C*8  prompts for replacement AMAR & JULMYR, by JDX & IV 
C  NAMRVP   C*8  name of RVP falling outside allowable range 
C  NAMTEM   C*4  'mini' & 'maxi' prefixes for fuel tank temp error messages 
C 
C  /PRINTC/: 
C 
C  NEWPAGE  C*1  Form feed (page break) printer control character.  
C 
C  /STRING/: 
C 
C  NOYES    C*4  character string for no/yes 
C  COMMA    C*4  character string for comma 
C  PERIOD   C*4  character string for period 
C  COLON    C*1  character string for colon 
C 
C  /FLAGS5/: 
C 
C  TCF_FLAG     temperature correction 
C  RVP_FLAG     RVP correction 
C  SOAK_FLAG    Soak time calculation 
C  SPD_FLAG     Speed Correction 
C  HDCGPM_FLAG  Unit Change 
C  WXRVP_FLAG   Weathered RVP 
C  EVAP_FLAG    Evaporative Emissions Flag 
C  HDR2007  I    1 = DEFAULT 2007 rule for heavy duty vehicle emissions 
C                    reductions is in effect 
C                0 = 2007 rule for heavy duty vehicle emissions 
C                    reductions is disabled; No 2007 rule for reductions of 
C                    heavy duty vehicle emissions reductions 
C 
C 
C  Notes: 
C 
C  BD 24 was modified in Version 01 to correct the definition order of USNAME 
C  and IUSD. 
C  7-May-1994 @ CSC-tmm Request #2-441 Modified code to permit the detergent gasoline 
C     parameters to be read in and stored in COMMON. Added DTGFLG to NAMES1. 
C 
C  May-1995 @ CSC-zp Req. 2-495 "Add MOB5a-H to MOBILE5b" 
C  Added the new  names for the flags from the special record when 
C  IMFLAG=6. 
C  These flags indicate the following:  
C   NUMIMP-the number of the I/M descriptive records;  
C   TI1FLG-Yes or No  for "Tier 1 " program; 
C   ITTCFG Yes or No for "TTC" program, if ITTCFG=3 another record should be 
C          in the  One-Time Data Section               
C   RSDFLG- Yes or No for "RSD" program 
C 
C  28 Aug 97 DynTel@MLA 2-622  Added the value 'SSFLAG' to the NAMFLG 
C            array for Mobile6. 
C  16 December 1997 Dyntel@MMS modified multispeed to single speed. 
C                   SPDFLG renamed to RNLFLG, and reduced by 1. 
C  24 Jun 1998 Dyntel@bag NAMFLG array arranged to represent only flags  
C              to be used in MOBILE6. 
C 
C  24 Jul 98 DynTel@MLA 2-000  Deleted NAMMMR (no longer ref'd in code.) 
C 
      IMPLICIT NONE 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'NAMES1.I' 
      INCLUDE 'PRINTC.I' 
      INCLUDE 'STRING.I' 
C 
C  /FLAGS5/: These names are used in input routines and in QUITER. 
C 
      DATA RVP_FLAG, SOAK_FLAG, TIER2_FLAG, AC_FLAG, 
     &     TCF_FLAG, SPD_FLAG, HDCGPM_FLAG,HDR2007, 
     &     SFTP_AC_FLAG, SFTP_SP_FLAG, WXRVP_FLAG, EVAP_FLAG  / 12*1/ 
      DATA FTP_INPUT_FLG, SPDFLG  / 2*0 / 
C 
C  /NAMES1/: These names are used in input routines and in QUITER. 
C 
      DATA NAMFLG /'  IREJN ',' ITTCFG ',' NUMIMP ',' FUEL P ', 
     &  ' RFGFLG ',' RSDFLG ','       ',' WINFLG ',' IMFLAG ', 
     &  ' RLFLAG '/ 
      DATA NAMRVP /'Period 1','Period 2','Refuelng','Weatherd'/ 
      DATA NAMTEM/'mini','maxi'/ 
C 
C  /PRINTC/: Printer control characters.  
C 
      DATA NEWPAGE /' '/ 
C 
C  /STRING1/: 
C 
      DATA NOYES/'No  ','Yes '/ 
      DATA COMMA/',   '/,PERIOD/'.   '/,COLON/':'/ 
C 
      END 
