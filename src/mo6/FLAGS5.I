C 
C  Changes: (Last change first) 
C 
C  19 Jan 01 @EPA-bag Added HDR2007 for no 2007 rule for heavy duty 
C            vehicle emissions reductions 
C  13 Oct 00 AIR Task 02: Bug Fix 273: Removed SOAKFLAG to allow for 
C            multi-scenario SOAK DISTRIBUTION commands. 
C  06 Oct 00 @EPA-bag BUG155&173 Remove FTP_FLAG and RN_FLAG 
C  26 Jun 00 AIR Task 04: Added EVAP_FLAG 
C  14 Jun 00 AIR Task 03: Added WXRVP_FLAG. Removed CO_FLAG. 
C  13 Apr 00 @EPA-ddj 1-013 Added SFTP_SP_FLAG and SFTP_AC_FLAG. 
C  13 Jan 00 @DynTel-ddj 1-013 Added AC_FLAG. 
C  30 Dec 99 @DynTel-ddj 1-021 Added TIER2_FLAG. 
C   7 Sep 99 @DynTel-ddj 2-710 Added FTP_INPUT_FLG. 
C   1 Jul 99 @DynTel-HXQ 2-705 Add HDCGPM_FLAG. 
C  28 May 99 @DynTel-ddj 2-700 Changed RPV_FLAG to RVP_FLAG. 
C  21 May 99 @DynTel-ddj, Added FTP_FLAG. 
C  27 Apr 1999 @DynTel-HQ Added SPD_FLG 
C   2 Apr 1999 @DynTel-DDJ NEW File 
C 
C The data is initialized in BD24. 
C 
      INTEGER AC_FLAG 
      INTEGER RVP_FLAG 
      INTEGER SOAK_FLAG 
      INTEGER TCF_FLAG 
      INTEGER SPD_FLAG 
      INTEGER HDCGPM_FLAG 
      INTEGER FTP_INPUT_FLG 
      INTEGER SPDFLG 
      INTEGER TIER2_FLAG 
      INTEGER SFTP_SP_FLAG 
      INTEGER SFTP_AC_FLAG 
      INTEGER WXRVP_FLAG 
      INTEGER EVAP_FLAG 
      INTEGER HDR2007 
C 
      COMMON /FLAGS5/RVP_FLAG, SOAK_FLAG, TIER2_FLAG, WXRVP_FLAG, 
     &            TCF_FLAG, SPD_FLAG, HDCGPM_FLAG, FTP_INPUT_FLG, 
     &            SPDFLG, AC_FLAG, SFTP_AC_FLAG, 
     &            SFTP_SP_FLAG,EVAP_FLAG,HDR2007 
