      SUBROUTINE INIFLG() 
C 
C  Initialize flags. Most of these were set in Mobile5 when the Control 
C  Section input records were read and processed. 
C 
C  Called by MOBILE. 
C 
C  Changes: (Last change first) 
C 
C  23 May 01 AIR Task 17: Added initialization of PM speed cf flag. 
C  09 May 01 @EPA-elg Deactivated POST96IM Flag by setting it to 1. 
C  16 Apr 01 @EPA-elg Added POST96IM Flag for I/M diagnostic label. 
C  19 Jan 01 @EPA-bag Added HDR2007 for no 2007 rule for heavy duty 
C            vehicle emissions reductions 
C  20 Nov 00 @EPA-elg Added NO CLEAN AIR ACT initialization of CAAFLG 
C 27 Nov 00 @EPA-bag bug#301 Removed FVMTFLG, HVMTFLG, SVMTFLG; initialize 
C           GLVMTMIX to VMTUSR, GLHVMT to HVMT, GLFVMT to FVMT, and GLSVMT 
C           to SVMT 
C 20 Nov 00 @EPA-bag bug#303 Initialized all values of the TLVMTU array to -1 
C           to fix weekly index check on user supplied trip length fraction 
C           values, i.e., the TLVMTU's 
C  13 Oct 00 AIR Task 02: Bug Fix 273: Removed SOAKFLAG to allow for 
C            multi-scenario SOAK DISTRIBUTION commands. 
C  06 Oct 00 @EPA-bag BUG155&137 Remove FTP_FLAG and RN_FLAG 
C  06 Oct 00 @EPA-BG  Bug 155, Removed IDLFLG, Removed INDUSTRY AVERAGE FUEL 
C           parameter GLIAFCNT 
C  21 Sep 00 AIR Task 07: Bug Fix 256 Added initialization of GSUNRISE, 
C            GSUNSET, GPEAKSUN, GCLD_CVR 
C  15 Sep 00 AIR Task 08: Included evap I/M initialization for EVIMFLAG. 
C  08 Aug 00 @EPA-elg Added NO_EVAPCALC to shut down evap calculations. 
C  09 Aug 00 @EPA-eg  bug#216; added initialization of GLFUEL_P 
C  25 Jul 00 @EPA-bag bug#54; removed TEMFLG 
C  21 Jul 00 @EPA-RAG running loss; removed pass/fail fractions, PFRATE's, 
C            for purge/pressure evap. tests; added NO_RSTLOS 
C 14 Jun 00 AIR Task 03: Removed CO_FLAG. 
C 08 Jun 00 AIR Task 03: Removed nonblock DO constructs. 
C 19 May 00 AIR Task 03: Added WXRVP_FLAG (Weathered RVP). 
C 13 Jun 00 @EPA-ddj Bug140, Changed GLDGPEY from 1998 to 1996 and 
C           changed GLFIDFC from 0.019 to 0.04. 
C 16 May 00 @EPA-ddj BUG87, Removed STARFG. 
C  2 May 00 @EPA-ddj D-28,  Added initialization of FUEL_P. 
C 28 Apr 00 @DynTel-ddj 1-020, Added SFTP_AC_FLAG and SFTP_SP_FLAG. 
C 16 Apr 00 AIR Task 02: Made VOC the default speciation for HC 
C 10 Apr 00 AIR Task 02: Removed PRTFLG 
C 04 Apr 00 AIR Task 02: Removed ASTM class as input. 
C 20 Feb 00 @DynTel-RAG 1-005 Added evaporative liquid leak flag. 
C 13 Jan 00 @DynTel-ddj 1-013, Removed code for wetbulb and drybulb 
C           temperatures and trailer and towing parameters.  Added 
C           AC_FLAG. 
C 30 Dec 99 @DynTel-ddj 1-021 Added Tier2 Flag. 
C 10 Nov 99 @DynTel-BG 1-006  Added grams/mile index, GPM, for start emissions. 
C           HCFLAG initialized to 1 for standard descriptive output (no 
C           expanded evap. emissions). 
C 05 Oct 99 @Dyntel-jwr 2-705 ITTCFG initialized to 2 so 
C           trained techs are assumed as the default setting.  
C 15 Sep 99 @DynTel-MLA 2-000  Removed GLSCNAME. 
C  7 Sep 99 @DynTel-ddj 2-710 Added FTP flags. 
C 16 Aug 99 @DynTel-bg  2-663 Renamed SIMFLG and SIMONLY to IMFLG_NP and 
C           IMONLY_NP, respectively. These are for the 'non-periodic' COIM 
C           and OBD programs. 
C   1 Jul 99 @DynTel-HXQ 2-705 Add HDCGPM_FLAG. 
C  28 May 99 @DynTel-ddj 2-700 Changed RPV_FLAG to RVP_FLAG. 
C  12 May 99 @DynTel-HQ initialize the flags to turn of the Speed correction  
C                        and running loss 
C  15 Apr 99 @DynTel-HQ initialize the flags for user input speed correction 
C                       related coefficient 
C  2  MAR 99 @DynTel-bag 2-664 Expanded GLVMTU and TLVMTU for hourly running losses; 
C            Deleted RNLFLG. 
C  16 Feb 99 @Dyntel-RJD 2-663  Reset default values for optional 
C            81-93 model year LDGV/T selective I/M program flags. 
C  11 Feb 99 @DynTel-MLA 2-663  Added initialization of STARFG. 
C  10 Dec 98 @DynTel-MLA 2-663  Added the IMONLY_NP flag. 
C  20 Nov 98 @DynTel-MLA 2-663  ITTCFG and TTCEFF were changed from integers 
C            to arrays (one value for each possible I/M program). 
C  20 Nov 98 @DynTel-ZK 2-000 Explicit Typing 
C  30 Oct 98 @DynTel-BAG 2-665  Added hourly temperature array, 
C            GLTHRLY(IHR) and GLDEL_TEM(IHR) for global declarations. 
C  28 Sep 98  @DynTel-HQ 2-678 Conversion Factor 
C 
C  Input on call: 
C 
C    common blocks: 
C    /ATPAR4/ EFFSAV 
C    /CHKLST/ MAXREQ,MAXIMPAR,MAXIMPGM 
C    /UNIT_PARAM/ GPM 
C 
C  Output on return: 
C    modules: 
C    LEAKYDAT - ILIQUID_LEAKS 
C    common blocks: 
C    /ALUHIN/ GLAC, GLAHUM, GSUNRISE, GSUNSET, GPEAKSUN, GCLD_CVR 
C    /ATPAR4/ EFF 
C    /CHKLST/ GLCHECK,IMCHECK 
C    /CITPAR/ SCNAME 
C    /CITRV1/ GLRVPBAS 
C    /DTGGAS/ GLDGPEY,GLDGPHSF,GLDGPSY,GLFIDEF,GLFIDFC,GLIVDEF,GLIVDFC 
C    /FLAGS1/ GLOXFLAG,TAMFLG,VMFLAG,OXYFLG,DSFLAG,DTGFLG, 
C             GLDSFLAG,GLDTGFLG,GLVMFLG,USR94P,CONVFLG 
C    /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,GLALHFLG 
C    /FLAGS3/ ATPFLG,RLFLAG 
C    /FLAGS4/ NMHFLG,HCFLAG,COLDFG 
C    /FLAGS5/ WXRVP_FLAG 
C    /IMPAR5/ ATPUEF, PPGUEF 
C    /IMPAR6/ IFILET,IMEFFR 
C    /IMPAR7/ SVALTC 
C    /LDGSIM/ IMFLG_NP,IMONLY_NP 
C    /OPCNTL/ OPTABL 
C    /OXY1/   GLIGAS,GLOXCNT,GLSHARE 
C    /P94INF/ P94AR 
C    /PART5/  PMSPDFLG 
C    /PRGCH2/ PRGFLG 
C    /PRSCH2/ PRSFLG 
C    /RESTLS1/ NO_RSTLOS 
C    /RFORM1/ GLRFG,GLWIN,RFGFLG,WINFLG 
C    /RSDINF/ RSDCLEAN,RSDDIRTY 
C    /RUNLS1/ GLVMTU,TLVMTU 
C    /TEMPS/  GLTHRLY, GLDEL_TEM 
C    /TTCINF/ ITTCFG,TTCEFF 
C 
C    Mobile5 Control Section flags are set to default values. 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  I        I    Loop control variable 
C  IH       I    Hourly loop control variable 
C 
C  Notes: 
C 
C  INIFLG was added for Mobile6. 
C 
C  1 Jul 98 @DynTel - MLA 2-000 Removed PROJID and its include file PROJEC. 
C           (PROJID is obsolete in Mobile6.) 
C 20 Jul 98 @DynTel  NH  2-000 
C 
C ************************************************************************** 
C ************************************************************************** 
C  This subroutine initializes the flags that were set as we read the 
C  Control Section data in earlier versions of the Mobile model. 
C ************************************************************************** 
C ************************************************************************** 
C 
      USE LEAKYDAT, ONLY : ILIQUID_LEAKS 
      IMPLICIT NONE 
C 
      INCLUDE 'ALUHIN.I' 
      INCLUDE 'ATPAR4.I' 
      INCLUDE 'CHKLST.I' 
      INCLUDE 'CITPAR.I' 
      INCLUDE 'CITRV1.I' 
      INCLUDE 'DTGGAS.I' 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'IMPAR4.I' 
      INCLUDE 'IMPAR5.I' 
      INCLUDE 'IMPAR6.I' 
      INCLUDE 'IMPAR7.I' 
      INCLUDE 'LDGSIM.I' 
      INCLUDE 'OPCNTL.I' 
      INCLUDE 'OXY1.I' 
      INCLUDE 'P94INF.I' 
      INCLUDE 'PART5.I' 
      INCLUDE 'PRGCH2.I' 
      INCLUDE 'PRSCH2.I' 
      INCLUDE 'RESTLS1.I' 
      INCLUDE 'RFORM1.I' 
      INCLUDE 'RSDINF.I' 
      INCLUDE 'RUNLS1.I' 
      INCLUDE 'SPEED9.I' 
      INCLUDE 'SULFUR.I' 
      INCLUDE 'TEMPS.I' 
      INCLUDE 'TTCINF.I' 
      INCLUDE 'UNIT_PARAM.I' 
      INCLUDE 'VMXCOM.I' 
C 
      INTEGER  IATP 
      INTEGER  ICOL 
      INTEGER  IROW 
      INTEGER, PARAMETER :: IYES = 1 
C 
C     Initialize a checklist of required inputs. This is implemented only 
C     for scenarios at the present time. No global (or One-Time) data is 
C     required to run this version of the Mobile6 code. The scenario inputs 
C     that are required in this version of the code to run Mobile6 are the 
C     the min/max temperatures, fuel RVP, region, calendar year, 
C     average speed, and ambient temperature. (This list will certainly 
C     change.) The first two items (min and max temps, and RVP), even  
C     though they are considered as scenario inputs, may also be 
C     entered globally. That's why the checklist is initialized here. 
C 
      GLCHECK = 0 
C 
C  /ALUHIN/ 
C  Initialize the ALH parameters. 
C 
      GLAC        = 1.0 
      GLAHUM      = 75.0 
      GSUNRISE    = 1 
      GSUNSET     = 16 
      GPEAKSUN(1) = 5 
      GPEAKSUN(2) = 11 
      GCLD_CVR    = 0.0 
C 
C     /ATPAR4/ 
C 
C     Reset EFF to default values. These may be changed in GETEFR. 
C 
      DO IROW = 1,15 
        DO IATP = 1,2 
          DO ICOL = 1,3 
            EFF(2,ICOL,IATP,IROW) = EFFSAV(2,ICOL,IATP,IROW) 
          END DO 
        END DO 
      END DO 
C 
C     /CHKLST/ 
C 
      IMCHECK = 1 
C 
C     /FLAGS1/ 
C 
      USR94P  = 1 
      TAMFLG  = 1 
      VMFLAG  = 1 
      SPCOFLG = 1 
      SCFFLG  = 1 
      GLVMFLG = 1 
      OXYFLG  = 1 
      GLOXFLAG= 1 
      DSFLAG  = 1 
      GLDSFLAG= 1 
      DTGFLG  = 1 
      GLDTGFLG= 1 
      CONVFLG = 1 
C 
C     /FLAGS2/ 
C 
      MYMRFG  = 1 
      NEWFLG  = 1 
      IMFLAG  = 1 
      GLALHFLG= 1 
      ALHFLG  = 1 
      EVIMFLAG= 1 
C 
C     /FLAGS3/ 
C 
      ATPFLG = 1 
      RLFLAG = 1 
      CAAFLG = 0 
      HDR2007= 1 
      POST96IM = 1 
C 
C     /FLAGS4/ 
C 
      NMHFLG = 3 
      HCFLAG = 1 
      COLDFG = 2 
C 
C     /FLAGS5/ 
C 
      AC_FLAG       = 1 
      RVP_FLAG      = 1 
      SOAK_FLAG     = 1 
      TCF_FLAG      = 1 
      SPD_FLAG      = 1 
      HDCGPM_FLAG   = 1 
      TIER2_FLAG    = 1 
      SFTP_AC_FLAG  = 1 
      SFTP_SP_FLAG  = 1 
      FTP_INPUT_FLG = 0 
      SPDFLG        = 0 
      WXRVP_FLAG    = 1 
C 
C     /IMPAR4/ 
C 
      ATPUEF = 0.0 
      PPGUEF = 0.50 
C 
C     /IMPAR6/ 
C 
      IFILET    = 1 
      IMEFFR    = 1 
C 
C     /IMPAR7/ 
C 
      SVALTC = 0 
C 
C    /LDGSIM,RSDINF/ 
C 
      IMFLG_NP = 0 
      RSDCLEAN = 0 
      RSDDIRTY = 0 
      IMONLY_NP= 1 
C 
C     /P94INF/ 
C 
      P94AR(1) = ' ' 
      P94AR(2) = ' ' 
      P94AR(3) = ' ' 
C 
C     /OPCNTL/ 
C 
      OPTABL = 1 
C 
C     /P94INF/ 
C 
      P94AR = ' ' 
C 
C     /PRGCH2/ 
C 
      PRGFLG=1 
C 
C     /PRSCH2/ 
C 
      PRSFLG=1 
C 
C    /TTCINF/ 
C 
C     Initialize ITTCFG to 2 to make fully trained 
C     technicians the default setting in IM240 
      ITTCFG = 2 
      TTCEFF = 1.0 
C 
C     The following values are the global defaults for values that might 
C     be selectively replaced for individual scenarios. 
C 
C     /CITPAR/ 
C 
      SCNAME = ' ' 
C 
C     /OXY1/     GLSHARE is the global data for SHRMKT. 
C                GLOXCNT is the global data for OXYCNT. 
C                GLIGAS is the global data for IGASHW. 
C 
      GLSHARE(1)=1.0 
      GLSHARE(2)=0.0 
      GLSHARE(3)=0.0 
C 
      GLOXCNT(1)=0.0 
      GLOXCNT(2)=0.0 
C 
      GLIGAS=1 
C 
C     /RFORM1/  GLRFG is the global version of RFGFLG 
C               GLWIN is the global version of WINFLG 
C 
      GLRFG   = 1 
      RFGFLG  = 1 
      GLWIN   = 0 
      WINFLG  = 1 
      FUEL_P  = 1 
C 
C     /RESTLS1/ 
C 
      NO_RSTLOS=0 
      NO_EVAPCALC=0 
C 
C     /RUNLS1/ 
C 
      GLVMTU = -1 
      TLVMTU = -1 
C 
C     /DTGGAS/ 
C 
      GLDGPSY    = 1995 
      GLDGPEY    = 1996 
      GLDGPHSF   = 0.50 
      GLFIDEF(1) = 0.025 
      GLFIDEF(2) = 0.0165 
      GLFIDEF(3) = 0.0 
      GLFIDFC    = 0.04 
      GLIVDEF(1) = 0.054 
      GLIVDEF(2) = 0.065 
      GLIVDEF(3) = 0.118 
      GLIVDFC    = 0.18 
C 
      GLRVPBAS = 0.0 
      GLTMIN   = 0.0 
      GLTMAX   = 0.0 
      GLTHRLY  = 0.0 
      GLDEL_TEM= 0.0 
C 
C     /SULFUR/ 
C 
      GLFUEL_P  = 1 
C 
C     module LEAKYDAT.FOR 
C     Initialize the liquid leak switch to include liquid 
C     leak emissions to be included in the evaporative sums 
C 
      ILIQUID_LEAKS=IYES 
C 
C     /SPEED9/ 
C 
      GLHVMT=HVMT 
      GLFVMT=FVMT 
      GLSVMT=SVMT 
C 
C    /VMXCOM/ 
C 
      GLVMTMIX=VMTUSR 
C 
C    /PART5/ 
C 
      PMSPDFLG = 0 
C 
 999  RETURN 
      END 
