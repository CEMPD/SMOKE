      SUBROUTINE INISCN(INPREC,NFILE,NSCEN,NRUN) 
C 
C  Print a header for the scenario. 
C  Initialize the flags and variables associated with scenario level 
C  and scenario+run level commands which need to be reset for each 
C  scenario. 
C 
C  Called by DSPTCH. 
C 
C  Calls LJCOMP. 
C 
C  Changes: (Last change first) 
C 
C  30 Aug 01 AIR Task 15-23: Added initialization of PM and Air Toxic checks. 
C  19 Jul 01 @EPA-djb Wish235, Initialize AVGFLG to zero and remove dead code. 
C  21 Feb 01 AIR Task 11: Included code to display scenario in PM runs. 
C  21 Nov 00 @EPA-bag bug300 add soak dist. scenario level initialization 
C  06 Oct 00 @EPA-bag BUG155&137 Remove FTP_FLAG and RN_FLAG 
C  06 Oct 00 @EPA-BG  Bug 155, Removed INDUSTRY AVERAGE FUEL 
C           parameters GLIAFCNT and IAFCNT 
C  21 Sep 00 AIR Task 07: Bug Fix 256, Added initialization of SUNRISE to 
C            GSUNRISE, SUNSET to GSUNSET, PEAKSUN to GPEAKSUN, and 
C            CLD_CVR to GCLD_CVR. 
C  12 Sep 00 @EPA-djb Bug 237, Initialized VMTUSR to GLVMTMIX.  Eliminated 
C            initialization of VMTMIX. 
C  17 Jul 00 AIR Task 06: Added SAWFLAG, HSACTFLAG, DSACTFLAG, FPPRFLAG 
C  21 Jul 00 @EPA-RAG running loss; removed pass/fail fractions, PFRATE's, 
C            for purge/pressure evap. tests 
C  10 Jul 00 @EPA-ddj Bug184 Changed initialization of FUEL_P from 1 to 
C            GLFUEL_P. 
C  14 Jun 00 AIR Task 03: Removed CO_FLAG. 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C   4 May 00 @EPA-ddj E-7, Added S_PPM to initialization. 
C  25 Feb 00 @DynTel-ddj 1-028  Changed index on RFGSYR from 3 to 2. 
C            Also removed obsolete ASTM code. 
C  13 Apr 00 AIR Task 02: Set default ALTITUDE (region) to Low (1) 
C  04 Apr 00 AIR Task 02: Removed ASTM class. Limited the Scenario Title 
C            output to 122 characters. 
C  17 Apr 00 @EPA-djb Add IVC loop and expand USRGSF (IVTEMP) loop to 14. 
C  20 Feb 00 @DynTel-bag 1-005  Added liquid leak switch parameter, ILIQUID_LEAKS 
C  13 Jan 00 @DynTel-ddj 1-013, Removed ALHFLG, wetbulb and drybulb 
C            temperature, and towing and trailer parameters. 
C  29 Nov 99 @DynTel-bag 1-006  Removed HOURLY.I and its associated variables for 
C            hourly evaporative emissions output 
C  15 Sep 99 @DynTel-MLA 2-000  Removed GLSCNAME. 
C  24 Aug 99 @DynTel-MLA 2-684 Initialize database output variables. 
C  13 Jul 99 @DynTel-ddj 2-699 Changed implementation logic of FTP_FLAG 
C            back so that FTP_FLAG is equal to 1 when invoked and 0 
C            otherwise. 
C  22 Jun 99 @DynTel-MLA 2-694 Removed LEV variables APPLEV, LEVYRS. 
C   1 Jun 99 @DynTel-ddj 2-700 Changed default value for FTP_FLAG from 
C            0 to 1 and RPV_FLAG changed to RVP_FLAG. 
C  26 May 99 @DynTel-HXQ 2-698 Delete code for operating mode and 
C                        average speed 
C  27 May 99 @DynTel-BG 2-000 Set default TEST FLAGS options to do the 
C            CO offset, RVP, soak time adjustment to start emissions,  
C            speed, running loss,... corrections. (The parameters for 
C            these switches are all initialized to 1 at the beginning of 
C            run.) 
C  19 May 99 @DynTel-BG 2-701 NEWPAGE[=CHAR(12),form feed] not written to 
C            output for first scenario. Default values of correction 
C            calculation flags are all set to 1, i.e., do the calculation 
C  15 Apr 99 @DynTel-HQ 2-698 operating mode data are initialized 
C  08 Apr 99 @DynTel-bag 2-000  Removed ambient temperature, AMBT 
C  11 Mar 99 @DynTel-MLA 2-701  Added INPREC to the parameter list. This 
C            is the SCENARIO RECORD labeled input record, which now holds 
C            the scenario title. Added optional page break before the header. 
C  2  Mar 99 @DynTel-BAG 2-664  Expanded GLVMTU and TLVMTU arrays for hourly 
C            running losses (weekly and hourly indices). Added weekday/ weekend 
C            index, IWEEK, to scenario intialization. Added initialization of 
C            TEMP_DIST to .FALSE. For temporary output of hourly emissions added 
C            HOURLY.I with the logical variable HOURLY_OUT. 
C  20 Nov 98 @DynTel-ZK 2-000 Explicit Typing 
C  30 Oct 98 @DynTel-BAG 2-665  Added hourly temperature initializations 
C           for the hourly temperature arrays TEMHRLY(24), GLTHRLY(24), DEL_TEM(24), 
C           and, GLDEL_TEM(24). 
C 
C  Input on call: 
C    /PRINTC/ NEWPAGE 
C    /ALUHIN/ GSUNRISE,GSUNSET,GPEAKSUN,GLAC,GLAHUM,GCLD_CVR 
C    /CITRV1/ GBLRVPBAS 
C    /DTGGAS/ GLDGPHSF 
C    /FLAGS1/ GLOXFLAG,GLVMFLG 
C    /FLAGS2/ GLALHFLG 
C    /OXY1/   GLSHARE,GLOXCNT,GLIGAS 
C    /REGISF/ GLUSGF 
C    /RFORM1/ GLRFG,GLWIN 
C    /RFORM2/ GLRFGSYR 
C    /RUNLS1/ GLVMTU 
C    /SAVE08/ VA021 
C    /SPEED9/ GLHVMT,GLFVMT,GLSVMT 
C    /TEMPS/  GLTMIN,GLTMAX,GLTHRLY,GLDEL_TEM 
C    /VMXCOM/ GLVMTMIX 
C 
C  Output on return: 
C    modules: 
C     LEAKYDAT : ILIQUID_LEAKS 
C     DATABASE : DBNSCEN, DBSCFLAGS, DBSCTITLE 
C    common blocks: 
C    /AC_TIMES/ SUNRISE, SUNSET, PEAKSUN(2) 
C    /AC_DATA/ CLD_CVR 
C    /ALUHIN/ AC,ABSHUM 
C    /CHKLST/ CHKLIST 
C    /CITPAR/ SCNAME 
C    /CITRV1/ RVPBAS 
C    /DTGGAS/ DGPISY,DGPIEY,DGPHSF,FIDEFF,FIDFRC,IVDEFF,IVDFRC 
C    /EVAL/   MEVAL 
C    /EVAPACT/SAWFLAG,HSACTFLAG,DSACTFLAG,FPPRFLAG 
C    /FLAGS1/ OXYFLG,VMFLAG 
C    /FLAGS2/ ALHFLG 
C    /FLAGS5/ RVP_FLAG, SOAK_FLAG, TCF_FLAG,SPD_FLAG, 
C    /IMPAR7/ IREAD 
C    /IOUCOM/ IOUREP 
C    /IVTYPE/ MAXVEH 
C    /OXY1/   SHRMKT,OXYCNT,IGASHW 
C    /REGION/ IREJN 
C    /REGISF/ USRGSF 
C    /RFORM1/ RFGFLG,WINFLG 
C    /RUNLS1/ TLVMTU 
C    /SKDIST/ SOAKDIST 
C    /SPEED0/ AVGFLG 
C    /SPEED9/ HVMT,FVMT,SVMT 
C    /TEMPS/  TEMMIN,TEMMAX,TEMHRLY,DEL_TEM 
C    /VMXCOM/ VMTMIX 
C    /WEEKLY/ IWEEK 
C 
C    Mobile5 Control Section flags are set to default values. 
C 
C  Local variable dictionary: 
C 
C   Name     Type              Description 
C  ------    ----  ------------------------------------------------------- 
C  I          I    Loop control variable 
C  IH         I    Hourly loop control variable 
C  ITL        I    running loss travel length loop control variable 
C  IWK        I    weekend/weekday loop control variable 
C  IWK_DAY    I    weekday index (=1) 
C  IWK_END    I    weekend index (=2) 
C 
C  Notes: 
C 
C  INIFLG was added for Mobile6. 
C 
C  23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C            certification naming conventions. 
C 
C 24-Aug-98 @Dyntel-nh Request 0-632 Removed dead code related to special I/M 
C 
C ************************************************************************** 
C ************************************************************************** 
C  This subroutine initializes the flags that were set as we read the 
C  Control Section data in earlier versions of the Mobile model. 
C ************************************************************************** 
C ************************************************************************** 
C 
      USE LEAKYDAT, ONLY : ILIQUID_LEAKS 
      USE DATABASE, ONLY: DBNSCEN, DBSCFLAGS, DBSCTITLE 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'AC.I' 
      INCLUDE 'ALUHIN.I' 
      INCLUDE 'ATOX1.I'       ! TOXFLG 
      INCLUDE 'CHKLST.I' 
      INCLUDE 'CITPAR.I' 
      INCLUDE 'CITRV1.I' 
      INCLUDE 'DTGGAS.I' 
      INCLUDE 'EVAL.I' 
      INCLUDE 'EVAPACT.I' 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'IMPAR4.I' 
      INCLUDE 'IMPAR7.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'OXY1.I' 
      INCLUDE 'PART1.I' 
      INCLUDE 'PRINTC.I' 
      INCLUDE 'REGION.I' 
      INCLUDE 'REGISF.I' 
      INCLUDE 'RFORM1.I' 
      INCLUDE 'RFORM2.I' 
      INCLUDE 'RUNLS1.I' 
      INCLUDE 'SAVE08.I' 
      INCLUDE 'SKDIST.I' 
      INCLUDE 'SPEED0.I'      ! AVGFLG 
      INCLUDE 'SPEED9.I' 
      INCLUDE 'SULFUR.I' 
      INCLUDE 'TEMPS.I' 
      INCLUDE 'VDATA.I' 
      INCLUDE 'VMXCOM.I' 
      INCLUDE 'WEEKLY.I' 
C 
      INTEGER, INTENT(IN) :: NFILE 
      INTEGER, INTENT(IN) :: NRUN 
      INTEGER, INTENT(INOUT) :: NSCEN 
C 
      INTEGER   I 
      INTEGER, PARAMETER      :: IYES=1 
      CHARACTER*(*)  INPREC 
      CHARACTER LINEA*80 
      CHARACTER LINEB*80 
      CHARACTER*1, PARAMETER  :: SPACE = ' ' 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
C     Update the scenario counter and initialize for the next scenario. 
C 
   5  NSCEN = NSCEN + 1 
C 
C     Write a scenario header to the report file. 
C 
      DBSCTITLE = INPREC(21:) 
      WRITE(UNIT=LINEA,FMT=110) NFILE, NRUN, NSCEN 
  110 FORMAT('File ',I8,', Run ',I8,', Scenario ',I8,'.') 
      CALL LJCOMP(LINEA,LINEB) 
      IF( NSCEN .GT. 1 ) THEN 
        WRITE(IOUREP,112) NEWPAGE, INPREC(21:142), LINEB 
        IF(PARTFLG) 
     *  WRITE(IOUPM ,112) NEWPAGE, INPREC(21:142), LINEB 
        IF(TOXFLG) 
     *  WRITE(IOUAT ,112) NEWPAGE, INPREC(21:142), LINEB 
      ELSE 
        WRITE(IOUREP,112) SPACE, INPREC(21:142), LINEB 
        IF(PARTFLG) 
     *  WRITE(IOUPM ,112) SPACE, INPREC(21:142), LINEB 
        IF(TOXFLG) 
     *  WRITE(IOUAT ,112) SPACE, INPREC(21:142), LINEB 
      ENDIF 
  112 FORMAT(A, 
     &      /'* # # # # # # # # # # # # # # # # # # # # # # # # #', 
     &      /'*',A, /'* ',A, 
     &      /'* # # # # # # # # # # # # # # # # # # # # # # # # #') 
C 
C     The first three items on the checklist may be entered either 
C     globally or locally. The last four must be entered for each 
C     scenario. 
C 
      CHKLIST = GLCHECK 
C 
C     DATABASE (MODULE) 
C 
      DO I = 1,8 
        DBSCFLAGS(I) = 'N' 
      END DO 
      DBNSCEN = DBNSCEN + 1 
C 
C     Initialize parameters for the new scenario. 
C 
C     /DTGGAS/ 
C 
      DGPISY    = GLDGPSY 
      DGPIEY    = GLDGPEY 
      DGPHSF    = GLDGPHSF 
      FIDEFF    = GLFIDEF 
      FIDFRC    = GLFIDFC 
      IVDEFF    = GLIVDEF 
      IVDFRC    = GLIVDFC 
C 
C     /FLAGS1,VMXCOM/ 
C 
      VMFLAG = GLVMFLG 
      VMTUSR = GLVMTMIX 
C 
C     /FLAGS5/ 
C 
      RVP_FLAG  = 1 
      SOAK_FLAG = 1 
      TCF_FLAG  = 1 
      SPD_FLAG  = 1 
C 
C     /OXY1/ 
C 
      SHRMKT = GLSHARE 
C 
      OXYCNT(1)=GLOXCNT(1) 
      OXYCNT(2)=GLOXCNT(2) 
C 
      IGASHW = GLIGAS 
C 
C     /RFORM1/ RFGFLG (Reformulated Gasoline) 
C              WINFLG (In-house only) 
C 
      RFGFLG = GLRFG 
      WINFLG = GLWIN 
C 
C     Alternate I/M credit file flags. 
C 
      IF (IREAD .EQ. 0) SMRTL1=' ' 
C 
C     Initialize the ALH parameters. 
C 
      AC      = GLAC 
      ABSHUM  = GLAHUM 
      SUNRISE = GSUNRISE 
      SUNSET  = GSUNSET 
      PEAKSUN = GPEAKSUN 
      CLD_CVR = GCLD_CVR 
C 
C  Initialize the values from the Scenario Record. 
C 
C     /REGION/ 
C 
      IREJN = 1 
C 
C     /WEEKLY/ 
C 
      IWEEK = 1 
C 
C     /EVAL/ 
C 
      MEVAL = 1 
C 
C     The Local Area Parameter Record values may be entered either 
C     globally (LOCFLG=1 in Mobile5) or locally (LOCFLG=2 in Mobile5). 
C     In Mobile6, global values can be replaced selectively for individual 
C     scenarios. Start by copying the default (global) values to the local 
C     variables. 
C 
      SCNAME = ' ' 
      TEMMIN = GLTMIN 
      TEMMAX = GLTMAX 
      TEMP_DIST = .FALSE. 
      TEMHRLY = GLTHRLY 
      DEL_TEM = GLDEL_TEM 
      RVPBAS = GLRVPBAS 
C 
      OXYFLG=GLOXFLAG 
C 
      DSFLAG=GLDSFLAG 
      IF (DSFLAG.GT.1) THEN 
         USRGSF = GLUSGF 
      ENDIF 
C 
      DTGFLG = GLDTGFLG 
C 
C     /RUNLS1/ 
C 
      TLVMTU = GLVMTU 
C 
C     /SPEED0/ 
C 
      AVGFLG = 0 
C 
C     /SPEED9/ 
C 
      HVMT=GLHVMT 
      FVMT=GLFVMT 
      SVMT=GLSVMT 
C 
      FUEL_P = GLFUEL_P 
      S_PPM = GLS_PPM 
C 
C     module LEAKYDAT.FOR 
C     Initialize the liquid leak switch so that liquid 
C     leak emissions are included in the evaporative sums 
C 
      ILIQUID_LEAKS=IYES 
C 
C     Evap Sawtooth, Hot Soak and Diurnal Activity, Fixed P/P Rates 
C 
      SAWFLAG   = 1 
      HSACTFLAG = 0 
      DSACTFLAG = 0 
      FPPRFLAG  = 0 
C 
C     Scenario level initialization of soak distribution to MOBILE6 
C     default values (user input via subroutine GETSOAK) 
C 
      SOAKDIST = VA021 
C 
C     Initialize the PM and Air Toxic checklists 
C 
      PRTCHK = 1 
      TOXCHK = 1 
C 
      RETURN 
      END 
