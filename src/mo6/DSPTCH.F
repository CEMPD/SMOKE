      SUBROUTINE DSPTCH(INPREC,STDLAB,SECTION,NFILE,NSCEN,NRUN,ICY, 
     &                  INERR,RC) 
C 
C  This subroutine calls the appropriate subroutine to process a single 
C  Mobile6 input record, according to the label on the record. 
C 
C  Called by MOBILE. 
C 
C  Changes: (Last change first) 
C 
C 06 Nov 01 AIR Task 27: Added MPG ESTIMATES for CO2. 
C 31 Aug 01 AIR Task 18: Added ADDITIONAL HAPS. 
C 24 Aug 01 AIR Task 23b: Commented out PM Speed Correction. 
C 16 Aug 01 EPA Task 15-21: Changed RVP WAIVER to RVP OXY WAIVER. 
C 19 Jul 01 @EPA-djb Wish235, Add "Average Speed" command. 
C 02 Jul 01 AIR Task 15: Added Air Toxic OXYGENATE and RVP WAIVER commands. 
C 25 May 01 AIR Task 17: Added PARTICULATE EF command. 
C 23 May 01 AIR Task 17: Added PM speed correction. 
C 30 Apr 01 AIR Task 12: Added Tier2 PM Phase-in. 
C 30 Jan 01 AIR Task 11: Added all particulate and air toxic lebels. 
C  19 Jul 01 @EPA-djb Wish235, Add "Average Speed" command. 
C  25 May 01 @EPA-elg edited label for post 1996 exhaust I&M 
C  16 Apr 01 @EPA-elg Added Diagostic label POST96IM flag. 
C  19 Jan 01 @EPA-bag Added NO 2007 HDDV RULE 
C  20 Nov 00 @EPA-elg Added NO CLEAN AIR ACT label 
C 22 Nov 00 @EPA-bag bug301 added NSCEN to GETHVMT, GETFVMT arguments to allow 
c           scenario and run level user input 
C  21 Nov 00 @EPA-bag Bug303, removed redundant and confusing QUITER calls for 
c            message #'s 606 and 607 which concerned running loss trip fractions 
C  25 Oct 00 @EPA-djb Bug292, Changed text of SECTION error message. 
C  12 Oct 00 AIR Task 02: Bug Fix 279: Changed GETNAT, GETPRS, GETPRG, GETATP 
C            routines to free format input. 
C  02 Oct 00 AIR Task 07: bug Fix 269: added INERR processing for Global 
C            data errors. 
C  06 Oct 00 @EPA-BG  BUG155 removed NEW EVAP TEST, IDLE EFS, INDUSTRY AVERAGE FUEL, 
C            FTP/ FEDERAL TEST PROC, NO RUNNING LOSS CF, OBD PERCENT EQUIPPED 
C  15 Sep 00 AIR Task 08: Added EVAP I/M FRAC BEN and GLL REDUCTION. 
C  15 Sep 00 @EPA-elg Added new label T2_EVAP PHASE-IN for BUG 204. 
C  07 Sep 00 @EPA-elg Changed the VMT BY HOUR, FACILITY label to VMT BY HOUR 
C  07 Sep 00 @EPA-elg Changed the WINFLG label to SEASON. 
C  08 Aug 00 @EPA-elg Added new label to shut down evap calculations. 
C  08 Aug 00 @EPA-eg  bug216: Removed ICY from GETFUEL parameter list 
C  11 Aug 00 @EPA-bag Task 06: Removed WARNMSG subroutine redundancy with 
C            calls to QUITER for evaporative diagnostic labels 
C  17 Jul 00 AIR Task 06: Added NO EVAP SAWTOOTH, FIXED P/P RATES, HOT SOAK 
C            ACTIVITY, DIURN SOAK ACTIVITY. 
C  18 Jul 00 @EPA-bag BUG#54  Removed AMBIENT TEMPERATURE label and 
C                     TEMFLG. 
C  21 Jul 00 @EPA-RAG running loss; removed call to GETPFR user supplied 
C            pass/fail fractions, PFRATE's, for purge/pressure evap. tests 
C            deleted labels CALCULATE PASS/FAIL RATES and PASS/FAIL RATES 
C            added NO RESTING LOSS label 
C  06 Jul 00 @EPA-DDJ BUG#184 Added ICY to parameter list of GETFUEL. 
C  07 Jul 00 AIR Task 05: Removed unused INPREC and INERR from NOXPUL 
C            parameter list. 
C   5 Jul 00 @EPA-ddj Bug46 Routed all error, warning and comments within 
C            CASE statement to QUITER. 
C  26 Jun 00 AIR Task 04: Enabled (removed from comments) EXPAND EVAPORATIVE 
C            commands and warnings. Added messages for STAGE II REFUELING, 
C            WEEKDAY and WEEKEND TRIP LENGTH DISTRIBUTIONS. Left commented 
C            out (rather than deleted) the EXPAND EVPORATIVE EFS 0 REF and 
C            EXPAND EVPORATIVE EFS (GMS) processing (these commands may be 
C            reinstated later.) 
C  14 Jun 00 AIR Task 03: Added NO COLDCO STD and warning message. Removed 
C            PHASE 2 CO STANDARDS, COLD CO CERTIFICATION DATA, and  
C            NO CO OFFSET. 
C  19 May 00 AIR Task 03: Added NO WEATHERED RVP and warning message. 
C   8 Jun 00 @EPA-ddj, Bug120 Removed labels associated 
C            spreadsheet output. 
C  08 Jun 00 @EPA-elg Added warning message for NO SPEED and TEMP correction 
C                     factor cases. 
C  01 Jun 00 @EPA-elg Added label SPEED VMT  replaces VMT BY HO, SP, FAC. 
C  30 May 00 @EPA-DJB Added label for start distribution input 
C  23 May 00 @EPA-elg Changed label 'NO DEFEAT' to NO DEFEAT DEVICE' 
C  16 May 00 @EPA-ddj BUG87, Removed STARTS AS GM/ST label. 
C  16 May 00 @EPA-bag    Removed NO DESC OUTPUT label 
C   2 May 00 @EPA-ddj D-28, Changed REFORMULATED GASOLINE label to 
C            FUEL PROGRAM. 
C  28 Apr 00 @DynTel-ddj 1-020  Added NO SFTP SPEED and NO SFTP AC 
C  27 Apr 00 @EPA-elg Added NGV Label code. 
C  24 APR 00 AIR Task 02: Removed spurious INPREC from GETVMT parameter list. 
C                         Removed spurious IMPEF from GET94P parameter list. 
C  15 Apr 00 AIR Task 02: Changed EXPAND LDGT to EXPAND LDT. Disabled all 
C            EXPAND EVAP ... type statements so that descriptive evap 
C            output will not occur. 
C  10 Apr 00 Removed CALULATE HC, CO, NOX ONLY. Removed PRTFLG. Added NODESC. 
C  04 Apr 00 AIR Task 02: Change REGION to ALTITUDE. Removed ASTM class. 
C  07 Apr 00 @EPA-djb Added INPREC arguement to call to GETVMT. 
C            Removed dead code. 
C  18 Feb 00 @DynTel-rag 1-005, Added labels for evaporative liquid leaks. 
C  13 Jan 00 @DynTel-ddj 1-013  Added labels for AC calculation and 
C            removed ALHFLG CORRECTION FACTOR label. 
C  25 Feb 00 @DynTel-bag  corrected misspelling of vehicle registration label 
C  30 Ded 99 @DynTel-ddj 1-021, Added labels for Tier2 data entry. 
C            (i. e. NO TIER2, T2 PHASE-IN and T2 CERT) 
C  29 Nov 99 @DynTel-bag 1-006  Removed call to GETHFL.FOR which enabled hourly 
C            evaporative emissions output. The include statement for HOURLY.I and 
C            its associated variables were also deleted.  Added grams/start index, 
C            GPS, for start emissions. Revised HCFLAG  
C  16 Dec 99 @DynTel-ddj 1-010 Changed QUITER message for DISABLE 
C            1990 CAAA from "1" to "526".  Disabled the rest of the code 
C            associated with the disablement of CAAA. 
C  01 Nov 99 @Dyntel-nh Added DEFEAT subroutine calls 
C  06 Oct 99 @Dyntel-jwr 2-705 Disabled the non-periodic I/M programs 
C            COIM and RSD by commenting out the relevant 27 lines.  
C 
C  04 Oct 99 @Dyntel-jwr 2-705 Changed tech training credit input label 
C            to NO I/M TTC CREDITS reflecting change in default setting 
C  15 Sep 99 @DynTel-MLA 2-000  Del. GLSCNAME; Read SCTITLE from SCEN. REC. 
C   7 Sep 99 @DynTel-ddj 2-710  Added three new labels for FTP input 
C            option. 
C  24 Aug 99 @DynTel-MLA 2-684  Save the scenario title for database and 
C            descriptive output. 
C  24 Aug 99 @DynTel-ddj 2-000  Added full label for VMT BY FAC and 
C            VMT BY HOUR, FAC, i.e., VMT BY FACILITY and VMT BY HOUR, FACILITY 
C  27 Jul 99 @DynTel-MLA 2-694  Added label 94+ LDG EMISSION FACTORS. 
C  13 Jul 99 @DynTel-ddj 2-699  Added labels for user input of gasoline 
C            sulfur content. 
C   6 Jul 99 @DynTel-NH  2-705 Add labels for user input of 81-93 LDGV 
C            emissions numbers. 
C   1 Jul 99 @DynTel-HXQ 2-705 Add HDCGPM_FLAG. 
C  22 Jun 99 @DynTel-MLA 2-694  Removed labels for LEV PROGRAM DESCRIPTION, 
C            94+ LDG RUNNING EFS, and 94+ LDG START EFS. 
C   2 Jun 99 @DynTel-MLA 2-694  Added OBD input routines GETMIL, GETORR, 
C            GETOIC, GETORE, and GETOPE. 
C  28 May 99 @DynTel-ddj 2-700  Made modifications for individual test 
C            flags. (CO_FLAG, RVP_FLAG, SOAK_FLAG, TCF_FLAG, SPC_FLAG 
C            and RN_FLAG) 
C  26 May 99 @DynTel-HXQ 2-698 Delete code for operating mode and 
C                        average speed. 
C  21 May 99 @DynTel-ddj 2-700, Added label for FTP Flag. 
C  12 May 99 @DynTel-MLA 2-694  Added cases for OBD I/M. 
C  15 Apr 99 @DynTel-HQ 2-698 Add label for speed correction. 
C  12 Apr 99 @DynTel-ddj 2-000  Added labels to allow bypassing of CO 
C            offset correction calculation, RPV correction, temperature 
C            correction and soak time adjustment to temperature 
C            correction factors; call to GETFLAGS. 
C  08 Apr 99 @DynTel-bag 2-000  Removed call to GETAMB (ambient temperature) 
C  24 Mar 99 @DynTel-DDJ 2-696  Added labels 94 and 95 for CO phase 2 data 
C  19 Mar 99 @DynTel-MLA 2-701 Added INPREC to the parm. list for INISCN. 
C  10 Mar 99 @DynTel-bag 2-664  Added label for weekend index; this will 
C            allow user to specify that MOBILE6 calculations use weekend 
C            activity distributions/emissions data. Added label and call 
C            to subroutine to read and open output file for evaporative 
C            emissions. 
C   3 Mar 99  @DynTel-MLA 2-659  Send error messages to IOUERR. 
C  11 Feb 99  @DynTel-MLA 2-663  Added case for label STARTS AS GM/MI (63). 
C   5 Jan 99  @DynTel-MLA 2-663  Added new labels for Mobile6 I/M programs. 
C  30 Oct 98  @DynTel-bag 2-677 added a call to GETMMT for user input of 
C             hourly temperatures 
C  28 Sep 98  @DynTel-HQ 2-678 Conversion Factor 
C  17 Sept 98 @DynTel-zk 2-000 Explicit Typing 
C  27 Jul 98   DynTel@NH 0-632 
C 
C  Input on call: 
C 
C    modules : 
C     LEAKYDAT - ILIQUID_LEAKS 
C    Common blocks: 
C    /FLAGS2/     NEWFLG 
C    /IOUCOM/     IOUERR 
C    /OPCNTL/     OPTABL 
C    /UNIT_PARAM/ GPS 
C 
C    ICY     - The cal. year of eval. 
C    INPREC  - The last input record read from the main input file. 
C    LABTYPE - An integer describing the input record label type. 
C 
C  Output on return: 
C 
C    Common blocks: 
C    /CITPAR/  SCNAME 
C    /FLAGS1/  DTGFLG,TAMFLG 
C    /FLAGS2/  NEWFLG 
C    /FLAGS3/  RLFLAG, CAAFLG, POST96IM 
C    /FLAGS4/  HCFLAG,NMHFLG 
C    /FLAGS5/  RVP_FLAG, SOAK_FLAG, TCF_FLAG, SPD_FLAG, 
C              WXRVP_FLAG, EVAP_FLAG 
C    /IMPAR6/  IFILET 
C    /IMPAR7/  TNAME 
C    /OPCNTL/  OPTABL, NODESC 
C    /PART5/   PMSPDFLG 
C    /RESTLS1/ NO_RSTLOS,NO_EVAPCALC 
C    /WEEKLY/  IWEEK 
C 
C    INERR - A running tally of the number of errors. 
C    RC    - The return code (zero if there are no errors). 
C          Return codes less than zero mean end-of-file conditions. 
C          Return codes greater than zero mean various error conditions. 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C   IWK     I     Determines whether weekday(IWK = 1) or weekend(IWK = 2) 
C                 trip length distributions are read by call to WTLRNL 
C  Notes: 
C 
C  DSPTCH was added for MOBILE6. It calls many of the subroutines that 
C  were called form ONESEC, CONSEC, and PARSEC in earlier versions of 
C  the Mobile model. 
C 
      USE LEAKYDAT, ONLY : ILIQUID_LEAKS 
      USE DATABASE, ONLY: DBSCTITLE 
C 
      IMPLICIT NONE 
      INCLUDE 'CHKLST.I' 
      INCLUDE 'CITPAR.I' 
      INCLUDE 'CITRV1.I' 
      INCLUDE 'EVAPACT.I' 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'IDLE.I' 
      INCLUDE 'IMPAR6.I' 
      INCLUDE 'IMPAR7.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'NGV.I' 
      INCLUDE 'OPCNTL.I' 
      INCLUDE 'PART5.I' 
      INCLUDE 'RESTLS1.I' 
      INCLUDE 'RFORM1.I' 
      INCLUDE 'RSDINF.I' 
      INCLUDE 'SPEED0.I'      ! AVGFLG 
      INCLUDE 'UNIT_PARAM.I' 
      INCLUDE 'WEEKLY.I' 
C 
      CHARACTER      BUFFER1*150 
      CHARACTER      BUFFER2*150 
      INTEGER        ICY 
      INTEGER        INERR 
      CHARACTER      INPREC*(*) 
      INTEGER        NFILE 
      INTEGER        NRUN 
      INTEGER        NSCEN 
      INTEGER        RC 
      INTEGER        SECTION 
      CHARACTER      STDLAB*40 
      INTEGER        IOS 
      INTEGER        IWK 
      INTEGER, PARAMETER :: IYES=1 
      INTEGER, PARAMETER :: INO=0 
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
C     Initialize. 
C 
      RC = 0 
C 
C     The section 0 labels belong in the header. They are processed in 
C     PROCHDR. This routine should only be called to process the run 
C     data section of the input file. 
C 
      IF (SECTION.EQ.0) THEN 
        WRITE(IOUERR,119,IOSTAT=IOS,ERR=95) INPREC 
  119   FORMAT(/'*** Error: the following command is only', 
     &            ' allowed in the header data section:',/A,/) 
        INERR = INERR + 1 
        GOTO 99 
      ELSE IF (SECTION.EQ.1.AND.NSCEN.GT.0) THEN 
        WRITE(IOUERR,120,IOSTAT=IOS,ERR=95) INPREC 
  120   FORMAT('*** Error: the following command is not', 
     &           ' allowed in the scenario',  
     &          /T12,'data section of the input file:',/A,/) 
        INERR = INERR + 1    
      ELSEIF (STDLAB.EQ.'END OF RUN'.AND.NSCEN.EQ.0) THEN 
        CALL QUITER(0.0,IOS,547,INERR) 
      ELSEIF (SECTION.EQ.3.AND.NSCEN.EQ.0) THEN 
        WRITE(IOUERR,122,IOSTAT=IOS,ERR=95) INPREC 
  122   FORMAT('*** Error: the following command is not', 
     &         ' allowed in the run data', 
     &         /T12,'section preceeding the scenario inputs:',/A,/) 
        INERR = INERR + 1 
      ENDIF 
C 
C     Data items in the first group must appear in the run data 
C     section. They cannot appear as scenario inputs. 
C 
      PROCLAB: SELECT CASE (STDLAB) 
        CASE ('TAMPERING RATES') 
          CALL GETTAM(INPREC,INERR,RC) 
        CASE ('NO TAMPERING') 
          TAMFLG = 0 
        CASE ('REGISTRATION DISTRIBUTIONS', 
     &        'VEHICLE REGISTRATIONS') 
          CALL GETREG(INPREC,INERR,RC) 
        CASE ('MILEAGE ACCUMULATION RATES') 
          CALL GETMAR(INPREC,INERR,RC) 
        CASE ('BASIC EFS') 
          CALL GETERQ(INPREC,INERR,RC) 
        CASE ('SPEED COEFFICIENT', 
     &        'SPD COEFFICIENT') 
          CALL GETSPCO(INPREC,INERR,RC) 
        CASE ('SCF BY FWAY AND AC') 
          CALL GETSCF(INPREC,INERR,RC) 
        CASE ('HD CONVERSION FACTORS') 
          CALL GETCONV(INPREC,INERR,RC) 
        CASE ('DISABLE 1990 CAAA') 
          CALL QUITER(5.0, NEWFLG, 526, INERR) 
        CASE ('NO HDCGPM') 
          HDCGPM_FLAG = 0 
        CASE ('1981-93 LDG EFS') 
          CALL GET_8193(INPREC,INERR,RC) 
C 
        CASE ('NGV FRACTION') 
          NGVFLG1 = 1 
          CALL GETNGVFR(INPREC,INERR,RC) 
        CASE ('NGV EF') 
          NGVFLG2 = 1 
          CALL GETNGVEF(INPREC,INERR,RC) 
C 
        CASE ('PARTICLE SIZE') 
          CALL GETPSC(INPREC,INERR,RC) 
        CASE ('DIESEL SULFUR') 
          CALL GETD_PPM(INPREC,INERR,RC) 
        CASE ('PARTICULATE EF') 
          CALL GETPMEF(INPREC,INERR,RC) 
        CASE ('ADDITIONAL HAPS') 
          CALL GETATHAPS(INPREC,INERR,RC) 
C 
        CASE ('E200') 
          CALL GETE200(INPREC,INERR,RC) 
        CASE ('E300') 
          CALL GETE300(INPREC,INERR,RC) 
        CASE ('GAS AROMATIC%') 
          CALL GETAROM(INPREC,INERR,RC) 
        CASE ('GAS OLEFIN%') 
          CALL GETOLEF(INPREC,INERR,RC) 
        CASE ('GAS BENZENE%') 
          CALL GETBENZ(INPREC,INERR,RC) 
        CASE ('OXYGENATE') 
          CALL GETATOXY(INPREC,INERR,RC) 
        CASE ('RVP OXY WAIVER') 
          CALL GETATWAV(INPREC,INERR,RC) 
C 
        CASE ('MPG ESTIMATES') 
          CALL GETMPG(INPREC,INERR,RC) 
C 
C       Start of I/M section. (Any changes in the I/M section 
C       must also be applied to subroutine GETXIM.) 
C 
        CASE ('I/M PROGRAM', 
     &        'I/M PGM') 
          CALL GETIPR(INPREC,INERR,RC) 
        CASE ('I/M MODEL YEARS') 
          CALL GETIMY(INPREC,INERR,RC) 
        CASE ('I/M EFFECTIVENESS') 
          CALL GETEFR(INPREC,INERR,RC) 
        CASE ('I/M CREDIT FILE') 
C         Save file name in common IMPAR7. 
ccs          CALL UCCOMP(INPREC(22:),TNAME) 
          CALL LJCOMP(INPREC(22:),TNAME)
          IFILET = 2 
        CASE ('I/M VEHICLES') 
          CALL GETIMV(INPREC,INERR,RC) 
        CASE ('I/M STRINGENCY') 
          CALL GETSGC(INPREC,INERR,RC) 
        CASE ('I/M COMPLIANCE') 
          CALL GETCPL(INPREC,INERR,RC) 
        CASE ('I/M WAIVER RATES') 
          CALL GETWAV(INPREC,INERR,RC) 
        CASE ('I/M WAIVER REDUCTION') 
          CALL GETWRR(INPREC,INERR,RC) 
        CASE ('I/M CUTPOINTS') 
          CALL GETCUT(INPREC,INERR,RC) 
        CASE ('I/M EXEMPTION AGE') 
          CALL GETXMP(INPREC,INERR,RC) 
        CASE ('I/M GRACE PERIOD') 
          CALL GETGRA(INPREC,INERR,RC) 
        CASE ('NO I/M TTC CREDITS') 
          CALL GETTTC(INPREC,INERR,RC) 
        CASE ('I/M TTC EFFECTIVENESS') 
          CALL GETTEF(INPREC,INERR,RC) 
        CASE ('OBD WORKING MIL') 
          CALL GETMIL(INPREC,INERR,RC) 
          CALL QUITER(0.0,0,693,INERR) 
        CASE ('OBD RESPONSE') 
          CALL GETORR(INPREC,INERR,RC) 
          CALL QUITER(0.0,0,691,INERR) 
        CASE ('OBD I/M COVERAGE') 
          CALL GETOIC(INPREC,INERR,RC) 
          CALL QUITER(0.0,0,692,INERR) 
        CASE ('OBD REPAIRED EFS') 
          CALL GETORE(INPREC,INERR,RC) 
        CASE ('EVAP I/M FRAC BEN') 
          CALL GETFBEN(INPREC,INERR,RC) 
        CASE ('I/M DESCRIPTION FILE') 
          CALL GETXIM(INPREC,INERR,RC) 
C 
C       End of I/M section. (Any changes in the I/M section 
C       must also be applied to subroutine GETXIM.) 
C 
        CASE ('ANTI-TAMPERING PROGRAM', 
     &        'ANTI-TAMPERING PGM', 
     &        'GET ATP') 
          CALL GETATP(INPREC,INERR,RC) 
        CASE ('PRESSURE CHECK PARAMETERS', 
     &        'PRESSURE CHK PARAMETERS') 
          CALL GETPRS(INPREC,INERR,RC) 
        CASE ('PURGE CHECK PARAMETERS', 
     &        'PURGE CHK PARAMETERS') 
          CALL GETPRG(INPREC,INERR,RC) 
        CASE ('STAGE II REFUELING', 
     &        'STG II REFUELING', 
     &        'STAGE 2 REFUELING', 
     &        'STG 2 REFUELING') 
          CALL GETNAT(INPREC,INERR,RC) 
          CALL QUITER(0.0,0,601,INERR) 
        CASE ('EXPAND LDT EFS') 
          OPTABL(1) = 2 
        CASE ('EXPAND HDGV EFS') 
          OPTABL(2) = 2 
        CASE ('EXPAND HDDV EFS') 
          OPTABL(3) = 2 
        CASE ('EXPAND BUS EFS') 
          OPTABL(4) = 2 
        CASE ('EXPRESS HC AS THC') 
          NMHFLG = 1 
        CASE ('EXPRESS HC AS NMHC') 
          NMHFLG = 2 
        CASE ('EXPRESS HC AS VOC') 
          NMHFLG = 3 
        CASE ('EXPRESS HC AS TOG') 
          NMHFLG = 4 
        CASE ('EXPRESS HC AS NMOG') 
          NMHFLG = 5 
        CASE ('EXPAND EVAPORATIVE') 
          IF (EVAP_FLAG.EQ.1) THEN 
          HCFLAG = 2 
          ELSE 
          CALL QUITER(0.0,0,602,INERR) 
          ENDIF 
        CASE ('NO REFUELING') 
          RLFLAG = 5 
          CALL QUITER(0.0,0,603,INERR) 
        CASE ('NO COLDCO STD') 
          COLDFG = 1 
          CALL QUITER(0.0,0,604,INERR) 
        CASE ('NO WEATHERED RVP') 
          WXRVP_FLAG = 0 
          CALL QUITER(0.0,0,699,INERR) 
        CASE ('94+ LDG IMPLEMENTATION') 
          CALL GET94P(INPREC,INERR,RC) 
        CASE ('NO TIER2') 
          CALL NO_TIER2() 
        CASE ('T2 EXH PHASE-IN') 
          CALL GETT2FRC(INPREC,INERR,RC) 
        CASE ('T2 EVAP PHASE-IN') 
          CALL GETT2EVP(INPREC,INERR,RC) 
        CASE ('T2 CERT') 
          CALL GETT2CRT(INPREC,INERR,RC) 
        CASE ('EXPAND EXHAUST') 
          EXPEXH = .TRUE. 
        CASE ('GLL REDUCTION') 
          CALL GETGLLRDC(INPREC,INERR,RC) 
          CALL QUITER(0.0,0,694,INERR) 
        CASE ('T2 PM PHASE-IN') 
          CALL GETT2PMFRC(INPREC,INERR,RC) 
        CASE ('IDLE PM EMISSIONS') 
          IDLEPMFLG = .TRUE. 
C 
C       The data items in the next group can appear in the run 
C       data section or in the scenario data section. That is, they 
C       can appear either before or after a scenario record has 
C       appeared in the input file. 
C 
        CASE ('SCENARIO RECORD') 
          DBSCTITLE = INPREC(22:) 
          SCNAME = INPREC(22:) 
          CALL INISCN(INPREC,NFILE,NSCEN,NRUN) 
        CASE ('WEEK DAY TRIP LENGTH DISTRIBUTIONS') 
          IWK = 1 
          CALL WTLRNL(IWK,INPREC,NSCEN,INERR,RC) 
        CASE ('WEEK END TRIP LENGTH DISTRIBUTIONS') 
          IWK = 2 
          CALL WTLRNL(IWK,INPREC,NSCEN,INERR,RC) 
        CASE ('VMT FRACTIONS', 
     &        'VMT MIX') 
          CALL GETVMT(NSCEN,INERR,RC) 
        CASE ('VMT BY FACILITY', 
     &        'VMT BY FAC') 
          CALL GETFVMT(INPREC,NSCEN,INERR,RC) 
        CASE ('VMT BY HOUR') 
          CALL GETHVMT(INPREC,NSCEN,INERR,RC) 
        CASE ('AVERAGE SPEED') 
          CALL GETAVG(INPREC,INERR,RC) 
          AVGFLG = 1 
        CASE ('SPEED VMT') 
          CALL GETSVMT(INPREC,NSCEN,INERR,RC) 
        CASE ('CLOUD COVER') 
          RC = 1 
          CALL GETAC(INPREC,NSCEN,INERR,RC) 
        CASE ('SUNRISE/SUNSET') 
          RC = 2 
          CALL GETAC(INPREC,NSCEN,INERR,RC) 
        CASE ('PEAK SUN') 
          RC = 3 
          CALL GETAC(INPREC,NSCEN,INERR,RC) 
        CASE ('NO AC') 
          RC = 4 
          CALL GETAC(INPREC,NSCEN,INERR,RC) 
        CASE ('ABSOLUTE HUMIDITY') 
          RC = 5 
          CALL GETAC(INPREC,NSCEN,INERR,RC) 
        CASE ('AC FRACTION') 
          RC = 6 
          CALL GETAC(INPREC,NSCEN,INERR,RC) 
        CASE ('DETERGENT GASOLINE') 
          CALL GETDTG(NSCEN,INERR,RC) 
        CASE ('NO DETERGENT GAS') 
          IF (DTGFLG.EQ.2) THEN 
            RC = 1 
          ENDIF 
          DTGFLG = 0 
          IF (NSCEN.EQ.0) THEN 
            GLDTGFLG = 0 
          ENDIF 
        CASE ('OXYGENATED FUELS') 
          CALL GETOXY(INPREC,NSCEN,INERR,RC) 
        CASE ('DIESEL FRACTIONS', 
     &        'DIESEL SALES') 
          CALL GETDSF(NSCEN,INERR,RC) 
        CASE ('MIN/MAX TEMPERATURES') 
          CALL GETMMT(INPREC,NSCEN,INERR,RC,0) 
          CHKLIST(1) = 1 
          IF (NSCEN.EQ.0) THEN 
            GLCHECK(1) = 1 
          ENDIF 
        CASE ('HOURLY TEMPERATURES') 
          CALL GETMMT(INPREC,NSCEN,INERR,RC,1) 
          CHKLIST(1) = 1 
          IF (NSCEN.EQ.0) THEN 
            GLCHECK(1) = 1 
          ENDIF 
        CASE ('FUEL RVP') 
          CALL GETRVP(INPREC,NSCEN,INERR,RC) 
          CHKLIST(2) = 1 
          IF (NSCEN.EQ.0) THEN 
            GLCHECK(2) = 1 
          ENDIF 
        CASE ('FUEL PROGRAM') 
          CALL GETFUEL(INPREC,NSCEN,INERR,RC) 
        CASE ('SEASON') 
          CALL GETWIN(INPREC,NSCEN,INERR,RC) 
        CASE('FTP MILES PER DAY') 
          FTP_INPUT_FLG = 1 
        CASE('STARTS PER DAY') 
          CALL GETSPD(INPREC,INERR,RC) 
        CASE('START DIST') 
          CALL GETSDIST(INPREC,INERR,RC) 
        CASE('SOAK DISTRIBUTION') 
          CALL GETSOAK(INPREC,INERR,RC) 
C        
C       The following data items can only appear inside of the Scenario 
C       Data Section. That is, there must be at least one Scenario record 
C       preceeding any of these records (including the End of Run record). 
C       These records are the Region, Calendar Year, Month, Average Speed, 
C       Ambient Temperature, Operating Mode Fractions, User-supplied 
C       Pass/Fail Rates, the Calculate Pass/Fail record, and the End of  
C       Run record. 
C 
        CASE ('ALTITUDE') 
          CALL GETRJN(INPREC,INERR,RC) 
        CASE ('CALENDAR YEAR') 
          CALL GETCYR(INPREC,ICY,INERR,RC) 
          CHKLIST(3) = 1 
        CASE ('EVALUATION MONTH') 
          CALL GETMTH(INPREC,INERR,RC) 
        CASE ('NO RESTING LOSS') 
          NO_RSTLOS = IYES 
          CALL QUITER(0.0,0,696,INERR) 
        CASE ('NO EVAP CALC') 
          NO_EVAPCALC = IYES 
          CALL QUITER(0.0,0,695,INERR) 
        CASE ('WEEKEND VEHICLE USAGE') 
          IWEEK = 2 
        CASE ('T1-LEV EFS') 
          CALL GET94E(INPREC,INERR,RC) 
        CASE ('NO RVP CORRECTION') 
          RVP_FLAG = 0 
          CALL QUITER(0.0,0,608,INERR) 
        CASE ('NO SOAK CORRECTION') 
          CALL QUITER(0.0,0,609,INERR) 
          SOAK_FLAG = 0 
        CASE ('NO TEMPERATURE CF', 
     &        'NO TCF', 
     &        'NO TEMP CORRECTION') 
          TCF_FLAG = 0 
          CALL QUITER(0.0,0,610,INERR) 
        CASE ('NO SPEED CORRECTION') 
          SPD_FLAG = 0 
          CALL QUITER(0.0,0,611,INERR) 
C        CASE ('PM SPEED CORRECTION') 
C          PMSPDFLG = 1 
C         CALL QUITER(0.0,0,622,INERR) 
        CASE ('NO SFTP SPEED') 
          SFTP_SP_FLAG = 0 
          CALL QUITER(0.0,0,612,INERR) 
        CASE ('NO SFTP AC') 
          SFTP_AC_FLAG = 0 
          CALL QUITER(0.0,0,613,INERR) 
        CASE('SULFUR CONTENT', 
     &       'SULFUR CORRECTION FACTOR', 
     &       'GASOLINE SULFUR') 
C         *** Input gasoline sulfur content (ppm) and 
C         *** User Supplied sulfur Correction Factors (USSCF) 
C         *** if applicable 
          CALL GET_SULFUR(INPREC,INERR,RC) 
C 
        CASE('NO NOX PULL AHEAD') 
            CALL NOXPUL(RC) 
        CASE ('REBUILD EFFECTS','NO REBUILD','NO DEFEAT DEVICE') 
            CALL RBEFFCTS(INPREC,INERR,RC) 
C 
        CASE('NO LIQUID LEAKS') 
          ILIQUID_LEAKS = INO 
        CASE ('NO EVP SAWTOOTH') 
          SAWFLAG = 0 
          CALL QUITER(0.0,0,697,INERR) 
        CASE ('FIXED PP RATES') 
          CALL GETFPPR(INPREC,INERR,RC) 
          IF(RC.EQ.0) CALL QUITER(0.0,0,698,INERR) 
        CASE ('HOT SOAK ACTIVITY') 
          CALL GETHSACT(INPREC,INERR,RC) 
        CASE ('DIURN SOAK ACTIVITY') 
          CALL GETDSACT(INPREC,INERR,RC) 
        CASE ('NO CLEAN AIR ACT') 
          CAAFLG = 1 
          WRITE(IOUREP,123) 
  123     FORMAT(/' Results Reflect No Clean Air Act of 1990',/) 
        CASE ('NO 2007 HDDV RULE') 
          HDR2007 = INO 
          CALL QUITER(0.0,0,621,INERR) 
        CASE ('1996+ EXHAUST IM') 
          POST96IM = 1 
          WRITE(IOUREP,124) 
  124     FORMAT(/' Exhaust I/M Requested on 1996+ MYs',/) 
        CASE ('END OF RUN', 
     &        '0000000000') 
          RC = 0 
        CASE DEFAULT 
C          *** Unrecognized Label Type! 
          WRITE(IOUERR,180,IOSTAT=IOS,ERR=95) INPREC(1:19) 
  180     FORMAT('*** ERROR 1 in DSPTCH: The label ', A19, 
     &           ' is unrecognized.') 
          RC = 1 
          INERR = INERR + 1 
      END SELECT PROCLAB 
C 
      IF (RC.NE.0) THEN 
        WRITE(BUFFER1,190,IOSTAT=IOS,ERR=95) RC 
  190   FORMAT('*** ERROR ',I11,' processing the following input', 
     &         ' record.') 
        CALL LJCOMP(BUFFER1,BUFFER2) 
        WRITE (IOUERR,191,IOSTAT=IOS,ERR=95) BUFFER2,INPREC 
  191   FORMAT(A,/'*** ',A/) 
      ENDIF 
C 
      RETURN 
C 
   95 RC = IOS 
      INERR = INERR + 1 
      WRITE(IOUOUT,195) IOS 
  195 FORMAT( 
     &  '********************************************************' 
     & /'*** Dispatch error writing to the report file: ',I8, 
     & /'********************************************************'/) 
C 
   99 RETURN 
      END 
