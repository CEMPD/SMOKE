      SUBROUTINE PRCSCN(ICY,ICYOLD,INERR,RC) 
C 
C  PRCSCN performs all of the processing for an individual scenario after 
C  all of the data has been read in. The current scenario is complete 
C  when an END OF RUN labeled input record is found or a SCENARIO RECORD  
C  labeled input record signals the start of the next scenario. 
C 
C  Called by MOBILE. 
C 
C  Calls LOCAL,REGMOD,EFCALX,OUTPUT,RDT12C 
C 
C  Changes: (last change first) 
C 
C 14 Jan 02 @EPA-bag  Bug375   Added INERR to parameter list of AVGSPD 
C  03 Dec 01 AIR Task 27: Added call for calculation of CO2. 
C  27 Nov 01 AIR Task 28: Added OUTSS for spreadsheet output. 
C  11 Jul 01 AIR Task 15: Added call to ATCALX. 
C  07 Feb 01 AIR Task 11: Added call to PMCALX. 
C  14 Dec 01 @EPA-BAG        Added NO_IM_PGM 
C  19 Jul 01 @EPA-djb Wish235, Add call to AVGSPD. 
C  17 Jul 01 @EPA-bag Bug , Changed the subroutine name of LOCAL to FUELPAR 
C            to avoid conflict with array name in BDSPEED9.I 
C  19 Oct 00 @EPA-djb Task W2, Added call to DB_AGGR. 
C  06 Oct 00 @EPA-BG  Bug 155, Removed call to IDLCAL and IDLFLG 
C  25 Sep 00 @EPA-djb Bug 238, Removed call to DBWSCN (database output). 
C  13 Sep 00 @EPA-djb Bug 237, Removed conditional on call to REGMOD. 
C  07 Jul 00 AIR Task 05: Removed unused ICY from IDLCAL call. 
C  08 Jun 00 @EPA-ddj     bug#120, replaced call to OUTPUT with 
C                                 a call to OUTDSC (OUTPUT was 
C                                 deleted along with the 
C                                 spreadsheet options) 
C   1 May 00 @EPA-ddj     1-020, removed DATA statements 
C  11 Apr 00 @AIR-dk AIR-02: Replaced all occurences IOUASK with IOUIN 
C            and IOUOUT so that Lahey Standard for keyboard/monitor I/O 
C            is followed. 
C   4 Apr 00 @AIR-dk AIR-02: Removed ASTM class and changed REGION to 
C            ALTITUDE. 
C  24 Aug 99 @DynTel-MLA 2-684 Call database output routines. 
C  26 May 99 @DynTel-HXQ 2-698 Delete 'AVE SPEED' and 'OPERATING MOD' 
C  23 Mar 99 @DynTel-DDJ 2-696  Call to HRLOOP removed and ambient temperature 
C            removed from CHKLIST 
C  12 Mar 99 @DynTel-bag 2-664  IWEEK is now set by a labelled input 
C            and is passed through WEEKLY.I(BD16.FOR). It is initialized 
C            (default value) to weekday (IWEEK = 1) for each scenario in 
C            INISCN.FOR. Deleted NRUN from parameter list. 
C   5 Mar 99 @DynTel-bag 2-664  Added IWEEK and NRUN to parameters in 
C            call to EFCALX and NRUN in PRCSCN parameter list. Deleted 
C            call to EVAPHOUR. 
C   3 Mar 99 @DynTel-MLA 2-659  Changed IOUREP to IOUERR in error messages. 
C  11 Feb 99 @DynTel-MLA 2-663  Switched order of calls to RDT12C and LOCAL 
C            for Mobile6. 
C  18 Dec 98 @DynTel-ZK 2-000 Explicit Typing 
C  30 Oct 98 @DynTel-bag 2-677  Renamed HSFULL to EVAPHOUR. 
C  24 Aug 98 @DynTel-MLA 2-670  Added call to HSFULL routine. 
C 
C  20 Aug 98 @DynTel-MLA 2-659  Added the IWEEK variable and a call 
C            to the new HRLOOP routine. 
C 
C  Input on call: 
C 
C    Constants:  
C 
C    argument list: ICY,ICYOLD 
C 
C    common blocks: 
C    /CHKLST/ CHKLIST,GLCHECK,MAXREQ 
C    /EVAL/   MEVAL 
C    /FLAGS2/ IMFLAG 
C    /FLAGS4/ HCFLAG 
C    /IOUCOM/ IOUOUT 
C    /REGION/ IREJN 
C 
C  Output on return: 
C 
C    common blocks: 
C    /EVAL/   MEVOLD 
C    /REGION/ IROLD 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C 
C  Notes: 
C 
C  This subroutine is new for Mobile6. 
C 
      USE DATABASE, ONLY: DBFLAG, DBAGGR, DBSELPOL 
C 
      IMPLICIT NONE 
      INCLUDE 'ATOX1.I' 
      INCLUDE 'CHKLST.I' 
      INCLUDE 'EVAL.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'OPCNTL.I' 
      INCLUDE 'PART1.I' 
      INCLUDE 'REGION.I' 
      INCLUDE 'SPEED0.I'      ! AVGFLG 
C     
      INTEGER  ICYOLD, ICY 
C             
      INTEGER ERRFLAG 
      INTEGER IDX 
      INTEGER INERR 
      INTEGER INPRC 
      INTEGER IOS 
      INTEGER RC 
      CHARACTER*20, DIMENSION(MAXREQ)  :: CHKTEXT = 
     &  (/  'min and max temps.  ', 
     &      'fuel RVP.           ', 
     &      'calendar year.      '/) 
      INTEGER, SAVE :: NO_IM_PGM = 1 
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
      RC = 0 
C 
C     Make sure that values were entered for all required items. 
C     Unique return codes (1,2,4,8,16,...) are set up for each 
C     item in the checklist. A return code of 9 means that items 
C     1 and 4 were missing, a return code of 10 means that items 
C     2 and 4 were missing, etc. (The first two items -- Min/Max 
C     Temps. and Fuel RVP -- may be entered as global or as local 
C     values.) 
C 
C     NOTE: This list is only preliminary. We may choose in the 
C     future to require other items, or we may choose to set up 
C     defaults for some of these items so that they will no longer 
C     be required user inputs. 
C 
      INPRC = 0 
      ERRFLAG = 1 
      DO 20 IDX = 1,MAXREQ 
        IF (CHKLIST(IDX).EQ.0) THEN 
          WRITE(IOUERR,120,IOSTAT=IOS,ERR=90) CHKTEXT(IDX) 
  120     FORMAT('*** Error: no input data was found for ',A20, 
     &          /'    User input is required for this item.') 
          INERR = INERR + 1 
          INPRC = INPRC + ERRFLAG 
        ENDIF 
        ERRFLAG = 2 * ERRFLAG 
   20 CONTINUE 
C 
C     Calculate some scenario-specific values. 
C 
      CALL FUELPAR(ICY,INERR) 
C 
C     Set FVMT and SVMT arrays is an average speed has been specified. 
C 
      IF(AVGFLG==1) CALL AVGSPD(INERR) 
C 
C     Read the Tech1&2 I/M credit file. 
C 
      IF(IMFLAG/=NO_IM_PGM) CALL RDT12C(INERR,RC) 
C 
      IF(IOUOUT.EQ.6.OR.IOUOUT.EQ.9) THEN 
        WRITE(IOUOUT,145) 
  145   FORMAT(' Performing calculations.') 
      ENDIF 
C 
C     Perform the emission factor calculations. 
C 
C     For efficiency, the call to REGMOD could be restricted to when 
C     conditions in the Scenario have changed (i.e., ICY, IREJN, 
C     MEVAL, VMTMIX, etc.) 
C 
      CALL REGMOD(ICY,INERR) 
      ICYOLD=ICY 
      IROLD=IREJN 
      MEVOLD=MEVAL 
C 
      IF (INERR.EQ.0) THEN 
        CALL EFCALX(ICY,INERR) 
      ENDIF 
C 
C     Perform CO2 emission factor calculations, if necessary. 
C 
      IF(INERR.EQ.0.AND.DBSELPOL(4).EQ.2) CALL CO2CALX(ICY,INERR) 
C 
C     Perform PM emission factor calculations, if necessary. 
C 
      IF(INERR.EQ.0 .AND. PARTFLG) CALL PMCALX(ICY,INERR,RC) 
C 
C     Print the output report. 
C 
      IF(IOUOUT.EQ.6.OR.IOUOUT.EQ.9) WRITE(IOUOUT,150) 
  150 FORMAT(' Preparing output.') 
      IF (INERR.EQ.0 .AND. (.NOT.NODESC) ) THEN 
        CALL OUTDSC(ICY) 
      ENDIF 
C 
      IF (INERR.EQ.0.AND.DBFLAG.EQ.2.AND.DBAGGR) THEN 
        CALL DB_AGGR(ICY,INERR) 
      ENDIF 
C 
      IF (INERR.EQ.0) CALL OUTSS(ICY) 
C 
C     Perform Air Toxic emission factor calculations and produce output, 
C     when necessary. Must perform these calculations AFTER the output  
C     of criteria pollutants because of computational and array overlap. 
C      
      IF(TOXFLG) CALL ATCALX(ICY,INERR) 
C 
C     Return the RC describing input errors in preference to any 
C     other error return code. 
C 
      IF (INPRC.NE.0) RC = INPRC 
      GOTO 99 
C 
   90 RC = IOS 
      INERR = INERR + 1 
C 
   99 RETURN 
      END 
