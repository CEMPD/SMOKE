      SUBROUTINE INILAB()
C
C  Initialize the LABLIST array and set the IDNUM field in the array.
C
C  Called by MOBILE.
C
C  Calls ...
C
C  Input on call:
C
C  Output on return:
C
C  Changes : (Last change first)
C
C  02 Sep 09 @EPA-bag bug22 Added ROADFE; allows ROADFE to be external data
C            file label as an alternate to MPG ESTIMATES
C  02 Jul 02 AIR Task 36: Added RELATIVE HUMIDITY and BAROMETRIC PRES.
C  20 Mar 02 AIR Task 32: Added NO HIGH EMITTERS, OBD EXH RESPONSE,
C            OBD EVAP RESPONSE, removed OBD RESPONSE
C  02 May 02 @EPA-bag Added HOT SOAK LQD LEAK
C  27 Nov 01 AIR Task 27: Changed SPREADSHEET to SPREADSHEET OUTPUT.
C  06 Nov 01 AIR Task 27: Add MPG ESTIMATES for CO2.
C  28 Sep 01 AIR Task 25: Changed TOXICS SPREADSHEET to SPREADSHEET to
C            accomodate both AT and PM outputs.
C  31 Aug 01 AIR Task 18: Added ADDITONAL HAPS (removed PM SPEED CORRECTION).
C  30 Aug 01 AIR Task 22: Added TOXICS SPREADSHEET 
C  16 Aug 01 AIR Task 15-21: Changed RVP WAIVER to RVP OXY WAIVER.
C  02 Jul 01 AIR Task 15: Added Air Toxic OXYGENATE and RVP WAIVER commands.
C  23 May 01 AIR Task 17: Added PM SPEED CORRECTION and PARTICULATE EF labels.
C  10 May 01 AIR Task 11/5: Added PM Idle.
C  30 Apr 01 AIR Task 12: Added T2 PM PHASE-IN.
C  17 Apr 01 AIR Task 11: Changed most particulate related input commands to
C            Scenario Level Only status.
C  30 Jan 00 AIR Task 11: Added all PARTICULATE and AIR TOXIC header, run
C            and scenario level labels.
C  19 Jul 01 @EPA-djb Wish235, Add Average Speed command (173).
C  16 Apr 01 @EPA-elg Added POST96IM Flag to control diagnostic I/M command.
C  19 Jan 01 @EPA-bag Added NO 2007 HDDV RULE
C  20 Nov 00 @EPA-elg Added NO CLEAN AIR ACT label
C  29 Nov 00 @EPA-BG  renamed database age to database ages
C  24 Oct 00 @EPA-BG  New labels for database age, model year and hour selection
C  19 Oct 00 @EPA-djb Task W2, Added AGGREGATED OUTPUT command.
C  27 Oct 00 @EPA-BG  BUG293 OBD REPAIRED EFS put back
C  06 Oct 00 @EPA-BG  Bug 155, Removed NEW EVAPORATIVE TEST, IDLE EFS,
C            INDUSTRY AVERAGE FUEL, STARTS AS GM/ST, FEDERAL TEST PROC,
c            NO RUNNING LOSS CF, OBD PERCENT EQUIPPED
C  15 Sep 00 AIR Task 08: Added EVAP I/M FRAC BEN and GLL Reduction.
C  07 Sep 00 @EPA-elg bug204 Added and changed labels to distinguish between Tier2
C            exhaust and evaporative emissions phase in labels
C  07 Sep 00 @EPA-elg Changed the VMT BY HOUR, FACILITY label to VMT BY HOUR
C  07 Sep 00 @EPA-elg Changed the WINFLG label to SEASON.
C  08 Aug 00 @EPA-elg Added new label to shut down evap calculations.
C  11 Jul 00 AIR Task 06: Added NO EVP SAWTOOTH, FIXED P/P RATES, HOT SOAK
C            ACTIVITY, and DIURN SOAK ACTIVITY.
C  18 Jul 00 @EPA-bag BUG#54  Removed AMBIENT TEMPERATURE label.
C  21 Jul 00 @EPA-RAG running loss; removed PASS/FAIL RATES label for
C            user suplied pass fail fractions, PFRATE's, for purge/pressure
C            evap. tests; remove CALCULATE PASS/FAIL RATES label; added
C            NO RESTING LOSS label
C  27 Jun 00 Added NO EVAPORATIVE EFS and renamed EXPAND EVPORATIVE EFS
C            to EXPAND EVPORATIVE; Commented out (rather than deleted) the
C            EXPAND EVPORATIVE EFS 0 REF and EXPAND EVPORATIVE EFS (GMS).
C  14 Jun 00 AIR Task 03: Added NO COLDCO STD. Removed PHASE 2 CO
C            STANDARDS, COLD CO CERTIFICATION DATA, and NO CO OFFSET.
C  19 May 00 AIR Task 03: Added NO WEATHERED RVP
C   8 Jun 00 @EPA-ddj, Bug120 Removed labels associated
C            spreadsheet output.
C  01 Jun 00 @EPA-elg Added label SPEED VMT  replaces VMT BY HO, SP, FAC.
C  30 May 00 @EPA-djb Added START DIST (180)
C  23 May 00 @EPA-elg Changed label 'NO DEFEAT' to NO DEFEAT DEVICE'
C  16 May 00 @EPA-bag   fixed label problem with NO DESC OUTPUT label;
C                      CHKLAB needs a match for each set of characters
C                      between the shortened and long version of a label;
C                      the shortened version is a lower character limit
C                      which cannot have a blank as the lower character
C                      limit for any of the character sets in a label
C                      (character sets are separated by a blank or space)
C  05 MAY 00 @EPA-hmc Added WITH FIELDNAMES and DAILY OUTPUT  labels
C   2 May 00 @EPA-ddj BUG7, Removed SCENARIO TITLE label.  Changed
C            REFORMULATED GASOLINE label to FUEL PROGRAM.
C  28 Apr 00 @EPA-elg Added two SFTP labels.
C  17 Apr 00 @EPA-elg Added NGV Label code.
C  13 Apr 00 AIR Task 02: Changed EXPAND LDGT to EXPAND LDT.
C  11 Apr 00 AIR Task 02: Modified code to ANSI Standard FORTRAN90, which
C            only allows up to 19 continuations in fixed format.
C  10 Apr 00 AIR Task 02: Removed DATABASE POLLUTANTS and all abbreviates
C            and installed POLLUTANTS and POL.
C  04 Apr 00 AIR Task 02: Changed REGION and REG to ALTITUDE. Removed ASTM
C            class. Added NO DESC OUTPUT and EXPAND EXHAUST.
C  20 Feb 00 @DynTel-bag 1-005 Added NO LIQUID LEAKS label
C  13 Jan 00 @DynTel-ddj 1-013  Removed ALHFLG CORRECTION FACTOR label
C            and added AC labels.
C  25 Feb 00 @DynTel-bag  corrected misspelling of vehicle registration label
C  30 Dec 99 @DynTel-ddj 1-021 Added NO TIER2, T2 CERT and
C            T2 PHASE-IN labels.
C  29 Nov 99 @DynTel-bag 1-006  Removed label for evaporative hourly
C            emissions in EVP_HRLY.FOR
C  01 Nov 99 @Dyntel-nh 1-004 Added 4 DEFEAT related labels.
C  04 Oct 99 @Dyntel-jwr 2-705 Changed imput tech training credit label 
C            to NO I/M TTC CREDITS to reflect new default setting for TTC
C  19 Sep 99 @DynTel-MLA 2-684  Added labels for new input methods and
C            database output
C   7 Sep 99 @DynTel-ddj 2-710  Added three new labels for the FTP
C            input option.
C  27 Jul 99 @DynTel-MLA 2-000  Removed INILAB from the module LABELS and
C            made it an ordinary external procedure (not a module proc.).
C  27 Jul 99 @DynTel-MLA 2-694  Added label 94+ LDG EMISSION FACTORS.
C  13 Jul 99 @DynTel-ddj 2-699  Added labels for user input of gasoline
C            sulfur content.
C   6 Jul 99 @DynTel-NH  2-705 Add labels for user input of 81-93 LDGV
C            emissions numbers.
C   1 Jul 99 @DynTel-HXQ 2-705 Add HDCGPM_FLAG.
C  22 Jun 99 @DynTel-MLA 2-694  Removed labels for 94+ START and 94+
C            RUNNING EFS and for LEV PROGRAM DESCRIPTION.
C  29 May 99 @DynTel-ddj 2-700  Modified TEST FLAGS so that each flag
C            would have its own label.
C  26 May 99 @DynTel-HXQ 2-698 Delete code for operating mode and
C                        average speed
C  20 May 99 @DynTel-ddj 2-700  Added label for FTP FLAG.
C  11 May 99 @DynTel-MLA 2-674  Created the module LABELS which contains
C            the CHKLAB subroutine and the new INILAB subroutine.
C  15 Apr 99 @DynTel-HQ 2-698 Add label for speed correction.
C  12 Apr 99 @DynTel-ddj 2-000  Added labels to allow bypassing of CO
C            offset correction calculation, RPV correction, temperature
C            correction and soak time adjustment to temperature
C            correction factors.
C  08 Apr 99 @DynTel-bag 2-000  Removed label for call to GETAMB (ambient
C            temperature)
C  24 Mar 99 @DynTel-DDJ 2-696  Added labels for CO phase 2 data
C  10 Mar 99 @DynTel-bag 2-664  Added label for weekend index; this will
C            allow user to specify that MOBILE6 calculations use weekend
C            activity distributions/emissions data. Also, added weekend
C            and weekday labels for running loss travel fraction input.
C            Added a label to enable output of hourly emissions.
C  11 Feb 99 @DynTel-MLA 2-663  Added label for start emissions in grams.
C   5 Jan 99 @DynTel-MLA 2-663  Added new labels for Mobile6 I/M programs.
C  30 Oct 98 @DynTel-BAG 2-677  Added label HOURLY TEMPERATURES for hourly
C            temperature determinations.
C  29 Sep 98 @DynTel-HQ 2-678 Conversion Factor
C  27 Aug 98 @DynTel-ZK 2-000 Explicit Typing
C  27 Jul 98 @DynTel-NH 0-632
C
C  Notes:
C
C     The derived type LABELDEF defines one variant form of an input
C     label. The structure contains the full and abbreviated forms of
C     the label text. It also contains a logical field, NEWLAB, which
C     is set to .TRUE. if there are no variant forms of the label. If
C     there are variant forms of a label, they should all be grouped
C     together in the LABELDEF list. (Values are assigned to LABELDEF
C     in the INILAB subroutine below.) The NEWLAB field should be set to
C     .TRUE. in the first record in the group and to .FALSE. in all the
C     remaining records in the group. The LABID field is set by the code.
C     (The INILAB subroutine uses the NEWLAB field to assign unique ID
C     numbers to different labels, but all of the variants of the same
C     label will be assigned the same ID number. The initial value of
C     the LABID field is not important.) The SECTION field should be set
C     to 0 if the label can appear in the header portion of an input file
C     and to some number greater than 0 if it can appear in the Run Data
C     section of the input file. The field should be set to 1 if the label
C     can appear only in the global data section of a test run, to 2 if 
C     the label can appear in either the global data or the scenario data
C     sections of a test run, or to 3 if the label can appear only in the
C     scenario section.
C
C     New labels can be added anywhere in the list. The only requirement
C     is that all variants of a single label must be grouped together in
C     the list, with the NEWLAB flag set to .TRUE. in the first variant
C     form to appear in the list and to .FALSE. in all the remaining forms.
C
C  Module variables:
C
C    LABLIST - Array of derived type LABELDEF (see module LABELS).
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  I        I    Loop control variable
C  LNUM     I    Unique token number for each label type.
C
      USE LABELS
      IMPLICIT NONE
C
      INTEGER   ::   I
      INTEGER   ::   LNUM
C
C. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C
C     Initialize the LABLIST array.
C
      LABLIST(1:7) = ( /
     &  LABELDEF('MOBILE6 INPUT FILE',
     &           'MOBILE6 INPUT FILE',.TRUE.,0,0),
     &  LABELDEF('MOBILE6 BATCH FILE',
     &           'MOBILE6 BATCH FILE',.TRUE.,0,0),
     &  LABELDEF('REPORT FILE',
     &           'REP FILE',.TRUE.,0,0),
     &  LABELDEF('RPT FILE',
     &           'RPT FILE',.FALSE.,0,0),
     &  LABELDEF('DATABASE OUTPUT',
     &           'DATABASE OUT',.TRUE.,0,0),
     &  LABELDEF('DBASE OUTPUT',
     &           'DB OUT',.FALSE.,0,0),
     &  LABELDEF('DATABASE OPTIONS',
     &           'DATABASE OPT',.TRUE.,0,0) /)
C
      LABLIST(8:16) = ( /
     &  LABELDEF('DBASE OPTIONS',
     &           'DB OPT',.FALSE.,0,0),
     &  LABELDEF('DATABASE EMISSIONS',
     &           'DATABASE EMIS',.TRUE.,0,0),
     &  LABELDEF('DBASE EMISSIONS',
     &           'DB EMIS',.FALSE.,0,0),
     &  LABELDEF('DATABASE FACILITIES',
     &           'DATABASE FAC',.TRUE.,0,0),
     &  LABELDEF('DBASE FACILITIES',
     &           'DB FAC',.FALSE.,0,0),
     &  LABELDEF('POLLUTANTS',
     &           'POL',.TRUE.,0,0),
     &  LABELDEF('DATABASE VEHICLES',
     &           'DATABASE VEH',.TRUE.,0,0),
     &  LABELDEF('DBASE VEHICLES',
     &           'DB VEH',.FALSE.,0,0),
     &  LABELDEF('EMISSIONS TABLE',
     &           'EMIS TAB',.TRUE.,0,0) /)
C
      LABLIST(17:25) = ( /
     &  LABELDEF('INPUTS TABLE',
     &           'INP TAB',.TRUE.,0,0),
     &  LABELDEF('DATAFILES TABLE',
     &           'DATA TAB',.TRUE.,0,0),
     &  LABELDEF('MESSAGES TABLE',
     &           'MESS TAB',.TRUE.,0,0),
     &  LABELDEF('RUNDATA TABLE',
     &           'RUND TAB',.TRUE.,0,0),
     &  LABELDEF('SCENARIOS TABLE',
     &           'SCEN TAB',.TRUE.,0,0),
     &  LABELDEF('ETYPES TABLE',
     &           'ETYP TAB',.TRUE.,0,0),
     &  LABELDEF('RUN DATA',
     &           'RUN DATA',.TRUE.,0,0),
     &  LABELDEF('TAMPERING RATES',
     &           'TAM RATE',.TRUE.,0,1),
     &  LABELDEF('NO TAMPERING',
     &           'NO TAMP',.TRUE.,0,1) /)
C
      LABLIST(26:34) = ( /
     &  LABELDEF('REGISTRATION DISTRIBUTIONS',
     &           'REG DIST',.TRUE.,0,1),
     &  LABELDEF('VEHICLE REGISTRATIONS',
     &           'VEH REG',.FALSE.,0,1),
     &  LABELDEF('MILEAGE ACCUMULATION RATES',
     &           'MILE ACCUM RATE',.TRUE.,0,1),
     &  LABELDEF('BASIC EFS',
     &           'BASIC EFS',.TRUE.,0,1),
     &  LABELDEF('SPEED COEFFICIENT',
     &           'SPEED COEF',.TRUE.,0,1),
     &  LABELDEF('SPD COEFFICIENT',
     &           'SPD COEF',.FALSE.,0,1),
     &  LABELDEF('SCF BY FWAY AND AC',
     &           'SCF BY FWAY AND AC',.TRUE.,0,1),
     &  LABELDEF('HD CONVERSION FACTORS',
     &           'HD CONV FACT',.TRUE.,0,1),
     &  LABELDEF('DISABLE 1990 CAAA',
     &           'DIS 1990 CAAA',.TRUE.,0,1) /)
C
      LABLIST(35:43) = ( /
     &  LABELDEF('I/M PROGRAM',
     &           'I/M PROG',.TRUE.,0,1),
     &  LABELDEF('I/M PGM',
     &           'I/M PGM',.FALSE.,0,1),
     &  LABELDEF('I/M MODEL YEARS',
     &           'I/M MOD YEAR',.TRUE.,0,1),
     &  LABELDEF('I/M EFFECTIVENESS',
     &           'I/M EFF',.TRUE.,0,1),
     &  LABELDEF('I/M CREDIT FILE',
     &           'I/M CRED FIL',.TRUE.,0,1),
     &  LABELDEF('I/M VEHICLES',
     &           'I/M VEH',.TRUE.,0,1),
     &  LABELDEF('I/M STRINGENCY',
     &           'I/M STRIN',.TRUE.,0,1),
     &  LABELDEF('I/M COMPLIANCE',
     &           'I/M COMP',.TRUE.,0,1),
     &  LABELDEF('I/M WAIVER RATES',
     &           'I/M WAIV RAT',.TRUE.,0,1) /)
C
      LABLIST(44:52) = ( /
     &  LABELDEF('I/M WAIVER REDUCTION',
     &           'I/M WAIV REDUC',.TRUE.,0,1),
     &  LABELDEF('I/M CUTPOINTS',
     &           'I/M CUT',.TRUE.,0,1),
     &  LABELDEF('I/M EXEMPTION AGE',
     &           'I/M EXEMP AGE',.TRUE.,0,1),
     &  LABELDEF('I/M GRACE PERIOD',
     &           'I/M GRA PER',.TRUE.,0,1),
     &  LABELDEF('I/M COIM PERCENTAGE',
     &           'I/M COIM PERCENT',.TRUE.,0,1),
     &  LABELDEF('I/M COIM PCT',
     &           'I/M COIM PCT',.FALSE.,0,1),
     &  LABELDEF('I/M INCLUDE COIM',
     &           'I/M INC COIM',.TRUE.,0,1),
     &  LABELDEF('NO I/M TTC CREDITS',
     &           'NO I/M TTC CRED',.TRUE.,0,1),
     &  LABELDEF('I/M TTC EFFECTIVENESS',
     &           'I/M TTC EFF',.TRUE.,0,1) /)
C
      LABLIST(53:61) = ( /
     &  LABELDEF('RSD FAILURE LIMIT',
     &           'RSD FAIL LIM',.TRUE.,0,1),
     &  LABELDEF('RSD TESTS PERFORMED',
     &           'RSD TEST PERF',.TRUE.,0,1),
     &  LABELDEF('RSD TOTAL VEHICLES',
     &           'RSD TOT VEH',.TRUE.,0,1),
     &  LABELDEF('RSD COVERAGE',
     &           'RSD COV',.TRUE.,0,1),
     &  LABELDEF('RSD NUMBER FAILURES',
     &           'RSD NUM FAIL',.TRUE.,0,1),
     &  LABELDEF('RSD NUMBER VEHICLES',
     &           'RSD NUM VEH',.TRUE.,0,1),
     &  LABELDEF('RSD CO CUTPOINT',
     &           'RSD CO CUT',.TRUE.,0,1),
     &  LABELDEF('RSD NOX CUTPOINT',
     &           'RSD NOX CUT',.TRUE.,0,1),
     &  LABELDEF('RSD LOSSES',
     &           'RSD LOSS',.TRUE.,0,1) /)
C
      LABLIST(62:70) = ( /
     &  LABELDEF('I/M RSD ONLY',
     &           'I/M RSD ONLY',.TRUE.,0,1),
     &  LABELDEF('I/M COIM ONLY',
     &           'I/M COIM ONLY',.TRUE.,0,1),
     &  LABELDEF('I/M DESCRIPTION FILE',
     &           'I/M DESC FILE',.TRUE.,0,1),
     &  LABELDEF('ANTI-TAMPERING PROGRAM',
     &           'ANTI-TAMP PROG',.TRUE.,0,1),
     &  LABELDEF('ANTI-TAMPERING PGM',
     &           'ANTI-TAMP PGM',.FALSE.,0,1),
     &  LABELDEF('GET ATP',
     &           'GET ATP',.FALSE.,0,1),
     &  LABELDEF('PRESSURE CHECK PARAMETERS',
     &           'PRES CHECK PARAM',.TRUE.,0,1),
     &  LABELDEF('PRESSURE CHK PARAMETERS',
     &           'PRES CHK PARAM',.FALSE.,0,1),
     &  LABELDEF('PURGE CHECK PARAMETERS',
     &           'PURGE CHECK PARAM',.TRUE.,0,1) /)
C
      LABLIST(71:79) = ( /
     &  LABELDEF('PURGE CHK PARAMETERS',
     &           'PURGE CHK PARAM',.FALSE.,0,1),
     &  LABELDEF('STAGE II REFUELING',
     &           'STAGE II REFU',.TRUE.,0,1),
     &  LABELDEF('STG II REFUELING',
     &           'STG II REFU',.FALSE.,0,1),
     &  LABELDEF('STAGE 2 REFUELING',
     &           'STAGE 2 REFU',.FALSE.,0,1),
     &  LABELDEF('STG 2 REFUELING',
     &           'STG 2 REFU',.FALSE.,0,1),
     &  LABELDEF('EXPAND LDT EFS',
     &           'EXP LDT EFS',.TRUE.,0,1),
     &  LABELDEF('EXPAND HDGV EFS',
     &           'EXP HDGV EFS',.TRUE.,0,1),
     &  LABELDEF('EXPAND HDDV EFS',
     &           'EXP HDDV EFS',.TRUE.,0,1),
     &  LABELDEF('EXPAND BUS EFS',
     &           'EXP BUS EFS',.TRUE.,0,1) /)
C
      LABLIST(80:87) = ( /
     &  LABELDEF('EXPAND EXHAUST',
     &           'EXP EXH',.TRUE.,0,1),
     &  LABELDEF('EXPRESS HC AS THC',
     &           'EXPR HC AS THC',.TRUE.,0,1),
     &  LABELDEF('EXPRESS HC AS NMHC',
     &           'EXPR HC AS NMHC',.TRUE.,0,1),
     &  LABELDEF('EXPRESS HC AS VOC',
     &           'EXPR HC AS VOC',.TRUE.,0,1),
     &  LABELDEF('EXPRESS HC AS TOG',
     &           'EXPR HC AS TOG',.TRUE.,0,1),
     &  LABELDEF('EXPRESS HC AS NMOG',
     &           'EXPR HC AS NMOG',.TRUE.,0,1),
     &  LABELDEF('EXPAND EVAPORATIVE',
     &           'EXP EVAP',.TRUE.,0,1),
     &  LABELDEF('NO EVAPORATIVE',
     &           'NO EVAP',.TRUE.,0,0) /)
C
      LABLIST(88:92) = ( /
     &  LABELDEF('NO EVP SAWTOOTH',
     &           'NO EVP SAWTOOTH',.TRUE.,0,3),
     &  LABELDEF('FIXED PP RATES',
     &           'FIX PP RATES',.TRUE.,0,3),
     &  LABELDEF('HOT SOAK ACTIVITY',
     &           'H S ACT',.TRUE.,0,3),
     &  LABELDEF('DIURN SOAK ACTIVITY',
     &           'D S ACT',.TRUE.,0,3),
     &  LABELDEF('NO REFUELING',
     &           'NO REFUEL',.TRUE.,0,1)/)
C
      LABLIST(93:100) = ( /
     &  LABELDEF('94+ LDG IMPLEMENTATION',
     &           '94+ LDG IMP',.TRUE.,0,1),
     &  LABELDEF('SCENARIO RECORD',
     &           'SCEN REC',.TRUE.,0,2),
     &  LABELDEF('WEEK DAY TRIP LENGTH DISTRIBUTIONS',
     &           'WE DA TRI LEN DI',.TRUE.,0,2),
     &  LABELDEF('WEEK END TRIP LENGTH DISTRIBUTIONS',
     &           'WE EN TRI LEN DI',.TRUE.,0,2),
     &  LABELDEF('VMT FRACTIONS',
     &           'VMT FRAC',.TRUE.,0,2),
     &  LABELDEF('VMT MIX',
     &           'VMT MIX',.FALSE.,0,2),
     &  LABELDEF('VMT BY FACILITY',
     &           'VMT BY FAC',.TRUE.,0,2),
     &  LABELDEF('VMT BY HOUR',
     &           'VMT BY HOUR',.TRUE.,0,2) /)
C
      LABLIST(101:107) = ( /
     &  LABELDEF('SPEED VMT',
     &           'SPEED VMT',.TRUE.,0,2),
     &  LABELDEF ('AC FRACTION',
     &            'AC FRACTION', .TRUE. , 0 , 2),
     &  LABELDEF('DETERGENT GASOLINE',
     &           'DET GAS',.TRUE.,0,2),
     &  LABELDEF('NO DETERGENT GAS',
     &           'NO DET GAS',.TRUE.,0,2),
     &  LABELDEF('OXYGENATED FUELS',
     &           'OXY FUEL',.TRUE.,0,2),
     &  LABELDEF('DIESEL FRACTIONS',
     &           'DIESEL FRAC',.TRUE.,0,2),
     &  LABELDEF('DIESEL SALES',
     &           'DIESEL SALES',.FALSE.,0,2) /)
C
      LABLIST(108:116) = ( /
     &  LABELDEF('NO LIQUID LEAKS',
     &           'NO LIQUID LEAKS',.TRUE.,0,3),
     &  LABELDEF('MIN/MAX TEMPERATURES',
     &           'MIN/MAX TEMP',.TRUE.,0,2),
     &  LABELDEF('HOURLY TEMPERATURES',
     &           'HOU TEMP',.TRUE.,0,2),
     &  LABELDEF('FUEL RVP',
     &           'FUEL RVP',.TRUE.,0,2),
     &  LABELDEF('FUEL PROGRAM',
     &           'FUEL PROGRAM',.TRUE.,0,2),
     &  LABELDEF('SEASON',
     &           'SEASON',.TRUE.,0,2),
     &  LABELDEF('ALTITUDE',
     &           'ALTITUDE',.TRUE.,0,3),
     &  LABELDEF('CALENDAR YEAR',
     &           'CAL YEAR',.TRUE.,0,3),
     &  LABELDEF('EVALUATION MONTH',
     &           'EVAL MON',.TRUE.,0,3) /)
C
      LABLIST(117:122) = ( /
     &  LABELDEF('NO RESTING LOSS',
     &           'NO REST LOSS',.TRUE.,0,2),
     &  LABELDEF('WEEKEND VEHICLE USAGE',
     &           'WE VEH US',.TRUE.,0,3),
     &  LABELDEF('T1-LEV EFS',
     &           'T1-LEV EFS',.TRUE.,0,3),
     &  LABELDEF('NO COLDCO STD',
     &           'NO COLDCO STD',.TRUE.,0,1),
     &  LABELDEF('START DIST',
     &           'START DIST',.TRUE.,0,1),
     &  LABELDEF('NO EVAP CALC',
     &           'NO EVAP CALC',.TRUE.,0,1)/ )
C
      LABLIST(123:130) = ( /
     &  LABELDEF('NO WEATHERED RVP',
     &           'NO W RVP',.TRUE.,0,1),
     &  LABELDEF('NO RVP CORRECTION',
     &           'NO RVP CORR',.TRUE.,0,3),
     &  LABELDEF('NO SOAK CORRECTION',
     &           'NO SOAK CORR',.TRUE.,0,3),
     &  LABELDEF('NO TEMPERATURE CF',
     &           'NO TEMP CF',.TRUE.,0,3),
     &  LABELDEF('NO TCF',
     &           'NO TCF',.FALSE.,0,3),
     &  LABELDEF('NO TEMP CORRECTION',
     &           'NO TEMP CORR',.FALSE.,0,3),
     &  LABELDEF('NO SPEED CORRECTION',
     &           'NO SPEED CORR',.TRUE.,0,3),
     &  LABELDEF('NO HDCGPM',
     &           'NO HDCGPM',.TRUE.,0,1)/)
C
      LABLIST(131:138) = ( /
     &  LABELDEF('NO RNL CORRECTION',
     &           'NO RNL CORR',.FALSE.,0,3),
     &  LABELDEF('FTP MILES PER DAY',
     &           'FTP MILES PER DAY',.TRUE.,0,3),
     &  LABELDEF('STARTS PER DAY',
     &           'STARTS PER DAY',.TRUE.,0,1),
     &  LABELDEF('SOAK DISTRIBUTION',
     &           'SOAK DIST',.TRUE.,0,3),
     &  LABELDEF('OBD WORKING MIL',
     &           'OBD WORK MIL',.TRUE.,0,1),
     &  LABELDEF('OBD EXH RESPONSE',
     &           'OBD EXH RESPONSE',.TRUE.,0,1),
     &  LABELDEF('OBD EVAP RESPONSE',
     &           'OBD EVAP RESPONSE',.TRUE.,0,1),
     &  LABELDEF('OBD I/M COVERAGE',
     &           'OBD I/M COV',.TRUE.,0,1) /)
C
      LABLIST(139:146) = ( /
     &  LABELDEF('OBD REPAIRED EFS',
     &           'OBD REP EF',.TRUE.,0,1),
     &  LABELDEF('NO HIGH EMITTERS',
     &           'NO HIGH EMIT',.TRUE.,0,1),
     &  LABELDEF('1981-93 LDG EFS',
     &           '1981-93 LDG EFS',.TRUE.,0,1),
     &  LABELDEF ('CLOUD COVER',
     &            'CLOUD COVER', .TRUE. , 0 , 2),
     &  LABELDEF ('SUNRISE/SUNSET',
     &            'SUNRISE/SUNSET', .TRUE. , 0 ,2 ),
     &  LABELDEF ('PEAK SUN',
     &            'PEAK SUN', .TRUE. , 0 , 2),
     &  LABELDEF ('NO AC',
     &            'NO AC', .TRUE. , 0 , 2),
     &  LABELDEF ('ABSOLUTE HUMIDITY',
     &            'ABSOLUTE HUMIDITY', .TRUE. , 0 , 2) /)
C
      LABLIST(147:155) = ( /
     &  LABELDEF('NO NOX PULL AHEAD',
     &           'NO NOX PULL AHEAD',.TRUE.,0,1),
     &  LABELDEF('REBUILD EFFECTS',
     &           'REBUILD EFFECTS',.TRUE.,0,1),
     &  LABELDEF('NO TIER2',
     &           'NO TIER2',.TRUE.,0,1),
     &  LABELDEF('T2 EXH PHASE-IN',
     &           'T2 EXH PHASE-IN',.TRUE.,0,1),
     &  LABELDEF('T2 CERT',
     &           'T2 CERT',.TRUE.,0,1),
     &  LABELDEF('NO REBUILD',
     &           'NO REBUILD',.TRUE.,0,1),
     &  LABELDEF('NO DEFEAT DEVICE',
     &           'NO DEFEAT DEVICE',.TRUE.,0,1),
     &  LABELDEF('SULFUR CONTENT',
     &           'SULFUR CONTENT',.TRUE.,0,3),
     &  LABELDEF('SULFUR CORRECTION FACTOR',
     &           'SULFUR CORR FACT',.FALSE.,0,3) /)
C
      LABLIST(156:164) = ( /
     &  LABELDEF('GASOLINE SULFUR',
     &           'GASOLINE SULFUR',.FALSE.,0,3),
     &  LABELDEF('NO DESC OUTPUT',
     &           'NO DESC OUT',.TRUE.,0,0),
     &  LABELDEF('NGV FRACTION',
     &           'NGV FRACTION',.TRUE.,0,1),
     &  LABELDEF('NGV EF',
     &           'NGV EF',.TRUE.,0,1),
     &  LABELDEF('NO SFTP SPEED',
     &           'NO SFTP SPEED',.TRUE.,0,2),
     &  LABELDEF('NO SFTP AC',
     &           'NO SFTP AC',.TRUE.,0,2),
     &  LABELDEF('DAILY OUTPUT',
     &           'DAILY OUTPUT',.TRUE.,0,0),
     &  LABELDEF('WITH FIELDNAMES',
     &           'WITH FIELDNAMES',.TRUE.,0,0),
     &  LABELDEF('T2 EVAP PHASE-IN',
     &           'T2 EVAP PHASE-IN',.TRUE.,0,1)/)
C
      LABLIST(165:173) = ( /
     &  LABELDEF('EVAP I/M FRAC BEN',
     &           'EVAP I/M FRAC BEN',.TRUE.,0,1),
     &  LABELDEF('GLL REDUCTION',
     &           'GLL REDUCTION',.TRUE.,0,1),
     &  LABELDEF('AGGREGATED OUTPUT',
     &           'AGGREGATED OUTPUT',.TRUE.,0,0),
     &  LABELDEF('DATABASE AGES',
     &            'DATABASE AGES',.TRUE.,0,0),
     &  LABELDEF('DATABASE YEARS',
     &            'DATABASE YEARS',.TRUE.,0,0),
     &  LABELDEF('DATABASE HOURS',
     &            'DATABASE HOURS',.TRUE.,0,0),
     &  LABELDEF('NO CLEAN AIR ACT',
     &           'NO CLEAN AIR ACT',.TRUE.,0,1),
     &  LABELDEF('NO 2007 HDDV RULE',
     &           'NO 2007 HDDV RULE',.TRUE.,0,1),
     &  LABELDEF('PARTICULATES',
     &           'PARTICULATES',.TRUE.,0,0)/)
C
      LABLIST(174:182) = ( /
     &  LABELDEF('PARTICLE SIZE',
     &           'PARTICLE SIZE',.TRUE.,0,3),
     &  LABELDEF('DIESEL SULFUR',
     &           'DIESEL SULFUR',.TRUE.,0,3),
     &  LABELDEF('NEW PM BER',
     &           'NEW PM BER',.TRUE.,0,3),
     &  LABELDEF('AIR TOXICS',
     &           'AIR TOXICS',.TRUE.,0,0),
     &  LABELDEF('E200',
     &           'E200',.TRUE.,0,3),
     &  LABELDEF('E300',
     &           'E300',.TRUE.,0,3),
     &  LABELDEF('GAS AROMATIC%',
     &           'GAS AROMATIC%',.TRUE.,0,3),
     &  LABELDEF('GAS OLEFIN%',
     &           'GAS OLEFIN%',.TRUE.,0,3),
     &  LABELDEF('GAS BENZENE%',
     &           'GAS BENZENE%',.TRUE.,0,3)/)
C
      LABLIST(183:191) = ( /
     &  LABELDEF('T2 PM PHASE-IN',
     &           'T2 PM PHASE-IN',.TRUE.,0,1),
     &  LABELDEF('IDLE PM EMISSIONS',
     &           'IDLE PM EMISSIONS',.TRUE.,0,1),
     &  LABELDEF('1996+ EXHAUST IM',
     &           '1996+ EXHAUST IM',.TRUE.,0,1),
     &  LABELDEF('ADDITIONAL HAPS',
     &           'ADDITIONAL HAPS',.TRUE.,0,3),
     &  LABELDEF('PARTICULATE EF',
     &           'PARTICULATE EF',.TRUE.,0,3),
     &  LABELDEF('OXYGENATE',
     &           'OXYGENATE',.TRUE.,0,3),
     &  LABELDEF('RVP OXY WAIVER',
     &           'RVP OXY WAIVER',.TRUE.,0,3),
     &  LABELDEF('AVERAGE SPEED',
     &           'AVERAGE SPEED',.TRUE.,0,3),
     &  LABELDEF('SPREADSHEET',
     &           'SPREADSHEET',.TRUE.,0,0)/)
C
      LABLIST(192:198) = ( /
     &  LABELDEF('MPG ESTIMATES',
     &           'MPG ESTIMATES',.TRUE.,0,3),
     &  LABELDEF('ROADFE',
     &           'ROADFE',.FALSE.,0,3),
     &  LABELDEF('HOT SOAK LQD LEAK',
     &           'HOT SOAK LQD LEAK',.TRUE.,0,3),
     &  LABELDEF('RELATIVE HUMIDITY',
     &           'REL HUM',.TRUE.,0,3),
     &  LABELDEF('BAROMETRIC PRES',
     &           'BARO PRES',.TRUE.,0,3),
     &  LABELDEF('END OF RUN',
     &           'END OF RUN',.TRUE.,0,3),
     &  LABELDEF('0000000000',
     &           '000',.FALSE.,0,3)/)
C
C
C     Loop through the label list. Set the IDNUM field, assigning a unique
C     value to each label that is not a variant form of a label that was
C     entered earlier.
C
      LNUM = 0
C
   10 DO 20 I = 1,NLABS
C
        IF (LABLIST(I)%NEWLAB) THEN
          LNUM = LNUM + 1
        ENDIF
C
        LABLIST(I)%LABID = LNUM
C
   20 CONTINUE
C
   90 RETURN
      END SUBROUTINE INILAB
