      SUBROUTINE QUITER(RVALUE,IVALUE,MESSAG,INERR) 
C 
C  QUITER prints error and warning messages on unit IOUERR. 
C 
C  Called by CALTEM, CHKPGM, CRANKC, EFCALX, EXMPTR, GET94P, 
C            GETAC,  GETATP, GETBYM, GETCOL, GETEFR, GETERQ, GETETP, 
C            GETFYR, GETHFL, GETIMC, GETIMP, GETIM2, GETNAT, 
C            GETPRG, GETPRS, GETSC1, GETTAM, GETVER, GETVMT, HCCALX, 
C            IBFPTR, IERPTR, IEVPTR, IGEPTR, IGLPTR, IGSFPT, IGTPTR, 
C            IGUPTR, IRNLPT, ISPPTR, ITAMPT, ITCPTR, LOCAL , MOBILE, 
C            REFUEL, RDCITY, STOERQ, TFCALX, WTLRNL, YRTEST. 
C 
C  Changes: (Last change first) 
C 
C  18 Jan 02 @EPA-djb Bug381: Changed text content of message 580, added 
C                     warnings 14, 15 & 16 to the list and changed their 
C                     device to IOUREP. Added warnings 514,515,581-586. 
C  15 Jan 02 AIR Task 15: Changed maximum alcohol oxy limit from 3.5% to 3.7% 
C                         (error message #6) 
C  14 Jan 02 @EPA-bag Bug377 , Changed "TO" to "T/O" in error message 82. 
C 14 Jan 02 @EPA-bag Bug375 , Added messages 577, 578, 579, &580 for the 
C                             AVERAGE SPEED command 
C  09 Jan 02 @EPA-djb Bug 356, No OBD I/M for pre-2007 HDGV (531). 
C  13 Dec 01 AIR Task 27-28: Fixed error message 83. 
C  04 Dec 01 AIR Task 27-28: Added CO2 error when NGV are selected. 
C  02 Oct 01 AIR Task 15-21: Extended the input AT Oxygenate from 5% to 20 vol% 
C  24 Aug 01 AIR Task 15-21: Added NAMEAT and aix toxic messages 115 and 116. 
C 23 May 01 ARI Task 17: Added message for PM speed correction. Expanded 
C           NAMEPART to 3 for the PARTICULATE EF (ZML/DR FILENAMES). 
C 30 Apr 01 AIR Task 12: Added messages for GETT2PMFRC. 
C 11 Apr 01 AIR Task 11: Changed messages that are true errors to GOTO 198 
C           and output via IOUERR file, while warnings/comments GOTO 199 
C           and output via IOUREP/IOUOUT.  
C 07 Feb 01 AIR Task 11: Added messages for Particulates and Air Toxics. 
C  19 Jul 01 @EPA-djb Wish235, Changed unused warnings 14, 15 and 16 
C            to average speed command warning messages. 
C  10 Jul 01 @EPA-djb Wish130, Changed errors 183 and 186 for use in GETSPD. 
C  25 May 01 @EPA-ELG message 567 for post 1996 I&M 
C  21 Feb 01 @EPA-bag QUITER message 174 for user supplied faclity mileage fractions 
C  19 Jan 01 @EPA-bag QUITER message for NO 2007 HEAVY DUTY DIESEL VEHICLE RULE; 
C            removed obsolete error message #83; error message #53 used for fatal 
C            divide by zero error in RSTLOS.FOR 
C  28 Dec 00 @EPA-bag bug294 QUITER message for IREJN out of bounds is now 
C            message #24 and changed message text; changed text in error #'s 153, 
C            35, 36, 82, 162; changed i&m model year out of bounds errors to 
C            QUITER messages 12 and 13 
C 03 Jan 00 @EPA-elg task-x6, added messages 565 and 566 for the NO CLEAN AIR ACT command 
C           (see CKOBDIM.FOR) 
C 21 Nov 00 @EPA-bag Bug303, removed redundant and confusing message #'s 606 and 607 
C           which concerned running loss trip fractions 
C 18 Oct 00 @EPA-bag Bug277, removed error message #111 
C  25 Oct 00 @EPA-djb Bug292, Changed text of ONE-TIME error message; ATP input error 
C           message #107 
C 11 Oct 00 AIR Task 02: Bug Fix 275: Added read error message #'s 562,563,564 
C 19 Oct 00 @EPA-bag  Removed obsolete error message #525 
C 02 Oct 00 AIR Task 07: Bug Fix 265, removed reference to M5b User Guide in 
C           message 112. 
C 06 Oct 00 @EPA-bag  BUG137&155 removed error message #12-16(cold CO; GETCOL), 531(GETFTP), 
C           115,116,140-144(GETFTP), 
C 02 Oct 00 @EPA-bag  BUG183 removed error message #524 
C 26 Sep 00 AIR Task 08: Added warning concerning evap I/M and NO EVP SAWTOOTH(# 520). 
C 21 Sep 00 AIR Task 07: Bug 192 Added messages for user-supplied AC input(#'s 618-620). 
C 22 Sep 00 AIR Task 07: Bug Fix 183 Trap incorrect number of flags entered 
C           for GETIMV (I&M vehicle selection, #129). 
C 25 Sep 00 @EPA-elg  BUG 257 Modified error code 68. 
C 27 Sep 00 @EPA-bag bug264&258 removed unused message #'s 129,174,181,182,184, 
C           512,514,515,520 
C 28 Aug 00 AIR Task 08: Added warning messages for several EVAP I/M situations. 
C 15 Sep 00 @EPA-elg  Added call to quiter for T2_EVAP diagnostic. 
C  11 Sep 00 @EPA-elg Added min/max check on sulfur level input. 
C 20 Aug 00 @EPA-bag  Added warning message for disabling of evaporative 
C           emissions calculations 
C 03 Aug 00 @EPA-elg   Minimum User suppplied I/M compliance rate is 
C                      now 50% changed from 0.%(message #34) 
C  11 Aug 00 @EPA-bag Task 06: Removed WARNMSG subroutine redundancy with 
C            calls to QUITER for evaporative diagnostic labels; error messages 
C            696-699 used for these EPA diagnostic warning messages. Added more 
C            descriptive text to diagnostic warning message at line labeled '195'. 
C 12 Jul 00 AIR Task 06: Added messages for Hot Soak and Diurnal Soak Activity 
C           and Fixed P/P Rates. 
C    4 Aug 00 @EPA-bag Bug62  Revised warning 154 for 1981-1993 LDG 
C             entered with label BASIC EFS. 
C   21 Jul 00 @EPA-BAG running loss  Added message 617 for resting losses disabled 
C   10 Jul 00 @EPA-ddj Bug169  Added message 616 for user supplied sulfur levels 
C    3 Jul 00 @EPA-ddj Bug169  Added messages 614 and 615 for user supplied 
C             diesel sales fractions and vehicle miles, VMT's 
C             entered with label BASIC EFS. 
C    3 Jul 00 @EPA-ddj Bug62  Added warning 154 for 1981-1993 LDG 
C             entered with label BASIC EFS. 
C  5 Jul 00 @EPA-ddj Bug46 Setup messages 601 through 699 to handle 
C           comments. 
C 31 May 00 @EPA-ddj Bug102, Removed QUITER warning message 170 
C           (detergent additive regulations). 
C 22 May 00 @EPA-elg  Added Exhaust I/M and OBD I/M Error Labels. 
C 01 May 00 @EPA- BG  Removed use of DATA statement initializations 
C 17 Apr 00 @EPA-elg Added NGV Error Labels. 
C 11 Apr 00 AIR Task 02: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C           so that Lahey Standard for keyboard/monitor I/O is followed. 
C 05 Apr 00 AIR Task 02: Removed messages 149, 154 and 170, but left 
C           code in so that the location can be reused. 
C 07 Apr 00 @EPA-djb Fixed NAMVEH labels. 
C 10 mar 00 @DynTel-MLA 1-011  Added messages 531 for database output 
C 13 Jan 00 @DynTel-ddj 1-013  Added messages to alert user 
C           to changes in values for AC parameters. 
C 03 Jan 00 @DynTel-ddj 1-021  Added msgs. 527-530 for tier2 phase in. 
C 16 Dec 99 @DynTel-ddj 1-010  Added message 526.  Message alerts 
C           users that the diablement of CAAA is not active with 
C           current version of MOBILE6. 
C 24 Aug 99 @DynTel-MLA 2-684  Added messages 524 and 525. 
C 13 Jul 99 @dynTel-ddj 2-699  Added QUITER(messages 11 and 177 for 
C           sulfur input file. 
C  7 Jul 99 @DynTel-MLA 2-694  Changed wording of messsage 175. 
C 22 Jun 99 @DynTel-MLA 2-694  Deleted messages 11 (LEVYRS out of bounds) 
C           and 177 (items in P94REF/P94SEF must be positive) 
C 21 Jun 99 @DynTel-bg  2-705  Added a warning for user input of OBD I/M 
C         vehicle models less than 1994. (message #523) 
C  2 Jun 99 @DynTel-MLA 2-694  Added message 522 for OBD I/M. 
C 17 May 99 @DynTel-bg 2-701 Replaced + and - character output with spaces 
C 13 Apr 99 @DynTel-bg 2-000 Changed message for error #10 
C 25 Mar 99 @DynTel-DDJ 2-696 Changed message for error #15 
C 16 Mar 99 @DynTel-BG 2-664 Removed message #97 concerning Wade index; 
C           CALUDI has been deleted from MOBILE6. Replaced #97 with 
C           with message for reading of user supplied running loss 
C           trip length fractions Removed message #146 for deleted diurnal 
C           code. Added output of read errors(#10 and #146) for hourly 
C           filename processing. 
C  4 Mar 99 @DynTel-MLA 2-659 Added GET94P to 'Called by' line. 
C 25 Feb 1999 @Dyntel-MLA 2-663  Added message 171 for new I/M code 
C 21 Dec 1998 @DynTel-ZK 2-000 Explicit Typing 
C 19 Nov 1998 @DynTel-MLA 2-633  Removed message 163 relating to RBH 
C             programs. Rewrote message for error 82 and made it non-fatal. 
C             Also made 162 non-fatal. 
C  9 Nov 1998 @Dyntel-BAG 2-000  Deleted message 10. 
C  9 Nov. 98 @DynTel-HXQ 2-680 to change CML BUS to SCH BUS. 
C 14 Sep 1998 @Dyntel-MLA 0-632  Deleted old message 129. Added new 
C             message 129 and mesages 153, 172, 178-190, and 
C             501-517. Expanded code to use message numbers greater 
C             than 500. 
C 14 Sep 1998 @Dyntel-RJD 0-632  Added message 152. 
C 
C  Input on call: 
C 
C    argument list: RVALUE,IVALUE,MESSAG,INERR 
C    common blocks: 
C    /ATPAR2/ ATPPGM,ATPFQT,CRATP 
C    /BASEQ4/ NEWPAR,MAXCT 
C    /BASEQ6/ MAXERU 
C    /FLAGS2/ IMFLAG 
C    /FLAGS3/ ATPFLG 
C    /IMPAR1/ CRIM 
C    /IMPAR6/ IFREQ,INTYP 
C    /IOUCOM/ IOUREP,IOUERR,IOUOUT 
C    /NAMES1/ NAMFLG,NAMRVP,NAMTEM 
C    /PRGCH2/ PRGPGM,PRGFQT,CRPRG 
C    /PRSCH2/ PRSPGM,PRSFQT,CRPRS 
C    /QUITXQ/ N1QUIT 
C    /YEARS4/ IYSTART,IY1960,IYEND 
C 
C  Output on return: 
C 
C    argument list: INERR 
C 
C  Warning: there is a conditional branch to STOP at statement number 198. 
C 
C  Local array subscripts: 
C 
C  IOOFBD(5)    -  IOOFBD ( ICH ) 
C  NAMEFF(2)    -  NAMEFF ( IVEFF ) 
C  NAMFNC(17)   -  NAMFNC ( ICH ) 
C  NAMFRQ(3)    -  NAMFRQ ( IMFREQ(1/2)/(ATP/PRS/PRG)FQT ) 
C  NAMPOL(3)    -  NAMPOL ( IP ) 
C  NAMTYP(2,4)  -  NAMTYP ( INTYP(1/2)/(ATP/PRS/PRG)PGM ) 
C  NAMVEH(8)    -  NAMVEH ( IV ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  IOOFBD   C*4  character string vector: 'out of bounds for' 
C  MAXMES   I    maximum message code value 
C  MCODE    I    = MESSAG = message code, unless < 0 or > MAXMES 
C  MCOM     C*7  message type = "Comment" prefix 
C  MERR     C*7  message type = "Error" prefix 
C  MTYPE    C*7  message type prefix, determined by MCODE 
C  MWARN    C*7  message type = "Warning" prefix 
C  NAMEFF   I    character key for Stage II vrs % efficiency IV label 
C  NAMDET   C*10 detergent gasoline names 
C  NAMEIM   C*19 I&M required input names  
C  NAMFNC   C*8  function names (to identify subprogram source of error) 
C  NAMFRQ   C*8  program inspection frequency cases 
C  NAMPOL   C*5  pollutant names 
C  NAMTYP   C*8  program inspection type cases 
C  NAMVEH   C*20 vehicle type names 
C  OBDTYPE  C*16 OBD input names 
C  QFREQ     I   index for NAMFRQ 
C  QINTYP    I   index for NAMTYP 
C  QAFREQ    I   index for NAMFRQ 
C  QFQTFN    I   index for NAMFRQ 
C  QPGMFN    I   index for NAMTYP 
C  QAPROG    I   index for NAMTYP 
C  MESSAG    I   integer message code 
C  INERR     I   error count 
C  IVALUE    I   value determined to be out of bounds for a particular variable  
C  I         I   used for indexed write statement 
C  IPAR      I   used for indexed write statement of idle coefficients 
C  VEHICLEGROUP C*19 vehicle grouping names for I&M programs 
C 
C  Notes: 
C 
C  Some specific things about QUITER are: 
C  All messages must fit in 80 columns 
C  A blank needs to start at the beginning of a line so it does 
C  not get cut off in MTS 
C  Each line of a warning is preceded by 15X and errors 13X 
C  Use available message slots before expanding MAXMES upper limit 
C  Messages should be placed with other messages that pertain to the same 
C  topic when possible.  A group of new messages that pertain to the same 
C  topic should be placed together rather than scattering them in the  
C  available message slots. 
C  The capability to print the offending value is available 
C  After a message prompt goto 199 if the message is not 
C  included in the final count of INERR's otherwise goto 198 
C 
C  June-14-1993 @ ARC-bsg Subtask 222 message 74 
C  July-13-1993 @ ARC-bsg Subtask 277 message 69 and 70 removed 
C  idle calculation has changed to no longer require these messages 
C 
C  TMM@ARC, Request 2-366, 1/4/94 
C    added message 112 for reformulated fuel. 
C 
C  April-4-1994 @ CSC-tmm, Request 2-411 
C    added messages 154, 155, 156 for on-board vapor recovery regulations 
C 
C  7-May-1994 @ CSC-tmm Request #2-441 Modified code to permit the 
C     detergent gasoline parameters to be read in and stored in COMMON. 
C     Added msgs. 69,70, & 85. 
C  Aug-18-1994 @ CSC-bsg Request 2-462 Added QUITER messages 157-160, for  
C  non-FTP related messages. 
C  30-August-1994 @CSC-pme, request 2-455. Prompt error message was updated. 
C  Sep-14-1994 @ CSC-bsg Request 2-467 New Non-FTP Effects message and the 
C  elimination of error message 157 
C 
C  May-1995 @ CSC-zp Req. 2-495 "Add MOB5a-H to MOBILE5b" 
C  New message 161 if error occurred in  reading PPEFF.D; 
C  New message 162 if one Pressure or Purge program is 3 and the other is not 
C  New message 163 if error occurred in reading TTC.D file  
C  New message 164 if error occurred in reading "I/M Control Flag Record" 
C  New message 165 if value out of bounds in the "TTC Record" 
C  New message 166 if I/M + TTC credit exceed 100%               
C 
C   9 Feb 96 @DynTel-MLA 2-611 Added messages 102, 113, 157, 167, and 168. 
C  13-Sep-1996 @DynTel-yc Request #2-630 Added messages 169 and 170 for 
C               detergent gasoline effects. 
C  25 Feb 97 DynTel@MLA 2-622  Added messages 171, 172, 173, and 174 for 
C            Mobile6 and expanded the vehicle type name array. 
C  10-Oct-1997 @DynTel-rwk 2-000 More Quiter guidance comments were added. 
C  15-Oct-1997 @DynTel-gjr 2-000 Changed "MOBILE5" to "MOBILE6" in error  
C              message 153           
C 
C  10-24-1997 @Dyntel-MMS 2-645 error message 152 was deleted 
C 
C  10-30-97   @Dyntel-gjr 2-645 error message 71 was deleted. 
C 
C   1/27/98   @Dyntel-MLA 2-645 Error message 172 was deleted. Message #1 
C             is now used instead (spreadsheet flag out-of-bounds). 
C 
C  10 Jun 98  @Dyntel-RJD 2-650 new error message 71 was added. 
C  
C  24 Jun 98  @Dyntel-bag 2-000 error messages 74, 153 and 174 were removed. 
C  
C  30 Jun 98  @Dyntel-bag 2-000 error message 171 was removed. IVPICK = 3 
C             is no longer an option. 
C  30 Jun 98  @Dyntel-MLA 2-659 error messages 175-177 were added. 
C   9 Jul 98  @Dyntel-MLA 2-659 Added message 74. 
C  23 Jul 98  @DynTel-MLA 2-000  Changed LDGT names to conform to 
C             certification naming conventions. 
C   4 Aug 98  @DynTel-MLA 2-000  Added message 178. 
C 
      IMPLICIT NONE 
      INCLUDE 'ATPAR2.I' 
      INCLUDE 'BASEQ4.I' 
      INCLUDE 'BASEQ6.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'IMPAR1.I' 
      INCLUDE 'IMPAR6.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'NAMES1.I' 
      INCLUDE 'PRGCH2.I' 
      INCLUDE 'PRSCH2.I' 
      INCLUDE 'YEARS4.I' 
C 
      INTEGER, INTENT(IN) :: MESSAG 
      INTEGER, INTENT(IN) :: IVALUE 
      INTEGER, INTENT(INOUT) :: INERR 
      REAL, INTENT(IN) :: RVALUE 
C 
      INTEGER QFREQ,QINTYP,QAFREQ,QFQTFN,QPGMFN, QAPROG, MCODE,I,IPAR 
      INTEGER IOS 
C 
      DIMENSION QFREQ(2),QINTYP(2),QFQTFN(2),QPGMFN(2) 
C 
      CHARACTER*1, DIMENSION(2) :: NAMEFF=(/'L','H'/) 
C 
      CHARACTER*4, DIMENSION(5) :: IOOFBD=(/ 
     * '  ou','t of',' bou','nds ','for '/) 
C 
      CHARACTER*8 :: MERR = 'Error:  ' 
      CHARACTER*8 :: MWARN = 'Warning:' 
      CHARACTER*8 :: MCOM = 'Comment:' 
      CHARACTER*8 :: MTYPE 
      CHARACTER*19, DIMENSION(3) :: VEHICLEGROUP =(/ 
     *   'LIGHT DUTY GAS    ', 
     *   'HEAVY DUTY GAS    ', 
     *   'BUS               '/) 
      CHARACTER*8, DIMENSION(28) :: NAMVEH =(/ 
     *   'LDGV    ', 
     *   'LDGT1   ', 
     *   'LDGT2   ', 
     *   'LDGT3   ', 
     *   'LDGT4   ', 
     *   'HDGV2b  ', 
     *   'HDGV3   ', 
     *   'HDGV4   ', 
     *   'HDGV5   ', 
     *   'HDGV6   ', 
     *   'HDGV7   ', 
     *   'HDGV8a  ', 
     *   'HDGV8b  ', 
     *   'LDDV    ', 
     *   'LDDT12  ', 
     *   'HDDV2b  ', 
     *   'HDDV3   ', 
     *   'HDDV4   ', 
     *   'HDDV5   ', 
     *   'HDDV6   ', 
     *   'HDDV7   ', 
     *   'HDDV8a  ', 
     *   'HDDV8b  ', 
     *   'MC      ', 
     *   'HDGB    ', 
     *   'HDDBT   ', 
     *   'HDDBS   ', 
     *   'LDDT34  '/) 
C 
      CHARACTER*8, DIMENSION(17) :: NAMFNC = (/ 
     *   ' IEVPTR ', 
     *   ' OLDTPD ', 
     *   'not used', 
     *   'not used', 
     *   ' CRANKC ', 
     *   ' IGEPTR ', 
     *   ' IGUPTR ', 
     *   ' IGLPTR ', 
     *   ' IGTPTR ', 
     *   ' IBFPTR ', 
     *   ' IERPTR ', 
     *   ' ISP C1 ', 
     *   ' ISP C2 ', 
     *   ' EXMPTR ', 
     *   ' ITCPTR ', 
     *   ' IGSFPT ', 
     *   ' IDLPTR '/) 
C 
      CHARACTER(LEN=9), DIMENSION(3) :: NAMFRQ=(/ 
     *             'Annual   ','Biennial ','Triennial'/) 
C 
      CHARACTER(LEN=5), DIMENSION(3)::NAMPOL=(/' HC  ',' CO  ',' NOx '/) 
C 
      CHARACTER(LEN=20), DIMENSION(2,4) :: NAMTYP=RESHAPE((/ 
     *            'Test Only           ','Test Only           ', 
     *            'Comp Test & Repair  ','Test & Repair       ', 
     *            'Manual Test & Repair','Test & Repair       ', 
     *            'Not Applicable      ','Not Applicable      '/), 
     *             (/2,4/)) 
C 
      CHARACTER(LEN=10), DIMENSION(9) :: NAMDET = (/ 
     *   'DGPHSF    ', 
     *   'FIDFRC    ', 
     *   'IVDFRC    ', 
     *   'FIDEFF(1) ', 
     *   'FIDEFF(2) ', 
     *   'FIDEFF(3) ', 
     *   'IVDEFF(1) ', 
     *   'IVDEFF(2) ', 
     *   'IVDEFF(3) '/) 
C 
      CHARACTER(LEN=19), DIMENSION(5) :: NAMEIM=(/ 
     *   'I/M MODEL YEARS ', 
     *   'I/M VEHICLES    ', 
     *   'I/M STRINGENCY  ', 
     *   'I/M COMPLIANCE  ', 
     *   'I/M WAIVER RATES'/) 
C 
      CHARACTER(LEN=16), DIMENSION(5) :: OBDTYPE=(/ 
     *   'OBD WORKING MIL ', 
     *   'OBD RESPONSE    ', 
     *   'OBD I/M COVERAGE', 
     *   'OBD REPAIRED EFS', 
     *   'OBD PCT EQUIPPED'/) 
C 
      CHARACTER(LEN=14), DIMENSION(3) :: NAMEPART=(/ 
     *   'PARTICLE SIZE ', 
     *   'DIESEL SULFUR ', 
     *   'PARTICULATE EF'/) 
C 
      CHARACTER(LEN=14), DIMENSION(6) :: NAMEAT=(/ 
     *   'E200          ', 
     *   'E300          ', 
     *   'GAS AROMATIC% ', 
     *   'GAS OLEFIN%   ', 
     *   'GAS BENZENE%  ', 
     *   'OXYGENATE     '/) 
C 
C  Aug-18-1994 @ CSC-bsg Request 2-462 Added QUITER messages 157-160, for  
C  non-FTP related error messages 
C 
C 
      INTEGER, PARAMETER :: MAXMES = 699 
C 
      MCODE=MESSAG 
C 
C  Trap unknown errors. 
C 
      IF(MCODE.LT.1.OR.(MCODE.GE.191.AND.MCODE.LE.500).OR. 
     &   MCODE.GT.MAXMES) MCODE=0 
C 
C  Identify type of message. 
C 
      MTYPE=MERR 
      IF(MCODE.EQ.10 .OR. MCODE.EQ.14 .OR. MCODE.EQ.15 
     *   .OR. MCODE.EQ.16 
     *   .OR. MCODE.EQ.21 .OR. MCODE.EQ.22 .OR. MCODE.EQ.31 
     *   .OR. MCODE.EQ.45 .OR. MCODE.EQ.47 
     *   .OR. MCODE.EQ.48 .OR. MCODE.EQ.49 
     *   .OR. MCODE.EQ.52 .OR. (MCODE.GE.56.AND.MCODE.LE.59) 
     *   .OR. MCODE.EQ.74 .OR. MCODE.EQ.90 
     *   .OR. MCODE.EQ.96 .OR. MCODE.EQ.98 
     *   .OR.(MCODE.GE.108.AND.MCODE.LE.112) 
     *   .OR. MCODE.EQ.114 .OR. (MCODE.GE.119.AND.MCODE.LE.121) 
     *   .OR. MCODE.EQ.135 .OR. MCODE.EQ.146 .OR. MCODE.EQ.152 
     *   .OR. (MCODE.GE.154.AND.MCODE.LE.156) 
     *   .OR. MCODE.EQ.158 .OR. MCODE.EQ.165 .OR. MCODE.EQ.166 
     *   .OR. MCODE.EQ.514 .OR. MCODE.EQ.515 
     *   .OR. MCODE.EQ.520 .OR. MCODE.EQ.521 
     *   .OR. (MCODE .GE. 581 .AND. MCODE.LE. 586) 
     *   .OR. (MCODE .GE. 691 .AND. MCODE.LE. 699)) 
     *   MTYPE=MWARN 
      IF(MCODE.EQ.44.OR.MCODE.EQ.83.OR.MCODE.EQ.517 
     *   .OR. (MCODE.GE.601 .AND. MCODE.LE.690) ) 
     *   MTYPE=MCOM 
C 
C  Write message type prefix: 
C 
      IF(MTYPE.EQ.MERR) THEN 
        WRITE(IOUERR,299) MCODE,MTYPE 
      ELSE 
        WRITE(IOUREP,299) MCODE,MTYPE 
      ENDIF 
  299 FORMAT(T3,'M',I3,1X,A8) 
C 
C  Branch to the appropriate error / warning / comment message. 
C 
      IF (MCODE.GT.500) GOTO 1000 
      GOTO( 
     *      1,   2,   3,   4,   5,   6,   7,   8,   9,  10, 
     *     11,  12,  13,  14,  15,  16,  17,  18,  19,  20, 
     *     21,  22,  23,  24,  25,  26,  27,  28,  29,  30, 
     *     31,  32,  33,  34,  35,  36,  37,  38,  39,  40, 
     *     41,  42,  43,  44,  45,  46,  47,  48,  49,  50, 
     *     51,  52,  53,  54,  55,  56,  57,  58,  59,  60, 
     *     61,  62,  63,  64,  65,  66,  67,  68,  69,  70, 
     *     71,  72,  73,  74,  75,  76,  77,  78,  79,  80, 
     *     81,  82,  83,  84,  85,  86,  87,  88,  89,  90, 
     *     91,  92,  93,  94,  95,  96,  97,  98,  99, 100, 
     *    101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 
     *    111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 
     *    121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 
     *    131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 
     *    141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 
     *    151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 
     *    161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 
     *    171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 
     *    181, 182, 183, 184, 185, 186, 187, 188, 189, 190) MCODE 
      GOTO 197 
C 
 1000 CONTINUE 
      GOTO( 
     *    501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 
     *    511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 
     *    521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 
     *    531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 
     *    541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 
     *    551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 
     *    561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 
     *    571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 
     *    581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 
     *    591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 
     *    601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 
     *    611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 
     *    621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 
     *    631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 
     *    641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 
     *    651, 652, 653, 654, 655, 656, 657, 658, 659, 660, 
     *    661, 662, 663, 664, 665, 666, 667, 668, 669, 670, 
     *    671, 672, 673, 674, 675, 676, 677, 678, 679, 680, 
     *    681, 682, 683, 684, 685, 686, 687, 688, 689, 690, 
     *    691, 692, 693, 694, 695, 696, 697, 698, 699) 
     &   (MCODE-500) 
      GOTO 197 
C 
C  Out of bounds errors in flags land here (flag # is in RVALUE) 
C 
    1 MCODE=IFIX(RVALUE) 
      WRITE(IOUERR,301) IVALUE,IOOFBD,NAMFLG(MCODE) 
  301 FORMAT(13X,I4,5A4,'flag',A8) 
      GOTO 198 
    2 WRITE(IOUERR,302) RVALUE,IOOFBD 
  302 FORMAT(13X,G10.3,5A4,'Ether Blend Market Share (0. to 1.)') 
      GOTO 198 
    3 WRITE(IOUERR,303) RVALUE,IOOFBD 
  303 FORMAT(13X,G10.3,5A4,'Alcohol Blend Market Share (0. to 1.)') 
      GOTO 198 
C 
    4 GOTO 199 
C 
    5 WRITE(IOUERR,305) RVALUE,IOOFBD 
  305 FORMAT(13X,G10.3,5A4,'Ether Blend Oxygen Weight Fraction', 
     *    ' (0 to 0.027)') 
      GOTO 198 
C 
    6 WRITE(IOUERR,306) RVALUE,IOOFBD 
  306 FORMAT(13X,G10.3,5A4,'Alcohol Blend Oxygen Weight Fraction', 
     *    ' (0 to 0.037)') 
      GOTO 198 
C 
    7 WRITE(IOUERR,307) IVALUE,IOOFBD 
  307 FORMAT(13X,I4,5A4,'Alcohol Blend RVP Waiver switch (1 or 2)') 
      GOTO 198 
C 
    8 WRITE(IOUERR,308) 
  308 FORMAT(13X,'Inconsistent Ether Blend input encountered') 
      GOTO 198 
C 
    9 WRITE(IOUERR,309) 
  309 FORMAT(13X,'Inconsistent Alcohol Blend input encountered') 
      GOTO 198 
C 
   10 WRITE(IOUREP,310) 
  310 FORMAT( 
     *15X,'Error occurred in reading the hourly output filename.' 
     */15X,'The filename specified for hourly output already exists.' 
     */15X,'No hourly output will be produced for this scenario.') 
      GOTO 199 
C 
   11 WRITE(IOUERR,311) RVALUE 
  311 FORMAT( 
     *'+',12X,' Error occurred in reading user supplied', 
     *        ' gasoline sulfur content.'/, 
     * 14x,'Sulfur content of ',F6.1,' is greater than upper allowed', 
     *' bound of 600 ppm.',/, 
     * 14X,'Sulfur content set equal to 600 ppm',/) 
      GOTO 198 
C 
   12 WRITE(IOUERR,312) IVALUE,IOOFBD,IY1960,IYEND 
  312 FORMAT(13X,I4,5A4,'MODYR1 (1941 to 2050)') 
      GOTO 198 
C 
   13 WRITE(IOUERR,313) IVALUE,IOOFBD,IY1960,IYEND 
  313 FORMAT(13X,I4,5A4,'MODYR2 (1941 to 2050)') 
      GOTO 198 
C 
   14 WRITE(IOUREP,314) RVALUE 
  314 FORMAT(13X,'The combined area wide average speed entered cannot ', 
     *  'be ',/,13X,'less than ',F4.1,' miles per hour.',/, 
     *  13X,'The average speed will be reset to this value.') 
      GOTO 199 
C 
   15 WRITE(IOUREP,315) RVALUE 
  315 FORMAT(13X,'The combined area wide average speed entered cannot ', 
     *  'be ',/,13X,'greater than ',F4.1,' miles per hour.',/, 
     *  13X,'The average speed will be reset to this value.') 
      GOTO 199 
C 
   16 WRITE(IOUREP,316) RVALUE 
  316 FORMAT( 
     * 13X,'The VMT distribution by facility (',F7.3,') provided',/, 
     * 13X,'with the AVERAGE SPEED command is not equal to one or',/, 
     * 13X,'to 100%. The user supplied values will be normalized.') 
      GOTO 199 
C 
   17 WRITE(IOUERR,317) RVALUE,IOOFBD 
  317 FORMAT(13X,G10.3,5A4,'VMTMIX (0. to 1.)') 
      GOTO 198 
   18 WRITE(IOUERR,318) RVALUE 
  318 FORMAT(13X,G10.3,' sum of VMTMIX is not equal to 1.') 
      GOTO 198 
   19 WRITE(IOUERR,319) RVALUE 
  319 FORMAT(13X,G10.3,' negative model year mileage') 
      GOTO 198 
   20 WRITE(IOUERR,320) RVALUE 
  320 FORMAT(13X,G10.3,' negative model year registration') 
      GOTO 198 
   21 WRITE(IOUREP,321) NAMVEH(IVALUE) 
  321 FORMAT(15X,' Mileage accumulation missing for one age of ',A8) 
      GOTO 199 
   22 WRITE(IOUREP,322) NAMVEH(IVALUE) 
  322 FORMAT(15X,' Age distribution is zero for class ',A8) 
      GOTO 199 
   23 WRITE(IOUERR,323) IVALUE,IOOFBD,MAXCT 
  323 FORMAT(13X,I4,5A4,'e.f. changes (1 to',I4,' pairs)') 
      GOTO 198 
   24 WRITE(IOUERR,324) IVALUE,IOOFBD 
  324 FORMAT(13X,I4,5A4,'altitude (1 or 2)') 
      GOTO 198 
   25 WRITE(IOUERR,325) IVALUE,IOOFBD 
  325 FORMAT(13X,I4,5A4,'vehicle type (1 to 28)') 
      GOTO 198 
   26 WRITE(IOUERR,326) IVALUE,IOOFBD 
  326 FORMAT(13X,I4,5A4,'pollutant (1 to 3)') 
      GOTO 198 
   27 WRITE(IOUERR,327) IVALUE,IOOFBD,IYSTART,IYEND 
  327 FORMAT(13X,I4,5A4,'year (',I4,' to ',I4,')') 
      GOTO 198 
   28 WRITE(IOUERR,328) 
      WRITE(IOUREP,328) 
  328 FORMAT('0',T2,'Excess data errors prevent further analysis.') 
      GOTO 199 
   29 WRITE(IOUERR,329) 
  329 FORMAT(13X,' last year cannot be less than first year') 
      GOTO 198 
   30 WRITE(IOUERR,330) RVALUE 
  330 FORMAT(13X,G10.3,' intercept must be positive') 
      GOTO 198 
   31 WRITE(IOUREP,331) RVALUE 
  331 FORMAT(15X,G10.3,' negative slope for ageing vehicle') 
      GOTO 199 
   32 WRITE(IOUERR,327) IVALUE,IOOFBD,IY1960,IYEND 
      GOTO 198 
   33 WRITE(IOUERR,333) IVALUE,IOOFBD 
  333 FORMAT(13X,I4,5A4,'stringency (10 to 50)') 
      GOTO 198 
   34 WRITE(IOUERR,334) RVALUE,IOOFBD 
  334 FORMAT(13X,F4.0,5A4,'I/M Program compliance rate', 
     *       ' ',' (50% to 100%)') 
      GOTO 198 
   35 WRITE(IOUERR,335) IVALUE,IOOFBD,IY1960,IYEND 
  335 FORMAT(13X,I4,5A4,'MODYR1 (',I4,' to ',I4,')') 
      GOTO 198 
   36 WRITE(IOUERR,336) IVALUE,IOOFBD,IY1960,IYEND 
  336 FORMAT(13X,I4,5A4,'MODYR2 (',I4,' to ',I4,')') 
      GOTO 198 
   37 WRITE(IOUERR,337) 
  337 FORMAT(13X,' MODYR1 cannot be greater than MODYR2') 
      GOTO 198 
   38 WRITE(IOUERR,338) IVALUE,IOOFBD,IYSTART,IYEND 
  338 FORMAT(13X,I4,5A4,'year (',I4,' to ',I4,')') 
      GOTO 198 
   39 WRITE(IOUERR,339) RVALUE 
  339 FORMAT(13X,G10.3,' speed must be positive') 
      GOTO 198 
   40 WRITE(IOUERR,340) RVALUE 
  340 FORMAT(13X,G10.3,' valid ambient temperature', 
     *                    ' is 0-110 deg. (F)') 
      GOTO 198 
   41 WRITE(IOUERR,341) RVALUE,IOOFBD 
  341 FORMAT(13X,G10.3,5A4,'PCCN (0. to 100.)') 
      GOTO 198 
   42 WRITE(IOUERR,342) RVALUE,IOOFBD 
  342 FORMAT(13X,G10.3,5A4,'PCHC (0. to 100.)') 
      GOTO 198 
   43 WRITE(IOUERR,343) RVALUE,IOOFBD 
  343 FORMAT(13X,G10.3,5A4,'PCCC (0. to 100.)') 
      GOTO 198 
   44 WRITE(IOUREP,344) RVALUE 
  344 FORMAT(15X,'Fraction of AC use set by user to be ',F5.3,'.') 
      GOTO 199 
   45 WRITE(IOUREP,345) RVALUE 
  345 FORMAT(15X,'Unacceptable value of',F5.3,' for cloud cover.  ', 
     &           'Fraction of cloud cover ',/,15X, 
     &           'must be between zero and one.  0% Cloud Cover used') 
      GOTO 199 
   46 WRITE(IOUERR,346) RVALUE,IOOFBD 
  346 FORMAT(13X,G10.3,5A4,'AC Demand fraction less than 0') 
      GOTO 198 
   47 WRITE(IOUREP,347) RVALUE 
  347 FORMAT(15X,'Unacceptable value of ',F5.1,' for humidity.  ', 
     &    'Must be between 20. and 528.',/15x,'Default value of 75.0', 
     &    ' used.') 
      GOTO 199 
   48 WRITE(IOUREP,348) NAMVEH(IVALUE) 
  348 FORMAT(13X,' there are no sales for vehicle class ',A8) 
      GOTO 199 
   49 WRITE(IOUREP,349) RVALUE 
  349 FORMAT(15X,G10.3,' MYR sum not = 1. (will normalize)') 
      GOTO 199 
   50 WRITE(IOUERR,350) RVALUE,IOOFBD 
  350 FORMAT(13X,G10.3,5A4,'tampering rate intercept (up to 1.0)') 
      GOTO 198 
   51 WRITE(IOUERR,351) IVALUE,MAXERU 
  351 FORMAT(13X,' equation',I3,' can not be stored.  Only',I3, 
     *  ' equations'/13X,' per pol/veh/reg case are allowed.'/ 
     *  '0',' ') 
      GOTO 198 
   52 WRITE(IOUREP,352) RVALUE 
  352 FORMAT(15X,G10.3,' speed increased to 2.5 mph minimum') 
      GOTO 199 
C 
   53 WRITE(IOUERR,353) NAMVEH(IVALUE),INT(RVALUE) 
  353 FORMAT(15X,'The subroutine RSTLOS.FOR determined a fraction'/ 
     *       15X,'for TIER2 phase in of ',A8,'in model year ',I4,/ 
     *       15X,'which is less than or equal to zero. This case'/ 
     *       15X,'has the 2007 heavy duty vehicle rule in effect.') 
      GOTO 198 
C 
   54 WRITE(IOUERR,354) 
  354 FORMAT(T1,'I/M credits for this program description ', 
     *       'not found.'/'0') 
      GOTO 198 
   55 WRITE(IOUERR,355) 
  355 FORMAT(T1,'Error return reading I/M credits.'/'0') 
      GOTO 198 
   56 WRITE(IOUREP,356) IVALUE 
  356 FORMAT(15X,' Unacceptable sunrise/sunset time of',I3,'.  Sunrise', 
     &           ' and sunset must be between 5 and 9.') 
      GOTO 199 
   57 WRITE(IOUREP,357) IVALUE 
  357 FORMAT(15X,' Unacceptable peak sun time of',I3,'.  Beginning of', 
     &           ' peak sun must be between 9 AM and 1 PM.',/,15x, 
     &           ' End of peak sun must be between 1 and 5 PM.') 
      GOTO 199 
   58 WRITE(IOUREP,358) IVALUE,RVALUE 
  358 FORMAT(15X,' Average miles/day is zero for vehicle class',I2, 
     *               ' and model year',F6.0) 
      GOTO 199 
   59 WRITE(IOUREP,359) IVALUE,(NEWPAR(IPAR,IVALUE),IPAR=1,3) 
  359 FORMAT(15X,' equation',I3,' zeroes all idle coefficients', 
     *   ' (and total) for IR=',I1,', IV=',I1,', IP=',I1) 
      GOTO 199 
   60 WRITE(IOUERR,360) INERR,NAMFNC(IVALUE) 
  360 FORMAT(13X,I4,' exceeds ub of highest myg in',A8) 
      GOTO 198 
   61 WRITE(IOUERR,361) INERR,NAMFNC(IVALUE) 
  361 FORMAT(13X,' default used for',I5,' in index function',A8) 
      GOTO 198 
   62 WRITE(IOUERR,362) IVALUE,IOOFBD 
  362 FORMAT(13X,I4,5A4,'new e.f. region (1 to 2)') 
      GOTO 198 
   63 WRITE(IOUERR,363) RVALUE,IOOFBD 
  363 FORMAT(13X,G10.3,5A4,'(0. <= PCHC+PCCC <= 100.)') 
      GOTO 198 
   64 WRITE(IOUERR,364) RVALUE,IOOFBD 
  364 FORMAT(13X,G10.3,5A4,'(0. <= PCHC+PCCC-PCCN <= 100.)') 
      GOTO 198 
   65 WRITE(IOUERR,365) IVALUE 
  365 FORMAT(13X,' MYR of LDDV not equal to LDGV for JDX=',I2) 
      GOTO 198 
   66 WRITE(IOUERR,366) IVALUE 
  366 FORMAT(13X,' MYR of LDDT not equal to LDGT1 for JDX=',I2) 
      GOTO 198 
   67 WRITE(IOUERR,367) NAMVEH(IVALUE) 
  367 FORMAT(13X,' EFFTP>0. and GSF=0. for vehicle ',A8) 
      GOTO 198 
   68 WRITE(IOUERR,368) NAMVEH(IVALUE) 
  368 FORMAT(13X,' EFFTP<0. and VMTMIX>0. for vehicle ',A8) 
      GOTO 198 
C 
   69 WRITE(IOUERR,369) 
  369 FORMAT( 
     *13X,' Error occurred in reading the', 
     *        ' detergent gasoline parameter record.') 
      GOTO 198 
C 
   70 WRITE(IOUERR,370) RVALUE*100, IOOFBD, NAMDET(IVALUE) 
  370 FORMAT( 
     *13X,' Error occurred in reading the', 
     *        ' detergent gasoline parameter record.'/ 
     *13X,F6.2,5A4,A10) 
      GOTO 198 
   71 WRITE(IOUERR,371)  
  371 FORMAT(13X,'Mode must be R or S in ber record') 
      GOTO 198 
   72 WRITE(IOUERR,372) IVALUE,IOOFBD 
  372 FORMAT(13X,I4,5A4,'I/M vehicle type (1 to 2)') 
      GOTO 198 
   73 WRITE(IOUERR,373) IVALUE,IOOFBD 
  373 FORMAT(13X,I4,5A4,'short test type flag (1 to 3)') 
      GOTO 198 
   74 WRITE(IOUREP,374) IVALUE 
  374 FORMAT(13X,'Data entered in equation ',I3,' for ', 
     * 'running emissions for a 94+ LDG' 
     */13X,'vehicle will not be used in emission ', 
     * 'factor calculations.  Those', 
     */13X,'values can only be changed by placing them ', 
     * 'in an external data ', 
     */13X,'file and reading them into the model using ', 
     * 'a 94P LDG RUNNING EF ', 
     */13X,'labeled input record.') 
      GOTO 199 
   75 WRITE(IOUERR,375) IVALUE,IOOFBD 
  375 FORMAT(13X,I4,5A4,'anti-tampering program start year') 
      GOTO 198 
   76 WRITE(IOUERR,376) IVALUE,IOOFBD 
  376 FORMAT(13X,I4,5A4,'anti-tampering program 1st model year') 
      GOTO 198 
   77 WRITE(IOUERR,377) IVALUE,IOOFBD 
  377 FORMAT(13X,I4,5A4,'anti-tampering program last model year') 
      GOTO 198 
   78 WRITE(IOUERR,378) IVALUE,IOOFBD 
  378 FORMAT(13X,I4,5A4,'anti-tampering vehicle type (1 to 2)') 
      GOTO 198 
   79 WRITE(IOUERR,379) RVALUE,IOOFBD 
  379 FORMAT(13X,G10.3,5A4,'evap ATP effectiveness rate (0. to 1.)') 
      GOTO 198 
   80 WRITE(IOUERR,380) RVALUE,IOOFBD 
  380 FORMAT(13X,G10.3,5A4,'exh ATP effectiveness rate (0. to 1.)') 
      GOTO 198 
   81 WRITE(IOUERR,381) RVALUE,IOOFBD 
  381 FORMAT(13X,G10.3,5A4,'AER matrix file column sum (0. to 1.)') 
      GOTO 198 
   82 WRITE(IOUERR,382) 
  382 FORMAT(13X,'An invalid I/M program type was specified.', 
     *       ' (Use ''T/O'' for', 
     *      /13X,'test-only, ''TRC'' for computerized test-and', 
     *       '-repair, or ''TRM''', 
     *      /13X,'for manual test-and-repair.)') 
      GOTO 198 
   83 WRITE(IOUERR,383)  
  383 FORMAT(13X,'CO2 cannot be calculated when NGVs are selected ',  
     *           'via the commands NGV FRACTION or NGV EF') 
      GOTO 198 
   84 WRITE(IOUERR,384) IVALUE,IOOFBD 
  384 FORMAT(13X,I4,5A4,'onboard VRS vehicle class flag (1 to 2)') 
      GOTO 198 
   85 WRITE(IOUERR,385) RVALUE,NAMPOL(IVALUE) 
  385 FORMAT(13X,G10.3,' is out of bounds for the combined ', 
     *        'effects of Detergent' / 
     *        13X,'Gasoline on',A5,'as entered by the user.') 
      GOTO 198 
   86 WRITE(IOUERR,386) IVALUE,IOOFBD 
  386 FORMAT(13X,I4,5A4,'no. of Stage II phase-in years (1 to 9)') 
      GOTO 198 
   87 WRITE(IOUERR,387) IVALUE,IOOFBD 
  387 FORMAT(13X,I4,5A4,'vrs or rvp controls start year') 
      GOTO 198 
   88 WRITE(IOUERR,388) RVALUE,IOOFBD 
  388 FORMAT(13X,G10.3,5A4,'diesel sales fraction (0. to 1.)') 
      GOTO 198 
   89 WRITE(IOUERR,389) IVALUE,IOOFBD 
  389 FORMAT(13X,I4,5A4,' flag PROMPT (1 to 6)') 
      GOTO 198 
   90 WRITE(IOUREP,390) NAMRVP(IVALUE),RVALUE 
  390 FORMAT(15X,A8,' RVP reset to ',F4.1) 
      GOTO 199 
   91 WRITE(IOUERR,391) RVALUE,IOOFBD,NAMTEM(IVALUE) 
  391 FORMAT(13X,F10.1,5A4,A4,'mum daily temperature') 
      GOTO 198 
   92 WRITE(IOUERR,392) RVALUE 
  392 FORMAT(13X,' max daily temp =',F4.0,' < min daily temp') 
      GOTO 198 
   93 WRITE(IOUERR,393) RVALUE,IOOFBD 
  393 FORMAT(13X,F10.2,5A4,'%TLVMT (set to zero.)') 
      GOTO 198 
   94 WRITE(IOUERR,394) RVALUE 
  394 FORMAT(15X,G10.3,' %TLVMT sum not = 100. (will normalize)') 
      GOTO 198 
   95 WRITE(IOUERR,395) RVALUE 
  395 FORMAT(13X,' ambient temp =',F4.0,' is < daily min temp', 
     *               ' or > daily max temp') 
      GOTO 198 
   96 WRITE(IOUREP,396) RVALUE 
  396 FORMAT(15X,G10.3,' speed reduced to 65 mph maximum') 
      GOTO 199 
   97 WRITE(IOUERR,397) 
  397 FORMAT( 
     *13X,' Error occurred in reading a running', 
     *        ' loss hourly trip distribution fraction.') 
      GOTO 198 
   98 WRITE(IOUREP,398) RVALUE 
  398 FORMAT( 
     *15X,'Diurnal temperature rise (max temp-min temp=',F5.1,')'/ 
     *15X,'is > 40F; diurnal evaporative emission factors will'/ 
     *15X,'be calculated, but may be inaccurate.') 
      GOTO 199 
   99 WRITE(IOUERR,399) RVALUE,IOOFBD,NAMEFF(IVALUE) 
  399 FORMAT(13X,F5.0,5A4,'Stage II efficiency for ',A1,'DGV ', 
     *          '(0 to 100 %)') 
      GOTO 198 
  100 WRITE(IOUERR,400) RVALUE,IOOFBD 
  400 FORMAT(13X,F3.0,5A4,'I/M Program waiver rate (0 to 50%)') 
      GOTO 198 
  101 WRITE(IOUERR,401) IVALUE,IOOFBD 
  401 FORMAT(13X,I4,5A4,'frequency of I/M inspection (1 (annual) ', 
     *                      'to 2 (biennial))') 
      GOTO 198 
  102 WRITE(IOUERR,402) 
  402 FORMAT(13X,'Invalid character in the effectiveness flag', 
     &  ' field of the I/M program', /13X,'descriptive record.', 
     &  ' The field should contain 1, 2, or a blank.') 
      GOTO 198 
  103 WRITE(IOUERR,403) IVALUE,IOOFBD 
  403 FORMAT(13X,I4,5A4,'type of ATP program (1 (Test only)', 
     *'to 2 (Test and Repair))') 
      GOTO 198 
  104 WRITE(IOUERR,404) IVALUE,IOOFBD 
  404 FORMAT(13X,I4,5A4,'frequency of ATP inspection (1 (annual) ', 
     *                      'to 2 (biennial))') 
      GOTO 198 
  105 WRITE(IOUERR,405) RVALUE,IOOFBD 
  405 FORMAT(13X,F4.0,5A4,'ATP compliance rate must be'/ 
     *       ' ',36X,'in the range 0.% to 100.%, inclusive') 
      GOTO 198 
  106 WRITE(IOUERR,406) IVALUE,IOOFBD 
  406 FORMAT(13X,I4,5A4,'ATP disablement inspections (1 (no) ', 
     *                      'to 2 (yes))') 
      GOTO 198 
  107 WRITE(IOUERR,407) 
  407 FORMAT(13X, ' Error reading ATP program description', 
     *                ' in the run data section') 
      GOTO 198 
  108 WRITE(IOUREP,408) 
      IF(ATPFLG.EQ.2) WRITE(IOUREP,908) CRATP 
      IF(IMFLAG.EQ.2) WRITE(IOUREP,909) CRIM(1) 
      IF(IMFLAG.EQ.3) WRITE(IOUREP,910) CRIM(1),CRIM(2) 
      IF(PRSFLG.EQ.2) WRITE(IOUREP,911) CRPRS 
      IF(PRGFLG.EQ.2) WRITE(IOUREP,912) CRPRG 
  408 FORMAT(' ', 15X, 'The user supplied ', 
     *                 'compliance rates are not equal.') 
  908 FORMAT(' ', 15X, 'ATP compliance rate:             ',9X,F5.1) 
  909 FORMAT(' ', 15X, 'I/M program #1 compliance rate:  ',9X,F5.1) 
  910 FORMAT(' ', 15X, 'I/M program #1 compliance rate:  ',9X,F5.1/ 
     *       ' ', 15X, 'I/M program #2 compliance rate:  ',9X,F5.1) 
  911 FORMAT(' ', 15X, 'Pressure Check compliance rate:  ',9X,F5.1) 
  912 FORMAT(' ', 15X, 'Purge Check compliance rate:     ',9X,F5.1) 
      GOTO 199 
  109 QFREQ(1) = IFREQ(1) 
      IF(IFREQ(1).GT.2) QFREQ(1) = 3 
      QFREQ(2) = IFREQ(2) 
      IF(IFREQ(2).GT.2) QFREQ(2) = 3 
C 
      QFQTFN(1) = PRSFQT 
      IF(PRSFQT.GT.2.OR.PRSFLG.EQ.1) QFQTFN(1) = 3 
      QFQTFN(2) = PRGFQT 
      IF(PRGFQT.GT.2.OR.PRGFLG.EQ.1) QFQTFN(2) = 3 
      QAFREQ = ATPFQT 
      IF(ATPFLG.EQ.1.OR.ATPFQT.GT.2) QAFREQ=3 
      IF(IMFLAG.EQ.2) QFREQ(2)=3 
C 
      WRITE(IOUREP,409) NAMFRQ(QAFREQ),(NAMFRQ(QFREQ(I)),I=1,2), 
     *                 (NAMFRQ(QFQTFN(I)),I=1,2) 
  409 FORMAT(' ', 15X, 'The user supplied inspection ', 
     *                 'frequencies are not equal.'/ 
     *       ' ', 15X, 'ATP inspection frequency:                 ',A9/ 
     *       ' ', 15X, 'I/M program #1 inspection frequency:      ',A9/ 
     *       ' ', 15X, 'I/M program #2 inspection frequency:      ',A9/ 
     *       ' ', 15X, 'Pressure Check inspection frequency:      ',A9/ 
     *       ' ', 15X, 'Purge Check inspection frequency:         ',A9) 
      GOTO 199 
C 
  110 QINTYP(1) = INTYP(1) 
      IF(INTYP(1).GT.2) QINTYP(1) = 4 
      QINTYP(2) = INTYP(2) 
      IF(INTYP(2).GT.2.OR.IMFLAG.EQ.2) QINTYP(2) = 4 
C 
      QPGMFN(1) = PRSPGM 
      IF(PRSPGM.GT.2.OR.PRSFLG.EQ.1) QPGMFN(1) = 4 
      QPGMFN(2) = PRGPGM 
      IF(PRGPGM.GT.2.OR.PRGFLG.EQ.1) QPGMFN(2) = 4 
C 
      QAPROG = ATPPGM 
      IF(ATPFLG.EQ.1.OR.ATPPGM.GT.2) QAPROG=4 
C 
      WRITE(IOUREP,410) NAMTYP(2,QAPROG),(NAMTYP(1,QINTYP(I)),I=1,2), 
     *                 (NAMTYP(2,QPGMFN(I)),I=1,2) 
  410 FORMAT(15X,'The user supplied inspection ', 
     *               'types are not equal.'/ 
     *     ' ', 15X, 'ATP inspection type:               ',A20/ 
     *     ' ', 15X, 'I/M program #1 inspection type:    ',A20/ 
     *     ' ', 15X, 'I/M program #2 inspection type:    ',A20/ 
     *     ' ', 15X, 'Pressure Check inspection type:    ',A20/ 
     *     ' ', 15X, 'Purge Check inspection type:       ',A20) 
      GOTO 199 
C 
  111 WRITE(IOUREP,411) RVALUE 
  411 FORMAT(15X,'The input dIesel sulfur level of',F6.1,' ppm exceeds'/ 
     *       15X,'the 2007 HDD Rule diesel sulfur limit of 15 ppm.') 
      GOTO 199 
C 
  112 WRITE(IOUREP,412) 
  412 FORMAT(15X,'Wintertime Reformulated Gasoline Rules Apply') 
          GOTO 199 
C 
  113 WRITE(IOUERR,413) 
  413 FORMAT(15X,'There was an error reading values from the ', 
     & 'program effectiveness',/15X,'record. Please check format.') 
      GOTO 198 
  114 WRITE(IOUREP,414) 
  414 FORMAT(15X,'Purge Check emission benefits assume the use ', 
     *               'of a dynomometer and'/15X,'the IM240 ', 
     *               'transient test procedure driving cycle.') 
      GOTO 199 
C 
  115 WRITE(IOUERR,415) 
  415 FORMAT(13X,'FUEL PROGRAM Option #2 (RFG) cannot be used with'/ 
     &       13X,'an AIR TOXICS run. Please see Section 2.8.10 of'/ 
     &       13X,'the MOBILE6 User Guide.') 
      GOTO 198 
C 
  116 WRITE(IOUERR,416) TRIM(NAMEAT(IVALUE)) 
  416 FORMAT(13X,'The ',A,' record for the', 
     &      /13X,'air toxic program is missing or invalid.', 
     &      /13X,'Please make sure a valid record is', 
     &      /13X,'input for the air toxic program.') 
      GOTO 198 
C 
  117 WRITE(IOUERR,417) IVALUE,IOOFBD 
  417 FORMAT(13X,I4,5A4,'by model year table pick (1 to 2)') 
      GOTO 198 
  118 WRITE(IOUERR,418) 
  418 FORMAT( 
     *15X,'Input of diesel sales fractions not allowed'/ 
     *15X,'in Run Data Section.') 
      GOTO 198 
  119 WRITE(IOUREP,419) 
  419 FORMAT( 
     *15X,'Zero tampering (TAMFLG=0) is for EPA in-house use ONLY.') 
      GOTO 195 
  120 WRITE(IOUREP,420) IVALUE,IOOFBD 
  420 FORMAT(15X,I4,5A4,'evaluation month of a reformulated gas'/ 
     *15X,'area. Must be 1 (Jan) or 7 (July). Season will default', 
     *' to summer.') 
      GOTO 199 
  121 WRITE(IOUREP,421) 
  421 FORMAT( 
     *15X,'NO emission credit from evaporative control', 
     *        ' system inspection'/ 
     *15X,'will be given unless gas cap inspection also required.') 
      GOTO 199 
  122 WRITE(IOUERR,422) 
  422 FORMAT( 
     *13X,' Error occurred in reading mileage', 
     *        ' accrual rates.') 
      GOTO 198 
  123 WRITE(IOUERR,423) 
  423 FORMAT( 
     *13X,' Error occurred in reading registration distributions.') 
      GOTO 198 
  124 WRITE(IOUERR,424) 
  424 FORMAT( 
     *13X,' Error occurred in reading a vertically input flag.') 
      GOTO 198 
  125 WRITE(IOUERR,425) 
  425 FORMAT( 
     *13X,' Error occurred in reading a horizontally input flag.') 
      GOTO 198 
  126 WRITE(IOUERR,426) 
  426 FORMAT( 
     *13X,' Error occurred in reading a', 
     *        ' replacement tampering rate.') 
      GOTO 198 
  127 WRITE(IOUERR,427) 
  427 FORMAT( 
     *13X,' Error occurred in reading the VMTmix record.') 
      GOTO 198 
  128 WRITE(IOUERR,428) 
  428 FORMAT( 
     *13X,' Error occurred in reading the', 
     *        ' replacement emission factor rates.') 
      GOTO 198 
C 
  129 WRITE(IOUERR,429) VEHICLEGROUP(IVALUE) 
  429 FORMAT(13X,'The number of vehicle flags entered ', 
     *13X,'for the I/M VEHICLES command in the ',/ 
     *13X, A19,/ 
     *13X,'vehicle group is incorrect.'/) 
      GOTO 198 
C 
  130 WRITE(IOUERR,430) 
  430 FORMAT( 
     *13X,' Error occurred in reading the', 
     *        ' Stage II or Onboard VRS input record.') 
      GOTO 198 
  131 WRITE(IOUERR,431) 
  431 FORMAT( 
     *13X,' Error occurred in reading the', 
     *        ' local area parameter record.') 
      GOTO 198 
  132 WRITE(IOUERR,432) 
  432 FORMAT( 
     *13X,' Error occurred in reading the', 
     *        ' scenario descriptive record.') 
      GOTO 198 
  133 WRITE(IOUERR,433) 
  433 FORMAT( 
     *13X,' Error occurred in reading the', 
     *        ' additional correction factor record.') 
      GOTO 198 
  134 WRITE(IOUERR,434) 
  434 FORMAT( 
     *13X,' Error occurred in reading the', 
     *        ' BY-Model-Year selection flags.') 
      GOTO 198 
  135 WRITE(IOUREP,435) 
  435 FORMAT( 
     *13X,' All effects of the 1990 Clean Air Act Amendments', 
     *' have been disabled.') 
      GOTO 199 
  136 WRITE(IOUERR,436) 
  436 FORMAT( 
     *13X,' Error occurred in reading the cold CO parameters.') 
      GOTO 198 
  137 WRITE(IOUERR,437) 
  437 FORMAT( 
     *13X,' Error occurred in reading the', 
     *        ' oxygenated fuels descriptive record.') 
      GOTO 198 
  138 WRITE(IOUERR,438) 
  438 FORMAT( 
     *13X,' Error occurred in reading diesel sales fractions.') 
      GOTO 198 
  139 WRITE(IOUERR,439) IVALUE,IOOFBD,IYEND 
  439 FORMAT(13X,I4,5A4,'ETP phase-in year (1995 to ',I4,')') 
      GOTO 198 
C 
  140 WRITE(IOUERR,440) RVALUE, IOOFBD 
  440 FORMAT(13X,G10.3,5A4,' aromatic percentage (0.0 to 100.0)') 
      GOTO 198 
C 
  141 WRITE(IOUERR,441) RVALUE, IOOFBD 
  441 FORMAT(13X,G10.3,5A4,' benzene percentage (0.0 to 100.0)') 
      GOTO 198 
C 
  142 WRITE(IOUERR,442) RVALUE, IOOFBD 
  442 FORMAT(13X,G10.3,5A4,' E200 percentage (0.0 to 100.0)') 
      GOTO 198 
C 
  143 WRITE(IOUERR,443) RVALUE, IOOFBD 
  443 FORMAT(13X,G10.3,5A4,' E300 percentage (0.0 to 100.0)') 
      GOTO 198 
C 
  144 WRITE(IOUERR,444) RVALUE, IOOFBD 
  444 FORMAT(13X,G10.3,5A4,' olefin percentage (0.0 to 100.0)') 
      GOTO 198 
C 
  145 WRITE(IOUERR,445) IVALUE,IOOFBD 
  445 FORMAT(13X,I4,5A4,'evaluation month 1 (Jan) or 7 (July)') 
      GOTO 198 
  146 WRITE(IOUREP,446) 
  446 FORMAT( 
     *13X,' Error occurred in reading the hourly ', 
     *        'output filename. The user supplied file already'/ 
     *        ' exists. No hourly output will be produced ', 
     *        'for this scenario.') 
      GOTO 199 
  147 WRITE(IOUERR,447) 
  447 FORMAT( 
     *13X,'User supplied market share of oxygenated fuels has ', 
     *        'to be 100% for'/13X,'a reformulated gas area.') 
      GOTO 198 
  148 WRITE(IOUERR,448) 
  448 FORMAT( 
     *13X,'User supplied market share weighted average of ', 
     *        'oxygen content has'/13X,'to be at least 2.1% ', 
     *        'for a reformulated gas area.') 
      GOTO 198 
  149 WRITE(IOUERR,449) IOOFBD 
  449 FORMAT( 
     *13X,'Value ',5A4,' ') 
      GOTO 198 
  150 WRITE(IOUERR,450) RVALUE,IOOFBD 
  450 FORMAT(13X,F5.0,5A4,'Pass/Purge/Pressure rate (0 to 100%).') 
      GOTO 198 
  151 WRITE(IOUERR,451) RVALUE 
  451 FORMAT(13X,'The Pass/Purge/Pressure rate sum (',F4.0,') ', 
     *               'does not equal 100%') 
      GOTO 198 
  152 WRITE(IOUREP,452) 
  452 FORMAT( 
     *15X, 'Normal and high average emission rates are',/ 
     *15X, 'the same. All vehicles have been treated as',/ 
     *15X, 'normal emitters.' ) 
      GOTO 199 
C 
  153 WRITE(IOUERR,453) IVALUE, IOOFBD 
  453 FORMAT(13X,I2,5A4,'the I/M program number (1 to 7).') 
      GOTO 198 
C 
  154 WRITE(IOUREP,454) IVALUE 
  454 FORMAT(15X,'User input of emission line parameters for ' 
     &  /15X,'light duty gasoline vehicles (1981-1993' 
     &  /15X,'model years) with the BASIC EFS label ' 
     &  /15X,'(entry #',I3,') is out of the model year' 
     &  /15X,'range (pre-1981 only). This entry will' 
     &  /15X,'be ignored. It must be entered with label' 
     &  /15X,'1981-1993 LDG EFS.') 
C 
      GOTO 199 
C 
  155 WRITE(IOUREP,455) 
  455 FORMAT( 
     *15X,'User input of Onboard Refueling Vapor Recovery',/ 
     *15X,'parameters has been ignored.') 
C 
      GOTO 199 
C 
  156 WRITE(IOUREP,456) 
  456 FORMAT( 
     *15X,'Onboard Refueling Vapor Recovery Regulation (1994)',/ 
     *15X,'effects have been disabled by the user.') 
C 
      GOTO 199 
C 
C  Aug-18-1994 @ CSC-bsg Request 2-462 Added QUITER messages 157-160, for  
C  non-FTP related error messages 
C  Sep-14-1994 @ CSC-bsg Request 2-467 New Non-FTP Effects message and the 
C  elimination of error message 157 
C 
  157 WRITE(IOUERR,457) RVALUE, IOOFBD 
  457 FORMAT(13X,F7.2,5A4,' the effectiveness of test-and-repair', 
     *      /13X,'I/M programs (0.00 to 1.00)') 
      GOTO 198 
  158 WRITE(IOUREP,458) IVALUE 
  458 FORMAT( 
     *15X,'Exhaust emission values include additional effects of ', 
     *        'in-use'/ 
     *15X,'driving behavior.  Model years beginning with ',I4, 
     *        ' will '/ 
     *15X,'have reduced effects due to a new certification ', 
     *        'procedure.') 
C 
      GOTO 199 
C 
  159 WRITE(IOUERR,459) IVALUE,IOOFBD 
  459 FORMAT(13X,I4,5A4,'First model year which will '/ 
     *13X,'not have non-FTP effects.') 
C 
      GOTO 198 
C 
  160 WRITE(IOUERR,460) 
  460 FORMAT( 
     *13X,'Error occurred in reading in-use driving behavior '/ 
     *13X,'record.') 
      GOTO 198 
C 
C 
  161 WRITE(IOUERR,461) 
  461 FORMAT( 
     *13X,'Error occurred in reading PPEFF.D file.') 
      GOTO 198 
C 
C  
  162 WRITE(IOUERR,462) 
  462 FORMAT(13X,'An invalid I/M inspection test type was specified.') 
      GOTO 198 
C       
C  Message 163 if error occurred in reading TTC.D file  
C  
  163 WRITE(IOUERR,463) 
  463 FORMAT( 
     *13X,'ERROR occurred in reading TTC.D file.') 
      GOTO 198 
C       
C  New message 164 if error occurred in reading "I/M Control Flag Record" 
C  
  164 WRITE(IOUERR,464) 
  464 FORMAT( 
     *13X,'ERROR occurred in I/M Control Flag Record') 
      GOTO 198 
C  
C  New message 165 if value out of bounds in the TTC Record"   
C 
  165 WRITE(IOUREP,465) RVALUE 
  465 FORMAT(15X,G7.2,' is out of bounds for the Effectiveness ', 
     *      'in the TTC record.',/,15X,'Reset to 100%') 
      GOTO 199 
C       
C  New message 166 if I/M + TTC credit exceed 100% 
C  
  166 WRITE(IOUREP,466) 
  466 FORMAT(15X,'I/M and TTC credit exceed 100%.', 
     *      ' Credit will be reset to 100%') 
      GOTO 199 
C       
C  New message 167: ATPUER is out of bounds.  
C  
  167 WRITE(IOUERR,467) RVALUE, IOOFBD  
  467 FORMAT( 
     *13X,F7.2, 5A4,' the effectiveness of test-and-repair',/ 
     *13X,'ATP programs (0.00 to 1.00)') 
      GOTO 198 
C       
C  New message 168: PPGUER is out of bounds.  
C  
  168 WRITE(IOUERR,468) RVALUE, IOOFBD 
  468 FORMAT( 
     *13X,F7.2, 5A4,' the effectiveness of test-and-repair',/ 
     *' ',12x,'purge/pressure programs (0.00 to 1.00)') 
      GOTO 198 
C       
C  New message 169: DGPIEY is out of bounds.  
C  
  169 WRITE(IOUERR,469) IVALUE, IOOFBD 
  469 FORMAT( 
     *13X,I4, 5A4,' the detergent gas phase-in end year ',/ 
     *13X,'(this should be later than 1995.)') 
      GOTO 198 
C 
  170 WRITE(IOUERR,470) IVALUE, INT(RVALUE) 
  470 FORMAT( 
     *13X,'Negative technology fraction calculated by PMFRACA'/ 
     *13X,'for Vehicle Class',I3,' and Model Year',I5,'.') 
      GOTO 198 
C 
  171 WRITE(IOUERR,471) 
  471 FORMAT( 
     &13X,'PCLEFT was called for a TECHIV+ vehicle.') 
      GOTO 198 
C 
  172 WRITE(IOUERR,472) IVALUE 
  472 FORMAT( 
     *15X,'FORTRAN error ',I8,' reading a numeric input field.') 
      GOTO 198 
C 
  173 WRITE(IOUERR,473) IVALUE,IOOFBD 
  473 FORMAT(13X,I4,5A4,' expanded output table pick (1 to 2)') 
      GOTO 198 
C 
  174 WRITE(IOUERR,474) IVALUE,IOOFBD 
  474 FORMAT(13X,I4,5A4,' vehicle type (1 to 28) used',/ 
     * 15X,'to supply facility vehicle mileage fractions.') 
      GOTO 198 
C 
  175 WRITE(IOUERR,475) RVALUE,IOOFBD 
  475 FORMAT(13X,F5.3,5A4,'a LEV implementation phase-in', 
     *       ' value (0.0 to 1.0).') 
      GOTO 198 
C 
  176 WRITE(IOUERR,476) RVALUE 
  476 FORMAT(13X,G10.3,' The sum of the items in LEVIMP must', 
     *       ' be 1.+/-0.0001') 
      GOTO 198       
C 
  177 WRITE(IOUERR,477) 
  477 FORMAT( 
     *'+',12X,' Error opening external file for User Supplied', 
     *        ' Sulfur Correction Factors.') 
      GOTO 198 
C 
  178 WRITE(IOUERR,478) IVALUE, IOOFBD 
  478 FORMAT(13X,I4,5A4,' vehicle type in EFCVT.') 
      GOTO 198 
C 
  179 WRITE(IOUERR,479) IVALUE, IOOFBD 
  479 FORMAT(13X,I4,5A4,' I/M grace period (1 to 25).') 
      GOTO 198 
C 
  180 WRITE(IOUERR,480) 
  480 FORMAT(13X,'The number of RSD test failures must be ', 
     *               'greater than zero.') 
      GOTO 198 
C 
  181 WRITE(IOUERR,481) RVALUE, IOOFBD 
  481 FORMAT(13X,G10.3,5A4,' diesel sulfur content (0.01 to 5000.)') 
      GOTO 198 
C 
  182 WRITE(IOUERR,482) RVALUE, IOOFBD 
  482 FORMAT(13X,G10.3,5A4,' particle size cutoff (1. to 10.)') 
      GOTO 198 
C 
  183 WRITE(IOUERR,483) RVALUE, IOOFBD 
  483 FORMAT(13X,G10.3,5A4,' RSD effectiveness (0.0 to 100.0)') 
      GOTO 198 
C 
  184 WRITE(IOUERR,484) TRIM(NAMEPART(IVALUE)) 
  484 FORMAT(13X,'The ',A,' record for the', 
     &      /13X,'particulate program is missing or invalid.', 
     &      /13X,'Please make sure a valid record is', 
     &      /13X,'input for the particulate program.') 
      GOTO 198 
C 
  185 WRITE(IOUERR,485) IVALUE, IOOFBD 
  485 FORMAT(13X,I4,5A4,' I/M exemption age (1 to 25).') 
      GOTO 198 
C 
  186 WRITE(IOUERR,486) IVALUE, IOOFBD 
  486 FORMAT(13X,I4,5A4,' for Weekday=1 or Weekend=2.') 
      GOTO 198 
C 
  187 WRITE(IOUERR,487) IVALUE 
  487 FORMAT(13X,'One or more items from the I/M PROGRAM', 
     &       ' record for I/M program', 
     &      /13X,'number ',I1,' are missing or invalid. Please', 
     &       ' check to see that a', 
     &      /13X,'valid I/M PROGRAM record was supplied for', 
     &       ' that I/M program.') 
      GOTO 198 
C 
  188 WRITE(IOUERR,488) IVALUE 
  488 FORMAT(13X,'Cutpoints must be supplied for I/M program', 
     &       ' ',I1,' because the test', 
     &      /13X,'type for that program is IM240. Make sure a', 
     &       ' valid I/M CUTPOINTS', 
     &      /13X,'record is included in the inputs for that', 
     &       ' I/M program.') 
      GOTO 198 
C 
  189 WRITE(IOUERR,489) NAMEIM(IFIX(RVALUE)),IVALUE 
  489 FORMAT(13X,'The ',A16,' record for I/M', 
     &       ' program ',I1,' is missing or', 
     &      /13X,'invalid. Please make sure a valid record', 
     &       ' is included in the', 
     &      /13X,'inputs for that I/M program.') 
      GOTO 198 
C 
  190 WRITE(IOUERR,490) IVALUE 
  490 FORMAT(13X,'The RSD FAILURE LIMIT, RSD TESTS PERFORMED', 
     &       ' and RSD TOTAL VEHICLES' 
     &      /13X,'records are all required when RSD coverage', 
     &       ' option 1 is selected.', 
     &      /13X,'One or more of these records is missing in', 
     &       ' the inputs for I/M', 
     &      /13X,'program ',I1,'.') 
      GOTO 198 
C 
C After this, new messages begin with message numbers starting at 500. 
C 
  501 WRITE(IOUERR,701) IVALUE 
  701 FORMAT(13X,'RSD coverage options 1 and 2 have both', 
     &       '  been selected for', 
     &      /13X,'I/M program ',I1,'. Only one of these', 
     &       ' options may be selected for', 
     &      /13X,'each I/M program that is modeled.') 
      GOTO 198 
C 
  502 WRITE(IOUERR,702) IVALUE 
  702 FORMAT(13X,'The RSD NUMBER FAILURES and RSD NUMBER', 
     &       ' VEHICLES records are', 
     &      /13X,'both required when RSD coverage option 3 is', 
     &       ' selected. One of', 
     &      /13X,'these records is missing from the inputs for', 
     &       ' I/M program ',I1,'.') 
      GOTO 198 
C 
  503 WRITE(IOUERR,703) IVALUE 
  703 FORMAT(13X,'The program type for I/M program ',I1, 
     &       ' is RSDONLY but no RSD coverage', 
     &      /13X,'option has been selected for that program.') 
      GOTO 198 
C 
  504 WRITE(IOUERR,704) IVALUE 
  704 FORMAT(13X,'Coverage option 1 or coverage option 2 was', 
     &       ' selected for I/M', 
     &      /13X,'program ',I1,', but no RSD CO or NOX', 
     &       ' CUTPOINT record was included in', 
     &      /13X,'the input for that program. One or', 
     &       ' both of these records must', 
     &      /13X,'be provided or the coverage option', 
     &       ' records should be deleted.') 
      GOTO 198 
C 
  505 WRITE(IOUERR,705) IVALUE 
  705 FORMAT(13X,'An RSD CO CUTPOINT record was provided for', 
     &       ' I/M program ',I1,', but', 
     &      /13X,'no RSD coverage option was selected. When an', 
     &       ' RSD CO CUTPOINT record is', 
     &      /13X,'provided, coverage options 1 or 2 must also', 
     &       ' be selected to model an RSD', 
     &      /13X,'dirty screening program.') 
      GOTO 198 
C 
  506 WRITE(IOUERR,706) IVALUE 
  706 FORMAT(13X,'Coverage option 3 was selected and an RSD CO', 
     &       ' CUTPOINT record', 
     &      /13X,'was provided for I/M program ',I1,'. The RSD', 
     &       ' CO CUTPOINT record', 
     &      /13X,'can not be used with coverage option 3.') 
      GOTO 198 
C 
  507 WRITE(IOUERR,707) IVALUE 
  707 FORMAT(13X,'A clean screening RSD program using coverage', 
     &       ' options 1 or 2 has', 
     &      /13X,'been implied for I/M program ',I1,', but no', 
     &       ' RSD NOX CUTPOINT record', 
     &      /13X,'has been provided for the clean screening', 
     &       ' program.') 
      GOTO 198 
C 
  508 WRITE(IOUERR,708) IVALUE 
  708 FORMAT(13X,'An RSD NOX CUTPOINT record was provided for', 
     &       ' I/M program ',I1,', but no', 
     &      /13X,'clean screening coverage option was selected', 
     &       ' for that program. Coverage option 1 or', 
     &      /13X,'coverage option 2 should be selected along', 
     &       ' with the RSD NOX CUTPOINT record.') 
      GOTO 198 
C 
  509 WRITE(IOUERR,709) IVALUE 
  709 FORMAT(13X,'Both the RSD NOX CUTPOINT record and the RSD', 
     &       ' LOSSES records have', 
     &      /13X,'been provided for I/M program ',I1,'. One of', 
     &       ' these records should be', 
     &      /13X,'deleted. The RSD LOSSES record may be', 
     &       ' used alone or the RSD', 
     &      /13X,'NOX CUTPOINT record may be used with coverage', 
     &       ' options 1 or 2.') 
      GOTO 198 
C 
  510 WRITE(IOUERR,710) IVALUE 
  710 FORMAT(13X,'The RSD LOSSES record is included for I/M', 
     &       ' program ',I1,' and coverage', 
     &      /13X,'option 3 is selected for the dirty screening', 
     &       ' program. Coverage', 
     &      /13X,'options 1 or 2 should not also be selected', 
     &       ' when the RSD LOSSES', 
     &      /13X,'record is included for the clean screening', 
     &       '  program.') 
      GOTO 198 
C 
  511 WRITE(IOUERR,711) RVALUE,IOOFBD 
  711 FORMAT(13X,F6.0,5A4,'TTC effectiveness rate', 
     &       ' (0.0 to 100.0).') 
      GOTO 198 
C 
  512 WRITE(IOUERR,712) 
  712 FORMAT(12X,'A valid AVERAGE SPEED command option was not used.') 
      GOTO 198 
C 
  513 WRITE(IOUERR,713) IVALUE 
  713 FORMAT(13X,'The RSD clean screening HC loss for I/M', 
     &       ' program ',I1,' must be', 
     &      /13X,'greater than zero.') 
      GOTO 198 
C 
  514 WRITE(IOUREP,714) RVALUE 
  714 FORMAT(13X,'The combined freeway and ramp average speed entered' 
     *      /12X,'cannot be less than ',F4.1,' miles per hour.' 
     *      /12X,'The average speed will be reset to this value.') 
      GOTO 199 
C 
  515 WRITE(IOUREP,715) RVALUE 
  715 FORMAT(12X,'The combined freeway and ramp average speed entered' 
     *      /12X,'cannot be greater than ',F4.1,' miles per hour.' 
     *      /12X,'The average speed will be reset to this value.') 
      GOTO 199 
C 
  516 WRITE(IOUERR,716) IVALUE 
  716 FORMAT(13X,'More than two I/M programs are defined.') 
      GOTO 198 
C 
  517 WRITE(IOUREP,717) IVALUE 
  717 FORMAT(13X,'I/M program ',I1,' operates on a selective', 
     &       ' basis only (as a result', 
     &      /13X,'of COIM or RSD testing). The program', 
     &       ' frequency data supplied', 
     &      /13X,'for that program will be ignored. The grace', 
     &       ' period and the', 
     &      /13X,'exemption age, if entered, will also be', 
     &       ' ignored.') 
      GOTO 199 
C 
  518 WRITE(IOUERR,718) IVALUE 
  718 FORMAT(13X,'Selective I/M (RSD only) was specified for', 
     &       ' I/M program ',I1,' but no' 
     &      /13X,'RSD program was included for that I/M', 
     &       ' program. An RSD program', 
     &      /13X,'must be included if the RSD ONLY option', 
     &       ' is specified.') 
      GOTO 198 
C 
  519 WRITE(IOUERR,719) IVALUE 
  719 FORMAT(13X,'Selective I/M (COIM only) was specified for', 
     &       ' I/M program ',I1,' but no' 
     &      /13X,'COIM program was included for that I/M', 
     &       ' program. A COIM program', 
     &      /13X,'must be included if the COIM ONLY option', 
     &       ' is specified.') 
      GOTO 198 
C 
  520 WRITE(IOUREP,720) 
  720 FORMAT(15X,' The NO EVP SAWTOOTH command cannot be used when'/ 
     &       15X,' Evap I/M programs exist. The NO EVP SAWTOOTH'/ 
     &       15X,' command has been disabled for this scenario.'/) 
      GOTO 199 
C 
  521 WRITE(IOUREP,721) IVALUE 
  721 FORMAT(13X,'A grace period was specified for I/M', 
     &       ' program ',I1,', and the', 
     &      /13X,'program covers pre-1981 model year', 
     &       ' vehicles. The grace', 
     &      /13X,'period will not be used in the calculation', 
     &       ' of I/M benefits', 
     &      /13X,'for pre-1981 model year vehicles.') 
      GOTO 199 
C 
  522 WRITE(IOUERR,722)  RVALUE, IOOFBD, OBDTYPE(IVALUE) 
  722 FORMAT(' ',12X,F7.2,5A4,A16,' (0.0 to 100.0).') 
      GOTO 198 
C 
  523 WRITE(IOUREP,723) IVALUE 
  723 FORMAT(13X,'OBD I/M tests are applicable for post-1993 model', 
     &      /13X,'year vehicles only. The initial vehicle model', 
     &      /13X,'year for the OBD I/M test (I/M program #',I1,')', 
     &      /13X,'has been reset to 1994.') 
      GOTO 199 
C 
  524 GOTO 199 
C 
  525 GOTO 199 
C 
  526 WRITE(IOUREP,726) 
  726 FORMAT(1X,/,12X,'User Input DISABLE 1990 CAAA not available ', 
     &       'in current version of MOBILE6.',/, 
     &   12X,'Input for this option has been disabled.',/) 
      GOTO 199 
C 
  527 WRITE(IOUERR,727) RVALUE 
  727 FORMAT(1X,/,12X,'Tier2 Phase-in Fraction of ',F7.3,' is less ', 
     &       'than zero.',/,12X,'Phase-in fraction must be between ', 
     &       'zero and one.',/,12X, 
     &       'Default Tier2 Phase-in Fractions will be used.' ) 
      GO TO 198 
C 
  528 WRITE(IOUERR,728) RVALUE 
  728 FORMAT(1X,/,12X,'Tier2 phase-in fraction of ',F7.3,' is greater ', 
     &       'than one.',/,12X,'Phase-in fraction must be between ', 
     &       'zero and one.',/,12X, 
     &       'Default Tier2 Phase-in Fractions will be used.' ) 
      GO TO 199 
C 
  529 WRITE(IOUERR,729) RVALUE 
  729 FORMAT(1X,/,12X,'Sum of 12 Tier2 phase-in bin fractions is ', 
     &       'equal to ',F7.3,/,12X,'Sum of 12 Tier2 bin fractions ', 
     &       'must be equal to 1.000 +/- 0.001',/,12X, 
     &       'Default Tier2 Phase-in Fractions will be used.' ) 
      GO TO 199 
C 
  530 WRITE(IOUREP,730) RVALUE, IVALUE 
  730 FORMAT(1X,/,12X,'Tier2 50K Certification Standard is ', 
     &       'equal to ',F7.3,' for pollutant',I2,'.',/, 
     &   12X,'Certification Standards must be greater than zero.', 
     & /,12X,'Default Tier2 50K Certification Standards will be used.') 
      GO TO 199 
C 
  531 WRITE(IOUERR,731) IVALUE 
  731 FORMAT( 
     *15X,'I/M Program Description',I2,' :',/, 
     *15X,'MOBILE6 Does NOT Support Evap OBD or OBD & GC I/M Programs '/ 
     *15X,'for Pre-2007 Model Years for HDGV.  Choose a Different '/ 
     *15X,'Test Type for these Model Year HDGVs.') 
      GOTO 198 
C 
  532 WRITE(IOUERR,732) IVALUE,IOOFBD 
  732 FORMAT(13X,I4,5A4,'NGV EF Vehicle Type Index.') 
      GOTO 198 
C 
  533 WRITE(IOUERR,733) IVALUE,IOOFBD 
  733 FORMAT(13X,I4,5A4,'NGV EF Pollutant Type Index.') 
      GOTO 198 
C 
  534 WRITE(IOUERR,734) IVALUE,IOOFBD 
  734 FORMAT(13X,I4,5A4,'NGV EF Run/Start Toggle.') 
      GOTO 198 
C 
  535 WRITE(IOUERR,735) RVALUE,IOOFBD 
  735 FORMAT(13X,F6.0,5A4,'NGV Fleet Fraction must be between', 
     &       ' (0.0 to 1.0).') 
      GOTO 198 
C 
  536 WRITE(IOUERR,736) IVALUE,IOOFBD 
  736 FORMAT(13X,I4,5A4,'Alternate EF Return Code Error.') 
      GOTO 198 
C 
  537 WRITE(IOUERR,737) IVALUE,IOOFBD 
  737 FORMAT(13X,I4,5A4,'Alternate EF V-Type Error.') 
      GOTO 198 
C 
  538 WRITE(IOUERR,738) IVALUE,IOOFBD 
  738 FORMAT(13X,I4,5A4,'Alternate EF Pollutant Type Error.') 
      GOTO 198 
C 
  539 WRITE(IOUERR,739) IVALUE,IOOFBD 
  739 FORMAT(13X,I4,5A4,'Alternate EF Emission Type Error.') 
      GOTO 198 
C 
  540 WRITE(IOUERR,740) IVALUE,IOOFBD 
  740 FORMAT(13X,I4,5A4,'Alternate EF Emitter Type Error.') 
      GOTO 198 
C 
  541 WRITE(IOUERR,741) IVALUE,IOOFBD 
  741 FORMAT(13X,I4,5A4,'Alternate EF Standard Type Error.') 
      GOTO 198 
C 
  542 WRITE(IOUERR,742) IVALUE,IOOFBD 
  742 FORMAT(13X,I4,5A4,'Alternate EF ZML < 0.0.') 
      GOTO 198 
C 
  543 WRITE(IOUERR,743) IVALUE,IOOFBD 
  743 FORMAT(13X,I4,5A4,'Alternate EF First Det. < 0.0.') 
      GOTO 198 
C 
  544 WRITE(IOUERR,744) IVALUE,IOOFBD 
  744 FORMAT(13X,I4,5A4,'Alternate EF Second Det. < 0.0.') 
      GOTO 198 
C 
  545 WRITE(IOUERR,745) 
  745 FORMAT( 
     *15X,'MOBILE6 Does NOT Support OBD I/M Programs for '/ 
     *15X,'Pre-1996 Model Years.  Choose a Different Test'/ 
     *15X,'Type for these Model Years') 
      GOTO 198 
C 
  546 WRITE(IOUERR,746) 
  746 FORMAT( 
     *15X,'MOBILE6 Does NOT Support Exhaust I/M Programs for '/ 
     *15X,'1996 and Later Model Years.  Choose a Different Test'/ 
     *15X,'Type for these Model Years') 
      GOTO 198 
C 
  547 CONTINUE 
        WRITE(IOUERR,747,IOSTAT=IOS,ERR=196) 
  747   FORMAT('*** Error: No scenario was found in this run.', 
     &        /T12,'No output will be produced.') 
      GO TO 198 
C 
  548 WRITE(IOUERR,748) IVALUE,IOOFBD 
  748 FORMAT(13X,I4,5A4,'Hot Soak Activity Fraction Minute Index.') 
      GOTO 198 
C 
  549 WRITE(IOUERR,749) RVALUE,IOOFBD,IVALUE 
  749 FORMAT(13X,F6.4,5A4,'Hot Soak Activity Fraction Summation', 
     &                    ' for Hour Period',I3,'.') 
      GOTO 198 
C 
  550 WRITE(IOUERR,750) IVALUE,IOOFBD 
  750 FORMAT(13X,I4,5A4,'Diurnal Soak Activity Fraction Hour Index.') 
      GOTO 198 
C 
  551 WRITE(IOUERR,751) RVALUE,IOOFBD 
  751 FORMAT(13X,F6.4,5A4,'Diurnal Soak Activity Fraction must', 
     &                    ' be in the range (0.0 to 1.0).') 
      GOTO 198 
C 
  552 WRITE(IOUERR,752) RVALUE,IOOFBD,IVALUE 
  752 FORMAT(13X,F6.4,5A4,'Fixed P/P Rate Summation.'/ 
     &                    ' Using Default/Calculated Values') 
      GOTO 198 
C 
  553 WRITE(IOUERR,753) 
  753 FORMAT( 
     &13X,' Error occurred in reading the Fixed P/P Rates record.') 
      GOTO 198 
C 
  554 WRITE(IOUERR,754) RVALUE 
  754 FORMAT(13X,' Sulfur concentration must range ', 
     *   'between 30 ppm and 1000 ppm. Entered value is ',F10.1) 
      GOTO 198 
C 
  555 WRITE(IOUERR,755) RVALUE 
  755 FORMAT(13X,' Sulfur maximum concentration ',F10.1,' is ', 
     *   'less than the average sulfur entered.') 
      GOTO 198 
C 
  556 WRITE(IOUERR,756) 
  756 FORMAT( 
     *15X,'MOBILE6 Does NOT Support Evap OBD or OBD & GC I/M Programs '/ 
     *15X,'for Pre-1996 Model Years.  Choose a Different Test Type'/ 
     *15X,'for these Model Years') 
      GOTO 198 
C 
  557 WRITE(IOUERR,757) RVALUE 
  757 FORMAT(1X,/,12X,'Tier2 EVAP Phase-in Fraction of ',F7.3,' is ', 
     &       'less than zero or greater than one.',/,12X, 
     &       'Phase-in fraction must be between zero and one.') 
      GO TO 198 
C 
  558 WRITE(IOUERR,758) RVALUE,IOOFBD 
  758 FORMAT(13X,F4.0,5A4,'GLL reduction rate (0 to 100%).') 
      GOTO 198 
C 
  559 WRITE(IOUERR,759) IVALUE 
  759 FORMAT(13X,'One or more items from the I/M PROGRAM', 
     &       ' record for Evap I/M program', 
     &      /13X,'number ',I1,' are missing or invalid. Please', 
     &       ' check to see that a', 
     &      /13X,'valid I/M PROGRAM record was supplied for', 
     &       ' that Evap I/M program.') 
      GOTO 198 
C 
  560 WRITE(IOUERR,760) NAMEIM(IFIX(RVALUE)),IVALUE 
  760 FORMAT(13X,'The ',A16,' record for Evap I/M', 
     &       ' program ',I1,' is missing or', 
     &      /13X,'invalid. Please make sure a valid record', 
     &       ' is included in the', 
     &      /13X,'inputs for that Evap I/M program.') 
      GOTO 198 
C 
  561 WRITE(IOUERR,761) RVALUE,IOOFBD 
  761 FORMAT( 
     *15X,F3.0,5A4,'Evap I/M Program Fractional Benefit'/ 
     *15X,'must be in the range (0.0 to 1.0)') 
      GOTO 198 
C 
  562 WRITE(IOUERR,762) 
  762 FORMAT( 
     *13X,' Error occurred in reading diurnal soak activity fractions.') 
      GOTO 198 
C 
  563 WRITE(IOUERR,763) 
  763 FORMAT( 
     *13X,' Error occurred in reading ASM I/M Test Credits', 
     *    ' from ASMDATA.D.') 
      GOTO 198 
C 
  564 WRITE(IOUERR,764) 
  764 FORMAT( 
     *13X,' Error opening external data file ASMDATA.D.') 
      GOTO 198 
C 
  565 WRITE(IOUERR,765) 
  765 FORMAT( 
     *15X,' The NO Clean Air Act Flag does NOT Support an '/ 
     *15X,' OBD I/M Program'/) 
      GOTO 198 
C 
  566 WRITE(IOUREP,766) 
  766 FORMAT( 
     *15X,' ** Warning ** This type of I/M Program  '/ 
     *15X,' Probably Did Not Exist Prior to the CAA of 1990'/) 
      GOTO 199 
C 
  567 WRITE(IOUERR,767) RVALUE 
  767 FORMAT(13X,'The Defeat Device Rebuild Fraction ',F6.3, 
     &       ' is out of bounds. It must be between 0.01 and 0.90') 
      GOTO 198 
C 
  568 WRITE(IOUERR,768) RVALUE 
  768 FORMAT(1X,/,12X,'Tier2 PM Phase-in Fraction of ',F7.3,' is less ', 
     &       'than zero.',/,12X,'PM Phase-in fraction must be between ', 
     &       'zero and one.',/,12X, 
     &       'Default Tier2 PM Phase-in Fractions will be used.' ) 
      GO TO 198 
C 
  569 WRITE(IOUERR,769) RVALUE 
  769 FORMAT(1X,/,12X,'Tier2 PM phase-in fraction of ',F7.3,' is ', 
     &       'greater than one.',/,12X,'PM Phase-in fraction must be ', 
     &       'between zero and one.',/,12X, 
     &       'Default Tier2 PM Phase-in Fractions will be used.' ) 
      GO TO 199 
C 
  570 WRITE(IOUERR,770) RVALUE 
  770 FORMAT(1X,/,12X,'Sum of 12 Tier2 PM phase-in bin fractions is ', 
     &       'equal to ',F7.3,/,12X,'Sum of 12 Tier2 PM bin fractions ', 
     &       'must be equal to 1.000 +/- 0.001',/,12X, 
     &       'Default Tier2 PM Phase-in Fractions will be used.' ) 
      GO TO 199 
C 
  571 WRITE(IOUERR,771) 
  771 FORMAT(1X/12X,'The Air Toxic OXYGENATE name must be MTBE,', 
     &       ' ETBE, ETOH (Ethanol), or TAME.') 
      GO TO 198 
C 
  572 WRITE(IOUERR,772) RVALUE 
  772 FORMAT(1X/12X,'The Air Toxic OXYGENATE volume percent ',F6.3, 
     &       ' is out of bounds. It must be between 0.0 and 20.0') 
      GO TO 198 
C 
  573 WRITE(IOUERR,773) RVALUE 
  773 FORMAT(1X/12X,'The Air Toxic OXYGENATE market share ',F6.3, 
     &       ' is out of bounds. It must be between 0.0 and 1.0') 
      GO TO 198 
C 
  574 WRITE(IOUERR,774) RVALUE 
  774 FORMAT(1X/12X,'The Air Toxic OXYGENATE market share total ',F6.3, 
     &       ' is out of bounds. It must be between 0.0 and 1.0') 
      GO TO 198 
C 
  575 WRITE(IOUERR,775) IVALUE 
  775 FORMAT(1X/12X,'The Air Toxic RVP OXY WAIVER flag ',I3, 
     &       ' is out of bounds. It must be 0 or 1.') 
      GO TO 198 
C 
  576 WRITE(IOUERR,776) INT(RVALUE),NAMVEH(IVALUE) 
  776 FORMAT(1X/12X,'The Air Toxic Technology Group Fraction Sum', 
     &       ' for MY ',I4,' ',A8,' is equal to zero.') 
      GO TO 199 
C 
  577 WRITE(IOUERR,777) 
  777 FORMAT(12X,'The entered values for Freeway VMT and Arterial' 
     *     /12X,'VMT fractions sum to zero. This is inconsistent with' 
     *     /12X,'the "AREAWIDE" option in the "AVERAGE SPEED" command.') 
      GO TO 198 
C 
  578 WRITE(IOUERR,778) 
  778 FORMAT(12X,'The entered values for Freeway VMT and Freeway' 
     *     /12X,'Ramp VMT fractions sum to zero. This is inconsistent' 
     *     /12X,'with the "FREEWAY" option in the "AVERAGE SPEED"' 
     *     /12X,'command.') 
      GO TO 198 
C 
  579 WRITE(IOUERR,779) 
  779 FORMAT(12X,'The combination of facility VMT values input in the' 
     *     /12X, 'in the "AVERAGE SPEED" command lead to singularities' 
     *     /12X, 'in the average speed calculation. Please check the' 
     *     /12X, 'the "FREEWAY" option in the "AVERAGE SPEED" command.' 
     *     /12X, 'The values have been reset to the MOBILE6 default' 
     *     /12X, 'values.') 
      GO TO 198 
C 
  580 WRITE(IOUERR,780) 
  780 FORMAT(12X,'The entered value for non-ramp Freeway VMT must' 
     *     /12X,'be greater than zero. A Freeway VMT value of zero is' 
     *     /12X,'inconsistent with the "FREEWAY" option of the ' 
     *     /12X,'"AVERAGE SPEED" command.') 
      GO TO 198 
C 
  581 WRITE(IOUREP,781) RVALUE 
  781 FORMAT(12X,'The user supplied freeway average speed of ',F4.1, 
     &      /12X,'will be used for all hours of the day. 100% of VMT' 
     &      /12X,'has been assigned to the freeway roadway type for' 
     &      /12X,'all hours of the day and all vehicle types.') 
      GO TO 199 
C 
  582 WRITE(IOUREP,782) RVALUE 
  782 FORMAT(12X,'The user supplied freeway average speed of ',F4.1, 
     &     /12X,'will be used for all hours of the day. 100% of VMT' 
     &     /12X,'has been assigned to a fixed combination of freeways' 
     &     /12X,'and freeway ramps for all hours of the day and all' 
     &     /12X,' vehicle types.') 
      GO TO 199 
C 
  583 WRITE(IOUREP,783) RVALUE 
  783 FORMAT(12X,'The user supplied arterial average speed of ',F4.1, 
     &     /12X,'will be used for all hours of the day.  100% of VMT' 
     &     /12X,'has been assigned to the arterial/collector roadway' 
     &    /12X,'type for all hours of the day and all vehicle types.') 
      GO TO 199 
C 
  584 WRITE(IOUREP,784) RVALUE 
  784 FORMAT(12X,'The user supplied area wide average speed of ',F4.1, 
     &     /12X,'will be used for all hours of the day.  100% of VMT' 
     &    /12X,'has been assigned to a fixed combination of freeways,' 
     &     /12X,'freeway ramps, arterial/collector and local roadways' 
     &     /12X,'for all hours of the day and all vehicle types.') 
      GO TO 199 
C 
  585 WRITE(IOUREP,785) RVALUE 
  785 FORMAT(12X,'100% of VMT has been assigned to the local roadway' 
     &    /12X,'type for all hours of the day for all vehicle types' 
     &    /12X,'with an average speed of ',F4.1,' mph.') 
      GO TO 199 
C 
  586 WRITE(IOUREP,786) RVALUE 
  786 FORMAT(12X,'100% of VMT has been assigned to the freeway ramp' 
     &      /12X,'roadway type for all hours of the day for all ' 
     &    /12X,'vehicle types with an average speed of ',F4.1,' mph.') 
      GO TO 199 
C 
  587 CONTINUE 
C 
  588 CONTINUE 
C 
  589 CONTINUE 
C 
  590 CONTINUE 
C 
  591 CONTINUE 
C 
  592 CONTINUE 
C 
  593 CONTINUE 
C 
  594 CONTINUE 
C 
  595 CONTINUE 
C 
  596 CONTINUE 
C 
  597 CONTINUE 
C 
  598 CONTINUE 
C 
  599 CONTINUE 
C 
  600 CONTINUE 
C 
      GO TO 197 
C 
  601 WRITE(IOUREP,801) 
  801 FORMAT(15X,'User has enabled STAGE II REFUELING.',/) 
      GO TO 199 
C 
  602 WRITE(IOUREP,802) 
  802 FORMAT(15X, 
     &    'The EXPAND EVAPORATIVE command was ignored because the'/15X, 
     &    'NO EVAPORATIVE command was issued in the Header Section.',/) 
      GO TO 199 
C 
  603 WRITE(IOUREP,803) 
  803 FORMAT(15X,'User has disabled the calculation ', 
     &                'of REFUELING emissions.',/) 
      GO TO 199 
C 
  604 WRITE(IOUREP,804) 
  804 FORMAT(15X,'User has disabled the COLD CO STANDARD', 
     &                ' corrections.',/) 
      GO TO 199 
C 
  605 CONTINUE 
C 
  606 GO TO 199 
C 
  607 GO TO 199 
C 
  608 WRITE(IOUREP,808) 
  808 FORMAT(15X,'User has disabled the RVP correction.',/) 
      GO TO 199 
C 
  609 WRITE(IOUREP,809) 
  809 FORMAT(15X,'User has disabled the SOAK correction.',/) 
      GO TO 199 
C 
  610 WRITE(IOUREP,810) 
  810 FORMAT(15X,'User has disabled the TEMPERATURE correction.',/) 
      GO TO 199 
C 
  611 WRITE(IOUREP,811) 
  811 FORMAT(15X,'User has disabled the SPEED correction.',/) 
      GO TO 199 
C 
  612 WRITE(IOUREP,812) 
  812 FORMAT(15X,'User has disabled the SFTP speed correction.',/) 
      GO TO 199 
C 
  613 WRITE(IOUREP,813) 
  813 FORMAT(15X,'User has disabled the SFTP air conditioning ', 
     &                   'correction.',/) 
      GO TO 199 
C 
  614 WRITE(IOUREP,814) 
  814 FORMAT(15X,'User supplied diesel sale fractions.') 
      GO TO 199 
C 
  615 WRITE(IOUREP,815) 
  815 FORMAT(15X,'User supplied VMT mix.') 
      GO TO 199 
C 
  616 WRITE(IOUREP,816) 
  816 FORMAT(15X,'User has supplied post-1999 sulfur levels.') 
      GO TO 199 
C 
  617 WRITE(IOUREP,817) RVALUE 
  817 FORMAT(15X,'User supplied alternate AC input: Cloud Cover ', 
     &           'Fraction set to',F5.2,'.') 
      GO TO 199 
C 
  618 WRITE(IOUREP,818) INT(RVALUE),IVALUE 
  818 FORMAT(15X,'User supplied alternate AC input: Sunrise at ',I2, 
     &           ' AM, Sunset at ',I2,' PM.') 
      GO TO 199 
C 
  619 WRITE(IOUREP,819) INT(RVALUE),IVALUE 
  819 FORMAT(15X,'User supplied alternate AC input: Peak Sun between ' 
     &          ,I2,' AM, and ',I2,' PM.') 
      GO TO 199 
C 
  620 WRITE(IOUREP,820) 
  820 FORMAT(15X,'User has disabled the air conditioning ', 
     &               'correction calculation.') 
      GO TO 199 
C 
  621 WRITE(IOUREP,821) 
  821 FORMAT(15X,' Results Reflect NO 2007 Heavy Duty Vehicle' 
     &            ' Emissions Reductions') 
      GOTO 199 
C 
  622 WRITE(IOUREP,822) 
  822 FORMAT(15X,'User has enabled the PM SPEED correction.',/) 
      GO TO 199 
C 
  623 CONTINUE 
C 
  624 CONTINUE 
C 
  625 CONTINUE 
C 
  626 CONTINUE 
C 
  627 CONTINUE 
C 
  628 CONTINUE 
C 
  629 CONTINUE 
C 
  630 CONTINUE 
C 
  631 CONTINUE 
C 
  632 CONTINUE 
C 
  633 CONTINUE 
C 
  634 CONTINUE 
C 
  635 CONTINUE 
C 
  636 CONTINUE 
C 
  637 CONTINUE 
C 
  638 CONTINUE 
C 
  639 CONTINUE 
C 
  640 CONTINUE 
C 
  641 CONTINUE 
C 
  642 CONTINUE 
C 
  643 CONTINUE 
C 
  644 CONTINUE 
C 
  645 CONTINUE 
C 
  646 CONTINUE 
C 
  647 CONTINUE 
C 
  648 CONTINUE 
C 
  649 CONTINUE 
C 
  650 CONTINUE 
C 
  651 CONTINUE 
C 
  652 CONTINUE 
C 
  653 CONTINUE 
C 
  654 CONTINUE 
C 
  655 CONTINUE 
C 
  656 CONTINUE 
C 
  657 CONTINUE 
C 
  658 CONTINUE 
C 
  659 CONTINUE 
C 
  660 CONTINUE 
C 
  661 CONTINUE 
C 
  662 CONTINUE 
C 
  663 CONTINUE 
C 
  664 CONTINUE 
C 
  665 CONTINUE 
C 
  666 CONTINUE 
C 
  667 CONTINUE 
C 
  668 CONTINUE 
C 
  669 CONTINUE 
C 
  670 CONTINUE 
C 
  671 CONTINUE 
C 
  672 CONTINUE 
C 
  673 CONTINUE 
C 
  674 CONTINUE 
C 
  675 CONTINUE 
C 
  676 CONTINUE 
C 
  677 CONTINUE 
C 
  678 CONTINUE 
C 
  679 CONTINUE 
C 
  680 CONTINUE 
C 
  681 CONTINUE 
C 
  682 CONTINUE 
C 
  683 CONTINUE 
C 
  684 CONTINUE 
C 
  685 CONTINUE 
C 
  686 CONTINUE 
C 
  687 CONTINUE 
C 
  688 CONTINUE 
C 
  689 CONTINUE 
C 
  690 CONTINUE 
C 
      GOTO 197 
C 
  691 WRITE(IOUREP,891) 
  891 FORMAT(/15X,'User has requested OBD RESPONSE.') 
      GOTO 195 
C 
  692 WRITE(IOUREP,892) 
  892 FORMAT(/15X,'User has requested OBD I/M COVERAGE.') 
      GOTO 195 
C 
  693 WRITE(IOUREP,893) 
  893 FORMAT(/15X,'User has requested OBD WORKING MIL.') 
      GOTO 195 
C 
  694 WRITE(IOUREP,894) 
  894 FORMAT(/15X,'User has requested GLL REDUCTION.') 
      GOTO 195 
C 
  695 WRITE(IOUREP,895) 
  895 FORMAT(15X,'User has disabled the EVAPORATIVE EMISSIONS ', 
     *           'calculations.') 
      GO TO 195 
C 
  696 WRITE(IOUREP,896) 
  896 FORMAT(15X,'User has disabled the calculation ', 
     *           'of RESTING LOSS emissions.') 
      GO TO 195 
C 
  697 WRITE(IOUREP,897) 
  897 FORMAT(/15X,'User has requested NO EVP SAWTOOTH.') 
      GOTO 195 
C 
  698 WRITE(IOUREP,898) 
  898 FORMAT(/15X,'User has requested FIXED PP RATES.') 
      GOTO 195 
C 
  699 WRITE(IOUREP,899) 
  899 FORMAT(15X,'User has disabled the WEATHERED RVP correction.') 
      GOTO 195 
C 
C  Output warning here if user stumbles upon/uses a "hidden" feature 
C 
  195 WRITE(IOUREP,495) 
  495 FORMAT( 
     *15X,'This feature is for EPA diagnostic use ONLY.'/ 
     *15X,'There is no user application of this feature'/ 
     *15X,'which would be acceptable for submission to EPA.'/) 
      GOTO 199 
C 
  196 IOS = IVALUE 
      WRITE(IOUOUT,496) IOS 
  496 FORMAT( 
     &  '********************************************************' 
     & /'*** Dispatch error writing to the report file: ',I8, 
     & /'********************************************************'/) 
      GO TO 198 
C 
C  Unknown error message numbers land here 
C 
  197 WRITE(IOUERR,497) MESSAG 
  497 FORMAT(13X,' Message code ',I4,' is unknown.') 
C  Messages 60 and 61 produce excessive output if not stopped here.  The 
C  other 
C  cases stopping execution here are considered "fatal" errors:  The user must 
C  fix the error-generating input and rerun. 
C 
  198 IF(MCODE.EQ.53.OR.MCODE.EQ.60.OR.MCODE.EQ.61.OR. 
     *   MCODE.EQ.85.OR.MCODE.EQ.89.OR.MCODE.EQ.97.OR. 
     *   MCODE.EQ.102.OR.MCODE.EQ.107.OR.MCODE.EQ.113.OR. 
     *   MCODE.EQ.118.OR.MCODE.EQ.161.OR.MCODE.EQ.163.OR. 
     *   MCODE.EQ.164.OR.MCODE.EQ.171.OR.MCODE.EQ.174.OR. 
     *   MCODE.EQ.178.OR.MCODE.EQ.536) 
     *   STOP 
      INERR=INERR+1 
C 
  199 RETURN 
      END 
