C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE QUITER(RVALUE,IVALUE,MESSAG,INERR)
C
C  QUITER prints error and warning messages on unit IOUERR.
C
C
C  Called by CALTEM, CALUDI, CCEVRT, CHKPGM, CONSEC, CRANKC, DIURNL,
C            EFCALX, EXMPTR, GETALH, GETATP, GETBYM, GETCOL, GETEFR,
C            GETERQ, GETETP, GETFYR, GETIMC, GETIMP, GETIM2, GETMRN,
C            GETNAT, GETPRG, GETPRS, GETSC1, GETTAM, GETVER, GETVMT, 
C            HCCALX, IBFPTR, IERPTR, IEVPTR, IGEPTR, IGLPTR, IGSFPT,
C            IGTPTR, IGUPTR, IRULPT, ISPPTR, ITAMPT, ITCPTR, LOCAL,
C            MOBILE, REFUEL, RDCITY, STOERQ, TFCALX, WTLRUL, YRTEST.
C
C  Input on call:
C
C    parameter list: RVALUE,IVALUE,MESSAG,INERR
C    common blocks:
C    /ATPAR2/ ATPPGM,ATPFQT,CRATP
C    /BASEQ4/ NEWPAR,MAXCT
C    /BASEQ6/ MAXERU
C    /FLAGS2/ IMFLAG
C    /FLAGS3/ ATPFLG
C    /IMPAR1/ CRIM
C    /IMPAR6/ IFREQ,INTYP
C    /IOUCOM/ IOUREP,IOUERR,IOUASK
C    /NAMES1/ NAMFLG,NAMRVP,NAMTEM
C    /PRGCH2/ PRGPGM,PRGFQT,CRPRG
C    /PRSCH2/ PRSPGM,PRSFQT,CRPRS
C    /QUITXQ/ N1QUIT
C    /YEARS4/ IYSTART,IY1960,IYEND
C
C  Output on return:
C
C    parameter list: INERR
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
C  NAMTMP(8)    -  NAMTMP ( IHRD+2 )
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
C  NAMFNC   C*8  function names (to identify subprogram source of error)
C  NAMFRQ   C*8  program inspection frequency cases
C  NAMPOL   C*5  pollutant names
C  NAMTMP   C*13 temperature names
C  NAMTYP   C*8  program inspection type cases
C  NAMVEH   C*20 vehicle type names
C
C
C  Notes:
C
C  Message slots available as of 02/09/94: NONE
C  Some specific things about QUITER are:
C  All messages must fit in 80 columns
C  A blank needs to start at the beginning of a line so it does
C  not get cut off in MTS
C  Each line of a warning is preceded by 14X and errors 12X
C  Use available message slots before expanding MAXMES upper limit
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
C  7-May-1994 @ CSC-tmm Request #2-441 Modified code to permit the detergent gasoline
C     parameters to be read in and stored in COMMON. Added msgs. 69,70, & 85.
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
C  9 Feb 96 @DynTel-MLA 2-611 Added messages 102, 113, 157, 167, and 168.
C  13-Sep-1996 @DynTel-yc Request #2-630 Added messages 169 and 170 for
C               detergent gasoline effects.
C 
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
      INCLUDE 'QUITXQ.I'
      INCLUDE 'YEARS4.I'
C
      INTEGER QFREQ,QINTYP,QAFREQ,QFQTFN,QPGMFN
      INTEGER QAPROG
      DIMENSION QFREQ(2),QINTYP(2),QFQTFN(2),QPGMFN(2)
C
      CHARACTER*1 NAMEFF(2)
      CHARACTER*4 IOOFBD(5)
      CHARACTER*8 MERR,MWARN,MCOM,MTYPE
      CHARACTER*8 NAMVEH(8),NAMFNC(17)
      CHARACTER*8 NAMFRQ(3)
      CHARACTER*5 NAMPOL(3)
      CHARACTER*20 NAMTYP(2,4)
      CHARACTER*13 NAMTMP(8)
      CHARACTER*10 NAMDET(9)
C
      DATA NAMEFF/'L','H'/
C
      DATA IOOFBD/
     *   '  ou','t of',' bou','nds ','for '/
C
      DATA MERR/'Error:  '/,MWARN/'Warning:'/,MCOM/'Comment:'/
C
      DATA NAMVEH/
     *   'LDGV    ',
     *   'LDGT1   ',
     *   'LDGT2   ',
     *   'HDGV    ',
     *   'LDDV    ',
     *   'LDDT    ',
     *   'HDDV    ',
     *   'MC      '/
C
      DATA NAMFNC/
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
     *   ' IDLPTR '/
C
      DATA NAMFRQ/'Annual  ','Biennial','NA      '/
C
      DATA NAMPOL/' HC  ',' CO  ',' NOx '/
C
      DATA NAMTMP/'ambient      ',
     *            'exhaust      ',
     *            'hot soak     ',
     *            'running loss ',
     *            'diurnal 9 AM ',
     *            'diurnal 11 AM',
     *            'diurnal 1 PM ',
     *            'resting loss '/
C
      DATA NAMTYP/'Test Only           ','Test Only           ',
     *            'Comp Test & Repair  ','Test & Repair       ',
     *            'Manual Test & Repair','Test & Repair       ',
     *            'Not Applicable      ','Not Applicable      '/
      DATA NAMDET/
     *   'DGPHSF    ',
     *   'FIDFRC    ',
     *   'IVDFRC    ',
     *   'FIDEFF(1) ',
     *   'FIDEFF(2) ',
     *   'FIDEFF(3) ',
     *   'IVDEFF(1) ',
     *   'IVDEFF(2) ',
     *   'IVDEFF(3) '/
C
C  Aug-18-1994 @ CSC-bsg Request 2-462 Added QUITER messages 157-160, for 
C  non-FTP related error messages
C
      DATA MAXMES/170/
C
      MCODE=MESSAG
C
C  Trap unknown errors.
C
      IF(MCODE.LT.1.OR.MCODE.GT.MAXMES) MCODE=0
C
C  Double space if this is first call to QUITER on current run.
C
      IF(N1QUIT.EQ.1) WRITE(IOUERR,298)
  298 FORMAT('0',' ')
      IF(N1QUIT.EQ.1) N1QUIT=0
C
C  Identify type of message.
C
      MTYPE=MERR
      IF(MCODE.EQ.21.OR.MCODE.EQ.22
     *   .OR.MCODE.EQ.31.OR.MCODE.EQ.49.OR.MCODE.EQ.52.OR.MCODE.EQ.58
     *   .OR.MCODE.EQ.59.OR.MCODE.EQ.90.OR.MCODE.EQ.96.OR.MCODE.EQ.98
     *   .OR.(MCODE.GE.108.AND.MCODE.LE.110)
     *   .OR.MCODE.EQ.112.OR.MCODE.EQ.114
     *   .OR.(MCODE.GE.119.AND.MCODE.LE.121).OR.MCODE.EQ.135
     *   .OR.MCODE.EQ.146.OR.(MCODE.GE.152.AND.MCODE.LE.156)
     *   .OR.MCODE.EQ.158.OR.MCODE.EQ.165.OR.MCODE.EQ.166
     *   .OR.MCODE.EQ.170)
     *   MTYPE=MWARN
      IF(MCODE.EQ.53.OR.MCODE.EQ.56.OR.MCODE.EQ.71.OR.MCODE.EQ.83)
     *   MTYPE=MCOM
C
C  Write message type prefix:
C
      WRITE(IOUERR,299) MCODE,MTYPE
  299 FORMAT('-',T2,'M',I3,1X,A8)
C
C  Branch to the appropriate error / warning / comment message.
C
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
     *    161, 162, 163, 164, 165, 166, 167, 168, 169, 170) MCODE
      GOTO 197
C
C  Out of bounds errors in flags land here (flag # is in RVALUE)
C
    1 MCODE=IFIX(RVALUE)
      WRITE(IOUERR,301) IVALUE,IOOFBD,NAMFLG(MCODE)
  301 FORMAT('+',12X,I4,5A4,'flag',A8)
      GOTO 198
    2 WRITE(IOUERR,302) RVALUE,IOOFBD
  302 FORMAT('+',12X,G10.3,5A4,'Ether Blend Market Share (0. to 1.)')
      GOTO 198
    3 WRITE(IOUERR,303) RVALUE,IOOFBD
  303 FORMAT('+',12X,G10.3,5A4,'Alcohol Blend Market Share (0. to 1.)')
      GOTO 198
    4 WRITE(IOUERR,304) RVALUE,IOOFBD
  304 FORMAT('+',12X,G10.3,5A4,'Total Oxy Fuel Market Share (0. to 1.)')
      GOTO 198
    5 WRITE(IOUERR,305) RVALUE,IOOFBD
  305 FORMAT('+',12X,G10.3,5A4,'Ether Blend Oxygen Content',
     *    ' (0 to 0.027)')
      GOTO 198
    6 WRITE(IOUERR,306) RVALUE,IOOFBD
  306 FORMAT('+',12X,G10.3,5A4,'Alcohol Blend Oxygen Content',
     *    ' (0 to 0.035)')
      GOTO 198
    7 WRITE(IOUERR,307) IVALUE,IOOFBD
  307 FORMAT('+',12X,I4,5A4,'Alcohol Blend RVP Waiver switch (1 or 2)')
      GOTO 198
    8 WRITE(IOUERR,308)
  308 FORMAT('+',12X,'Inconsistent Ether Blend input encountered')
      GOTO 198
    9 WRITE(IOUERR,309)
  309 FORMAT('+',12X,'Inconsistent Alcohol Blend input encountered')
      GOTO 198
   10 WRITE(IOUERR,310) IVALUE
  310 FORMAT('+',12X,I1,
     *' is not an allowable value for the Maximum I/M credit flag.')
      GOTO 198
   11 WRITE(IOUERR,311) IVALUE
  311 FORMAT('+',12X,I4,
     *' for LEV start year is less than 1994.')
      GOTO 198
   12 WRITE(IOUERR,312) IVALUE,IOOFBD
  312 FORMAT('+',12X,I1,5A4,'number of cold CO standards (1 to 2)')
      GOTO 198
   13 WRITE(IOUERR,313) IVALUE,IOOFBD,IYEND
  313 FORMAT('+',12X,I4,5A4,'cold CO standard start MY (1980-',I4,')')
      GOTO 198
   14 WRITE(IOUERR,314) RVALUE,IOOFBD,IVALUE
  314 FORMAT('+',12X,F5.1,5A4,'cold CO standard #',I1,' (2.0-15.0)')
      GOTO 198
   15 WRITE(IOUERR,315)
  315 FORMAT('+',12X,'cold CO standard #1 must start before #2')
      GOTO 198
   16 WRITE(IOUERR,316) RVALUE,IOOFBD,IVALUE
  316 FORMAT('+',12X,F8.3,5A4,'% IV under CCO Std #',I1,' (0.0-100.0)')
      GOTO 198
   17 WRITE(IOUERR,317) RVALUE,IOOFBD
  317 FORMAT('+',12X,G10.3,5A4,'VMTMIX (0. to 1.)')
      GOTO 198
   18 WRITE(IOUERR,318) RVALUE
  318 FORMAT('+',12X,G10.3,' sum of VMTMIX is not equal to 1.')
      GOTO 198
   19 WRITE(IOUERR,319) RVALUE
  319 FORMAT('+',12X,G10.3,' negative model year mileage')
      GOTO 198
   20 WRITE(IOUERR,320) RVALUE
  320 FORMAT('+',12X,G10.3,' negative model year registration')
      GOTO 198
   21 WRITE(IOUERR,321) RVALUE
  321 FORMAT('+',14X,G10.3,' registration with zero mileage')
      GOTO 199
   22 WRITE(IOUERR,322) RVALUE
  322 FORMAT('+',14X,G10.3,' mileage with zero registration')
      GOTO 199
   23 WRITE(IOUERR,323) IVALUE,IOOFBD,MAXCT
  323 FORMAT('+',12X,I4,5A4,'e.f. changes (1 to',I4,' pairs)')
      GOTO 198
   24 WRITE(IOUERR,324) IVALUE,IOOFBD
  324 FORMAT('+',12X,I4,5A4,'region (1 to 2)')
      GOTO 198
   25 WRITE(IOUERR,325) IVALUE,IOOFBD
  325 FORMAT('+',12X,I4,5A4,'vehicle type (1 to 8)')
      GOTO 198
   26 WRITE(IOUERR,326) IVALUE,IOOFBD
  326 FORMAT('+',12X,I4,5A4,'pollutant (1 to 3)')
      GOTO 198
   27 WRITE(IOUERR,327) IVALUE,IOOFBD,IYSTART,IYEND
  327 FORMAT('+',12X,I4,5A4,'year (',I4,' to ',I4,')')
      GOTO 198
   28 WRITE(IOUERR,328)
  328 FORMAT('0',T2,'Excess data errors prevent further analysis.')
      GOTO 199
   29 WRITE(IOUERR,329)
  329 FORMAT('+',12X,' last year cannot be less than first year')
      GOTO 198
   30 WRITE(IOUERR,330) RVALUE
  330 FORMAT('+',12X,G10.3,' intercept must be positive')
      GOTO 198
   31 WRITE(IOUERR,331) RVALUE
  331 FORMAT('+',14X,G10.3,' negative slope for ageing vehicle')
      GOTO 199
   32 WRITE(IOUERR,327) IVALUE,IOOFBD,IY1960,IYEND
      GOTO 198
   33 WRITE(IOUERR,333) IVALUE,IOOFBD
  333 FORMAT('+',12X,I4,5A4,'stringency (10 to 50)')
      GOTO 198
   34 WRITE(IOUERR,334) RVALUE,IOOFBD
  334 FORMAT('+',12X,F4.0,5A4,'I/M Program compliance rate',
     *       ' ',' (0 to 100%)')
      GOTO 198
   35 WRITE(IOUERR,335) IVALUE,IOOFBD,IYSTART,IYEND
  335 FORMAT('+',12X,I4,5A4,'MODYR1 (',I4,' to ',I4,')')
      GOTO 198
   36 WRITE(IOUERR,336) IVALUE,IOOFBD,IYSTART,IYEND
  336 FORMAT('+',12X,I4,5A4,'MODYR2 (',I4,' to ',I4,')')
      GOTO 198
   37 WRITE(IOUERR,337)
  337 FORMAT('+',12X,' MODYR1 cannot be greater than MODYR2')
      GOTO 198
   38 WRITE(IOUERR,327) IVALUE,IOOFBD,IY1960,IYEND
      GOTO 198
   39 WRITE(IOUERR,339) RVALUE
  339 FORMAT('+',12X,G10.3,' speed must be positive')
      GOTO 198
   40 WRITE(IOUERR,340) RVALUE
  340 FORMAT('+',12X,G10.3,' valid ambient temperature',
     *                    ' is 0-110 deg. (F)')
      GOTO 198
   41 WRITE(IOUERR,341) RVALUE,IOOFBD
  341 FORMAT('+',12X,G10.3,5A4,'PCCN (0. to 100.)')
      GOTO 198
   42 WRITE(IOUERR,342) RVALUE,IOOFBD
  342 FORMAT('+',12X,G10.3,5A4,'PCHC (0. to 100.)')
      GOTO 198
   43 WRITE(IOUERR,343) RVALUE,IOOFBD
  343 FORMAT('+',12X,G10.3,5A4,'PCCC (0. to 100.)')
      GOTO 198
   44 WRITE(IOUERR,344) RVALUE,IOOFBD
  344 FORMAT('+',12X,G10.3,5A4,'AC (0. to 1.)')
      GOTO 198
   45 WRITE(IOUERR,345) RVALUE,IOOFBD
  345 FORMAT('+',12X,G10.3,5A4,'extra load (0. to 1.)')
      GOTO 198
   46 WRITE(IOUERR,346) RVALUE,IOOFBD
  346 FORMAT('+',12X,G10.3,5A4,'trailers (0. to 1.)')
      GOTO 198
   47 WRITE(IOUERR,347) RVALUE,IOOFBD
  347 FORMAT('+',12X,G10.3,5A4,'humidity (20. to 140.)')
      GOTO 198
   48 WRITE(IOUERR,348) NAMVEH(IVALUE)
  348 FORMAT('+',12X,' there are no sales for vehicle ',A8)
      GOTO 198
   49 WRITE(IOUERR,349) RVALUE
  349 FORMAT('+',14X,G10.3,' MYR sum not = 1. (will normalize)')
      GOTO 199
   50 WRITE(IOUERR,350) RVALUE,IOOFBD
  350 FORMAT('+',12X,G10.3,5A4,'tampering rate intercept (up to 1.0)')
      GOTO 198
   51 WRITE(IOUERR,351) IVALUE,MAXERU
  351 FORMAT('+',12X,' equation',I3,' can not be stored.  Only',I3,
     *  ' equations'/' ',12X,' per pol/veh/reg case are allowed.'/
     *  '0',' ')
      GOTO 198
   52 WRITE(IOUERR,352) RVALUE
  352 FORMAT('+',14X,G10.3,' speed increased to 2.5 mph minimum')
      GOTO 199
   53 WRITE(IOUERR,353)
  353 FORMAT('+',14X,' All flags must be corrected before rerunning.')
      GOTO 198
   54 WRITE(IOUERR,354)
  354 FORMAT('-','I/M credits for this program description ',
     *       'not found.'/'0')
      GOTO 198
   55 WRITE(IOUERR,355)
  355 FORMAT('-','Error return reading I/M credits.'/'0')
      GOTO 198
   56 WRITE(IOUERR,356)
  356 FORMAT('+',14X,' A/C correction factor will be calculated.'
     *   /' ',T15,' Value of inputted AC usage parameter is ignored.')
      GOTO 199
   57 WRITE(IOUERR,357)
  357 FORMAT('+',12X,'WB temp cannot be greater than DB temp')
      GOTO 198
   58 WRITE(IOUERR,358) IVALUE,RVALUE
  358 FORMAT('+',14X,' Average miles/day is zero for vehicle class',I2,
     *               ' and model year',F6.0)
      GOTO 199
   59 WRITE(IOUERR,359) IVALUE,(NEWPAR(IPAR,IVALUE),IPAR=1,3)
  359 FORMAT('+',14X,' equation',I3,' zeroes all idle coefficients',
     *   ' (and total) for IR=',I1,', IV=',I1,', IP=',I1)
      GOTO 199
   60 WRITE(IOUERR,360) INERR,NAMFNC(IVALUE)
  360 FORMAT('+',12X,I4,' exceeds ub of highest myg in',A8)
      GOTO 198
   61 WRITE(IOUERR,361) INERR,NAMFNC(IVALUE)
  361 FORMAT('+',12X,' default used for',I5,' in index function',A8)
      GOTO 198
   62 WRITE(IOUERR,362) IVALUE,IOOFBD
  362 FORMAT('+',12X,I4,5A4,'new e.f. region (1 to 2)')
      GOTO 198
   63 WRITE(IOUERR,363) RVALUE,IOOFBD
  363 FORMAT('+',12X,G10.3,5A4,'(0. <= PCHC+PCCC <= 100.)')
      GOTO 198
   64 WRITE(IOUERR,364) RVALUE,IOOFBD
  364 FORMAT('+',12X,G10.3,5A4,'(0. <= PCHC+PCCC-PCCN <= 100.)')
      GOTO 198
   65 WRITE(IOUERR,365) IVALUE
  365 FORMAT('+',12X,' MYR of LDDV not equal to LDGV for JDX=',I2)
      GOTO 198
   66 WRITE(IOUERR,366) IVALUE
  366 FORMAT('+',12X,' MYR of LDDT not equal to LDGT1 for JDX=',I2)
      GOTO 198
   67 WRITE(IOUERR,367) NAMVEH(IVALUE)
  367 FORMAT('+',12X,' EFFTP>0. and GSF=0. for vehicle ',A8)
      GOTO 198
   68 WRITE(IOUERR,368) NAMVEH(IVALUE)
  368 FORMAT('+',12X,' EFFTP<=0. and VMTMIX>0. for vehicle ',A8)
      GOTO 198
C
   69 WRITE(IOUERR,369)
  369 FORMAT(
     *'+',12X,' Error occurred in reading the',
     *        ' detergent gasoline parameter record.')
      GOTO 198
C
   70 WRITE(IOUERR,370) RVALUE*100, IOOFBD, NAMDET(IVALUE)
  370 FORMAT(
     *'+',12X,' Error occurred in reading the',
     *        ' detergent gasoline parameter record.'/
     *' ',12X,F6.2,5A4,A10)
      GOTO 198
   71 WRITE(IOUERR,371) IOUREP,IOUERR,IOUASK
  371 FORMAT('+',14X,' Current output unit numbers are',
     *   ' IOUREP=',I1,' IOUERR=',I1,' IOUASK=',I1)
      GOTO 199
   72 WRITE(IOUERR,372) IVALUE,IOOFBD
  372 FORMAT('+',12X,I4,5A4,'I/M vehicle type (1 to 2)')
      GOTO 198
   73 WRITE(IOUERR,373) IVALUE,IOOFBD
  373 FORMAT('+',12X,I4,5A4,'short test type flag (1 to 3)')
      GOTO 198
   74 MCODE=IFIX(RVALUE)
      WRITE(IOUERR,374) IVALUE,IOOFBD,NAMFLG(MCODE)
  374 FORMAT('+',12X,I4,5A4,'flag',A8/
     *' ',12X,'this value can only be used with .CSV output files')
      GOTO 198
   75 WRITE(IOUERR,375) IVALUE,IOOFBD
  375 FORMAT('+',12X,I4,5A4,'anti-tampering program start year')
      GOTO 198
   76 WRITE(IOUERR,376) IVALUE,IOOFBD
  376 FORMAT('+',12X,I4,5A4,'anti-tampering program 1st model year')
      GOTO 198
   77 WRITE(IOUERR,377) IVALUE,IOOFBD
  377 FORMAT('+',12X,I4,5A4,'anti-tampering program last model year')
      GOTO 198
   78 WRITE(IOUERR,378) IVALUE,IOOFBD
  378 FORMAT('+',12X,I4,5A4,'anti-tampering vehicle type (1 to 2)')
      GOTO 198
   79 WRITE(IOUERR,379) RVALUE,IOOFBD
  379 FORMAT('+',12X,G10.3,5A4,'evap ATP effectiveness rate (0. to 1.)')
      GOTO 198
   80 WRITE(IOUERR,380) RVALUE,IOOFBD
  380 FORMAT('+',12X,G10.3,5A4,'exh ATP effectiveness rate (0. to 1.)')
      GOTO 198
   81 WRITE(IOUERR,381) RVALUE,IOOFBD
  381 FORMAT('+',12X,G10.3,5A4,'AER matrix file column sum (0. to 1.)')
      GOTO 198
   82 WRITE(IOUERR,382) IVALUE,IOOFBD
  382 FORMAT('+',12X,I4,5A4,'type of I/M inspection (1 (Test only ',
     *'or 2 (Test and Repair (Computerized)) or'/12X,
     *'3 (Test and Repair (Manual)))')
      GOTO 198
   83 WRITE(IOUERR,383)
  383 FORMAT(
     *'+',14X,'One or more evaporative temperatures (input daily'/
     *' ',14X,'maximum, input ambient, calculated hot soak, and/or'/
     *' ',14X,'calculated running loss) is 40F or less, or input'/
     *' ',14X,'daily minimum is 25F or less; no evaporative emission'/
     *' ',14X,'factors (hot soak, diurnal, running loss, or resting'/
     *' ',14X,'loss) will be calculated.')
      GOTO 199
   84 WRITE(IOUERR,384) IVALUE,IOOFBD
  384 FORMAT('+',12X,I4,5A4,'onboard VRS vehicle class flag (1 to 2)')
      GOTO 198
   85 WRITE(IOUERR,385) RVALUE,NAMPOL(IVALUE)
  385 FORMAT('+',12X,G10.3,' is out of bounds for the combined ',
     *        'effects of Detergent' /
     *        14X,'Gasoline on',A5,'as entered by the user.')
      GOTO 198
   86 WRITE(IOUERR,386) IVALUE,IOOFBD
  386 FORMAT('+',12X,I4,5A4,'no. of Stage II phase-in years (1 to 9)')
      GOTO 198
   87 WRITE(IOUERR,387) IVALUE,IOOFBD
  387 FORMAT('+',12X,I4,5A4,'vrs or rvp controls start year')
      GOTO 198
   88 WRITE(IOUERR,388) RVALUE,IOOFBD
  388 FORMAT('+',12X,G10.3,5A4,'diesel sales fraction (0. to 1.)')
      GOTO 198
   89 WRITE(IOUERR,389) IVALUE,IOOFBD
  389 FORMAT('+',12X,I4,5A4,' flag PROMPT (1 to 6)')
      GOTO 198
   90 WRITE(IOUERR,390) NAMRVP(IVALUE),RVALUE
  390 FORMAT('+',14X,A8,' RVP reset to ',F4.1)
      GOTO 199
   91 WRITE(IOUERR,391) RVALUE,IOOFBD,NAMTEM(IVALUE)
  391 FORMAT('+',12X,F10.1,5A4,A4,'mum daily temperature')
      GOTO 198
   92 WRITE(IOUERR,392) RVALUE
  392 FORMAT('+',12X,' max daily temp =',F4.0,' < min daily temp')
      GOTO 198
   93 WRITE(IOUERR,393) RVALUE,IOOFBD
  393 FORMAT('+',12X,F10.2,5A4,'%TLVMT (set to zero.)')
      GOTO 198
   94 WRITE(IOUERR,394) RVALUE
  394 FORMAT('+',14X,G10.3,' %TLVMT sum not = 100. (will normalize)')
      GOTO 198
   95 WRITE(IOUERR,395) RVALUE
  395 FORMAT('+',12X,' ambient temp =',F4.0,' is < daily min temp',
     *               ' or > daily max temp')
      GOTO 198
   96 WRITE(IOUERR,396) RVALUE
  396 FORMAT('+',14X,G10.3,' speed reduced to 65 mph maximum')
      GOTO 199
   97 WRITE(IOUERR,397) RVALUE,IOOFBD,IVALUE
  397 FORMAT('+',12X,F10.3,5A4,'Wade Index calc (CALUDI) of A, pass ',
     *                      I1,': must be >= 0.0')
      GOTO 198
   98 WRITE(IOUERR,398) RVALUE
  398 FORMAT(
     *'+',14X,'Diurnal temperature rise (max temp-min temp=',F5.1,')'/
     *' ',14X,'is > 40F; diurnal evaporative emission factors will'/
     *' ',14X,'be calculated, but may be inaccurate.')
      GOTO 199
   99 WRITE(IOUERR,399) RVALUE,IOOFBD,NAMEFF(IVALUE)
  399 FORMAT('+',12X,F5.0,5A4,'Stage II efficiency for ',A1,'DGV ',
     *          '(0 to 100 %)')
      GOTO 198
  100 WRITE(IOUERR,400) RVALUE,IOOFBD
  400 FORMAT('+',12X,F3.0,5A4,'I/M Program waiver rate (0 to 50%)')
      GOTO 198
  101 WRITE(IOUERR,401) IVALUE,IOOFBD
  401 FORMAT('+',12X,I4,5A4,'frequency of I/M inspection (1 (annual) ',
     *                      'to 2 (biennial))')
      GOTO 198
  102 WRITE(IOUERR,402)
  402 FORMAT('+',12X,'Invalid character in the effectiveness flag',
     &  ' field of the I/M program', /' ',12X,'descriptive record.',
     &  ' The field should contain 1, 2, or a blank.')
      GOTO 198
  103 WRITE(IOUERR,403) IVALUE,IOOFBD
  403 FORMAT('+',12X,I4,5A4,'type of ATP program (1 (Test only)',
     *'to 2 (Test and Repair))')
      GOTO 198
  104 WRITE(IOUERR,404) IVALUE,IOOFBD
  404 FORMAT('+',12X,I4,5A4,'frequency of ATP inspection (1 (annual) ',
     *                      'to 2 (biennial))')
      GOTO 198
  105 WRITE(IOUERR,405) RVALUE,IOOFBD
  405 FORMAT('+',12X,F4.0,5A4,'ATP compliance rate must be'/
     *       ' ',36X,'in the range 0.% to 100.%, inclusive')
      GOTO 198
  106 WRITE(IOUERR,406) IVALUE,IOOFBD
  406 FORMAT('+',12X,I4,5A4,'ATP disablement inspections (1 (no) ',
     *                      'to 2 (yes))')
      GOTO 198
  107 WRITE(IOUERR,407)
  407 FORMAT('+',12X, ' Error reading ATP program description',
     *                ' in the ONE-TIME data section')
      GOTO 198
  108 WRITE(IOUERR,408)
      IF(ATPFLG.EQ.2) WRITE(IOUERR,908) CRATP
      IF(IMFLAG.EQ.2) WRITE(IOUERR,909) CRIM(1)
      IF(IMFLAG.EQ.3) WRITE(IOUERR,910) CRIM(1),CRIM(2)
      IF(PRSFLG.EQ.2) WRITE(IOUERR,911) CRPRS
      IF(PRGFLG.EQ.2) WRITE(IOUERR,912) CRPRG
  408 FORMAT('+', 14X, 'The user supplied ',
     *                 'compliance rates are not equal.')
  908 FORMAT(' ', 14X, 'ATP compliance rate:             ',9X,F5.1)
  909 FORMAT(' ', 14X, 'I/M program #1 compliance rate:  ',9X,F5.1)
  910 FORMAT(' ', 14X, 'I/M program #1 compliance rate:  ',9X,F5.1/
     *       ' ', 14X, 'I/M program #2 compliance rate:  ',9X,F5.1)
  911 FORMAT(' ', 14X, 'Pressure Check compliance rate:  ',9X,F5.1)
  912 FORMAT(' ', 14X, 'Purge Check compliance rate:     ',9X,F5.1)
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
      WRITE(IOUERR,409) NAMFRQ(QAFREQ),(NAMFRQ(QFREQ(I)),I=1,2),
     *                 (NAMFRQ(QFQTFN(I)),I=1,2)
  409 FORMAT('+', 14X, 'The user supplied inspection ',
     *                 'frequencies are not equal.'/
     *       ' ', 14X, 'ATP inspection frequency:                 ',A8/
     *       ' ', 14X, 'I/M program #1 inspection frequency:      ',A8/
     *       ' ', 14X, 'I/M program #2 inspection frequency:      ',A8/
     *       ' ', 14X, 'Pressure Check inspection frequency:      ',A8/
     *       ' ', 14X, 'Purge Check inspection frequency:         ',A8)
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
      WRITE(IOUERR,410) NAMTYP(2,QAPROG),(NAMTYP(1,QINTYP(I)),I=1,2),
     *                 (NAMTYP(2,QPGMFN(I)),I=1,2)
  410 FORMAT('+',14X,'The user supplied inspection ',
     *               'types are not equal.'/
     *     ' ', 14X, 'ATP inspection type:               ',A20/
     *     ' ', 14X, 'I/M program #1 inspection type:    ',A20/
     *     ' ', 14X, 'I/M program #2 inspection type:    ',A20/
     *     ' ', 14X, 'Pressure Check inspection type:    ',A20/
     *     ' ', 14X, 'Purge Check inspection type:       ',A20)
CC
      GOTO 199
  111 WRITE(IOUERR,411) NAMTMP(IVALUE),RVALUE
  411 FORMAT('+',12X,' The calculated ',A13,'temperature',F6.1,
     *               ' is < daily min temp or > daily max temp')
      GOTO 198
CC
  112 WRITE(IOUERR,412)
  412 FORMAT('+',14X,'Wintertime Reformulated Gasoline Rules Apply'/
     *       ' ',14X,'See Section 2.2.16 of the User Guide')
          GOTO 199
C
  113 WRITE(IOUERR,413)
  413 FORMAT('+',14X,'There was an error reading values from the ',
     & 'program effectiveness',/' ',14X,'record. Please check format.')
      GOTO 198
  114 WRITE(IOUERR,414)
  414 FORMAT('+',14X,'Purge Check emission benefits assume the use ',
     *               'of a dynomometer and'/' ',14X,'the IM240 ',
     *               'transient test procedure driving cycle.')
      GOTO 199
  115 WRITE(IOUERR,415) IVALUE,IOOFBD
  415 FORMAT('+',12X,I2,5A4,'number of ETP phase-in years (1 to 4).')
      GOTO 198
  116 WRITE(IOUERR,416) RVALUE,IOOFBD
  416 FORMAT('+',12X,F4.0,5A4,'ETP phase-in rate (0 to 100%).')
      GOTO 198
  117 WRITE(IOUERR,417) IVALUE,IOOFBD
  417 FORMAT('+',12X,I4,5A4,'by model year table pick (1 to 2)')
      GOTO 198
  118 WRITE(IOUERR,418)
  418 FORMAT(
     *'+',14X,'Input of diesel sales fractions not allowed'/
     *' ',14X,'in One-time Data Section.')
      GOTO 198
  119 WRITE(IOUERR,419)
  419 FORMAT(
     *'+',14X,'Zero tampering (TAMFLG=0) is for EPA in-house use ONLY.')
      GOTO 195
  120 WRITE(IOUERR,420) IVALUE,IOOFBD
  420 FORMAT('+',14X,I4,5A4,'evaluation month of a reformulated gas'/
     *' ',14X,'area. Must be 1 (Jan) or 7 (July). Season will default',
     *' to summer.')
      GOTO 199
  121 WRITE(IOUERR,421)
  421 FORMAT(
     *'+',14X,'NO emission credit from evaporative control',
     *        ' system inspection'/
     *' ',14X,'will be given unless gas cap inspection also required.')
      GOTO 199
  122 WRITE(IOUERR,422)
  422 FORMAT(
     *'+',12X,' Error occurred in reading mileage',
     *        ' accrual rates.')
      GOTO 198
  123 WRITE(IOUERR,423)
  423 FORMAT(
     *'+',12X,' Error occurred in reading registration distributions.')
      GOTO 198
  124 WRITE(IOUERR,424)
  424 FORMAT(
     *'+',12X,' Error occurred in reading a vertically input flag.')
      GOTO 198
  125 WRITE(IOUERR,425)
  425 FORMAT(
     *'+',12X,' Error occurred in reading a horizontally input flag.')
      GOTO 198
  126 WRITE(IOUERR,426)
  426 FORMAT(
     *'+',12X,' Error occurred in reading a',
     *        ' replacement tampering rate.')
      GOTO 198
  127 WRITE(IOUERR,427)
  427 FORMAT(
     *'+',12X,' Error occurred in reading the VMTmix record.')
      GOTO 198
  128 WRITE(IOUERR,428)
  428 FORMAT(
     *'+',12X,' Error occurred in reading the',
     *        ' replacement emission factor rates.')
      GOTO 198
  129 WRITE(IOUERR,429)
  429 FORMAT(
     *'+',12X,' Error occurred in reading the',
     *        ' I/M program descriptive record.')
      GOTO 198
  130 WRITE(IOUERR,430)
  430 FORMAT(
     *'+',12X,' Error occurred in reading the',
     *        ' Stage II or Onboard VRS input record.')
      GOTO 198
  131 WRITE(IOUERR,431)
  431 FORMAT(
     *'+',12X,' Error occurred in reading the',
     *        ' local area parameter record.')
      GOTO 198
  132 WRITE(IOUERR,432)
  432 FORMAT(
     *'+',12X,' Error occurred in reading the',
     *        ' scenario descriptive record.')
      GOTO 198
  133 WRITE(IOUERR,433)
  433 FORMAT(
     *'+',12X,' Error occurred in reading the',
     *        ' additional correction factor record.')
      GOTO 198
  134 WRITE(IOUERR,434)
  434 FORMAT(
     *'+',12X,' Error occurred in reading the',
     *        ' BY-Model-Year selection flags.')
      GOTO 198
  135 WRITE(IOUERR,435)
  435 FORMAT(
     *'+',12X,' All effects of the 1990 Clean Air Act Amendments',
     *' have been disabled.')
      GOTO 199
  136 WRITE(IOUERR,436)
  436 FORMAT(
     *'+',12X,' Error occurred in reading the cold CO parameters.')
      GOTO 198
  137 WRITE(IOUERR,437)
  437 FORMAT(
     *'+',12X,' Error occurred in reading the',
     *        ' oxygenated fuels descriptive record.')
      GOTO 198
  138 WRITE(IOUERR,438)
  438 FORMAT(
     *'+',12X,' Error occurred in reading diesel sales fractions.')
      GOTO 198
  139 WRITE(IOUERR,439) IVALUE,IOOFBD,IYEND
  439 FORMAT('+',12X,I4,5A4,'ETP phase-in year (1995 to ',I4,')')
      GOTO 198
  140 WRITE(IOUERR,440) RVALUE,IOOFBD
  440 FORMAT('+',12X,F4.0,5A4,'ETP reduction rate (0 to 100%).')
      GOTO 198
  141 WRITE(IOUERR,441)
  441 FORMAT(
     *'+',12X,' The ETP rates for full and',
     *        ' partial diurnals must be equal.')
      GOTO 198
  142 WRITE(IOUERR,442)
  442 FORMAT(
     *'+',12X,' The ETP running loss rates for purge and',
     *        ' pressure failures must be equal.')
      GOTO 198
  143 WRITE(IOUERR,443)
  443 FORMAT(
     *'+',12X,' The ETP resting loss rates must all be equal.')
      GOTO 198
  144 WRITE(IOUERR,444)
  444 FORMAT(
     *'+',12X,' Error occurred in reading ETP parameters.')
      GOTO 198
  145 WRITE(IOUERR,445) IVALUE,IOOFBD
  445 FORMAT('+',12X,I4,5A4,'evaluation month 1 (Jan) or 7 (July)')
      GOTO 198
  146 WRITE(IOUERR,446)
  446 FORMAT(
     *'+',12X,' Diurnal emissions are zero when Tmax-Tmin<1. The'/
     *' ',12X,' minimum and maximum temperatures of the day must'/
     *' ',12X,' be used to get the proper diurnal output.')
      GOTO 199
  147 WRITE(IOUERR,447)
  447 FORMAT(
     *'+',12X,'User supplied market share of oxygenated fuels has ',
     *        'to be 100% for'/' ',12X,'a reformulated gas area.')
      GOTO 198
  148 WRITE(IOUERR,448)
  448 FORMAT(
     *'+',12X,'User supplied market share weighted average of ',
     *        'oxygen content has'/' ',12X,'to be at least 2.1% ',
     *        'for a reformulated gas area.')
      GOTO 198
  149 WRITE(IOUERR,449) IOOFBD
  449 FORMAT(
     *'+',12X,'Value ',5A4,'ASTM Volatility Class (A - E).')
      GOTO 198
  150 WRITE(IOUERR,450) RVALUE,IOOFBD
  450 FORMAT('+',12X,F5.0,5A4,'Pass/Purge/Pressure rate (0 to 100%).')
      GOTO 198
  151 WRITE(IOUERR,451) RVALUE
  451 FORMAT('+',12X,'The Pass/Purge/Pressure rate sum (',F4.0,') ',
     *               'does not equal 100%')
      GOTO 198
  152 WRITE(IOUERR,452)
  452 FORMAT('+',14X,
     *    'The I/M descriptive record field elements have changed',/,
     *' ','from MOBILE4.1 to MOBILE5. The three year values, 1) I/M ',/,
     *' ','start year, 2)Purge start year, and 3)Pressure start year',/,
     *' ','have been replaced. Please consult the MOBILE5 Users Guide.')
      GOTO 199
C
  153 WRITE(IOUERR,453)
  453 FORMAT(
     *'+',14X,'Refueling emissions in grams-per-gallon are only',/
     *' ',14X,'available using the 120 column descriptive output',/
     *' ',14X,'option (OUTFMT =3 or 5).  See MOBILE5 Users Guide',/
     *' ',14X,'chapters 2.1.15, 2.1.19 and 2.1.20 for more information.'
     *)
C
      GOTO 199
C
  154 WRITE(IOUERR,454)
  454 FORMAT(
     *'+',14X,'Refueling emissions for LDGV and LDGT after 1998',/
     *' ',14X,'model year have been reduced as a result of the',/
     *' ',14X,'Onboard Refueling Vapor Recovery Regulations (1994).')
C
      GOTO 199
C
  155 WRITE(IOUERR,455)
  455 FORMAT(
     *'+',14X,'User input of Onboard Refueling Vapor Recovery',/
     *' ',14X,'parameters has been ignored.')
C
      GOTO 199
C
  156 WRITE(IOUERR,456)
  456 FORMAT(
     *'+',14X,'Onboard Refueling Vapor Recovery Regulation (1994)',/
     *' ',14X,'effects have been disabled by the user.')
C
      GOTO 199
C
C  Aug-18-1994 @ CSC-bsg Request 2-462 Added QUITER messages 157-160, for 
C  non-FTP related error messages
C  Sep-14-1994 @ CSC-bsg Request 2-467 New Non-FTP Effects message and the
C  elimination of error message 157
C
  157 WRITE(IOUERR,457) RVALUE, IOOFBD
  457 FORMAT('+',12X,F7.2,5A4,' the effectiveness of test-and-repair',
     *      /' ',12X,'I/M programs (0.00 to 1.00)')
      GOTO 198
  158 WRITE(IOUERR,458) IVALUE
  458 FORMAT(
     *'+',14X,'Exhaust emission values include additional effects of ',
     *        'in-use'/
     *' ',14X,'driving behavior.  Model years beginning with ',I4,
     *        ' will '/
     *' ',14X,'have reduced effects due to a new certification ',
     *        'procedure.')
C
      GOTO 199
C
  159 WRITE(IOUERR,459) IVALUE,IOOFBD
  459 FORMAT('+',12X,I4,5A4,'First model year which will '/
     *' ',12X,'not have non-FTP effects.')
C
      GOTO 198
C
  160 WRITE(IOUERR,460)
  460 FORMAT(
     *'+',12X,'Error occurred in reading in-use driving behavior '/
     *' ',12X,'record.')
      GOTO 198
C
C
  161 WRITE(IOUERR,461)
  461 FORMAT(
     *'+',12X,'Error occurred in reading PPEFF.D file.')
      GOTO 198
C
C   Message 162 if one Pressure or Purge program is 3 and the other is not
C 
  162 WRITE(IOUERR,462)
  462 FORMAT(
     *'+',12X,'Pressure or Purge program have incompatable test types.',
     *   /, 13X, 'One is Hybrid and the other is not')
      GOTO 198
C      
C  Message 163 if error occurred in reading TTC.D file 
C 
  163 WRITE(IOUERR,463)
  463 FORMAT(
     *'+',12X,'ERROR occurred in reading TTC.D file.')
      GOTO 198
C      
C  New message 164 if error occurred in reading "I/M Control Flag Record"
C 
  164 WRITE(IOUERR,464)
  464 FORMAT(
     *'+',12X,'ERROR occurred in I/M Control Flag Record')
      GOTO 198
C 
C  New message 165 if value out of bounds in the TTC Record"  
C
  165 WRITE(IOUERR,465) RVALUE
  465 FORMAT('+',14X,G7.2,' is out of bounds for the Effectiveness ',
     *      'in the TTC record.',/,' ',14X,'Reset to 100%') 
      GOTO 199
C      
C  New message 166 if I/M + TTC credit exceed 100%
C 
  166 WRITE(IOUERR,466)
  466 FORMAT('+',14X,'I/M and TTC credit exceed 100%.',
     *      ' Credit will be reset to 100%')
      GOTO 199
C      
C  New message 167: ATPUER is out of bounds. 
C 
  167 WRITE(IOUERR,467) RVALUE, IOOFBD 
  467 FORMAT(
     *'+',12X,F7.2, 5A4,' the effectiveness of test-and-repair',/
     *' ',12X,'ATP programs (0.00 to 1.00)')
      GOTO 198
C      
C  New message 168: PPGUER is out of bounds. 
C 
  168 WRITE(IOUERR,468) RVALUE, IOOFBD
  468 FORMAT(
     *'+',12X,F7.2, 5A4,' the effectiveness of test-and-repair',/
     *' ',12x,'purge/pressure programs (0.00 to 1.00)')
      GOTO 198
C      
C  New message 169: DGPIEY is out of bounds. 
C 
  169 WRITE(IOUERR,469) IVALUE, IOOFBD
  469 FORMAT(
     *'+',12X,I4, 5A4,' the detergent gas phase-in end year ',/
     *' ',12x,'(this should be later than 1995.)')
      GOTO 198
C
  170 WRITE(IOUERR,470)
  470 FORMAT(
     *'+',14X,'Exhaust emissions for gasoline fueled vehicles',/
     *' ',14X,'beginning in 1995 have been reduced as a result of',/
     *' ',14X,'Gasoline Detergent Additive Regulations (1994).')
C
      GOTO 199
C
C  Output warning here if user stumbles upon/uses a "hidden" feature
C
  195 WRITE(IOUERR,495)
  495 FORMAT(
     *' ',14X,'There is no user application of this feature'/
     *' ',14X,'which would be acceptable for submission to EPA.')
      GOTO 199
C
C  Unknown error message numbers land here
C
  197 WRITE(IOUERR,497) MESSAG
  497 FORMAT('+',12X,' Message code ',I4,' is unknown.')
C  Messages 60 and 61 produce excessive output if not stopped here.  The other
C  cases stopping execution here are considered "fatal" errors:  The user must
C  fix the error-generating input and rerun.
C
  198 IF(MCODE.EQ.53.OR.MCODE.EQ.60.OR.MCODE.EQ.61.OR.MCODE.EQ.82.OR.
     *   MCODE.EQ.85.OR.MCODE.EQ.89.OR.MCODE.EQ.97.OR.MCODE.EQ.102.OR.
     *   MCODE.EQ.107.OR.MCODE.EQ.113.OR.MCODE.EQ.118.OR.
     *   MCODE.EQ.161.OR.MCODE.EQ.162.OR.MCODE.EQ.163.OR.MCODE.EQ.164)
     *   STOP
      INERR=INERR+1
C
  199 RETURN
      END
