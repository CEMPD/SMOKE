      SUBROUTINE GETERQ(INPREC,INERR,RC) 
C 
C  GETERQ reads in, performs range checks on, converts to gpm if 
C  IV = HDV and stores replacement emission factor rates & their 
C  associated parameters. Emission line parameters for light duty 
C  gas vehicles IV=1-5 with pre-1981 model years; all model year 
C  light duty diesel vehicles, IV=14(light duty diesel vehicles), 
C  15(light duty diesel trucks1&2), and 28(light duty diesel 
C  trucks3&4); all model year motorcycles, IV=24; and all model 
C  year heavy duty vehicles, IV=6-13(heavy duty gas), 16-23(heavy 
C  duty diesel), and 25-27(diesel and gas buses). (The subroutine 
C  EF_CASE.FOR is the vehicle type/ model year/ technology/ ... routing 
C  method. It directs calls to a set of subroutines which have access 
C  to the appropriate emission line parameters which are used to calculate 
C  vehicle "basic" exhaust emissions.) 
C 
C  Called by DSPTCH. 
C 
C  Calls CHKLAB, HDCGPM, QUITER, STOERQ, UCCOMP and YRTEST. 
C 
C  Changes: (Last change first) 
C 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C  23 Aug 00 @EPA-bag Bug62  conditional for light duty gas vehicles with 
C            model years greater than 1993 changed to include start emission 
C            parameters 
C   2 Aug 00 @EPA-ddj Bug62  Changed parameter sent to QUITER for 
C            1981-1993 LDG entered with label BASIC EFS. 
C   3 Jul 00 @EPA-ddj Bug62  Added warning for 1981-1993 LDG entered 
C            with label BASIC EFS. 
C  31 May 00 @EPA-elg   Changed fixed READ format to Free READ format. 
C  11 Apr 00 AIR Task 01: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C            so that Lahey Standard for keyboard/monitor I/O is followed. 
C  20 Jan 00 @DynTel-ddj 1-010 Revised comments regarding values of 
C            5 and 6 for NEWFLG. 
C  28 Jul 99 @DynTel-MLA 2-000 Removed the USE LABELS statement. The 
C            module subroutine CHKLAB was removed from that module 
C            and converted to an ordinary external procedure. 
C   1 Jul 99 @DynTel-HXQ 2-705 Add HDCGPM_FLAG. 
C  15 Jun 99 @DynTel-HXQ  Changing the format of user input. 
C  11 May 99 @DynTel-MLA 2-674  Modified calls to CHKLAB because the 
C            calling sequence for that subroutine was changed. 
C   2 Mar 99 @DynTel-MLA 2-659  Cleaned up comments. Changed IOUREP to 
C            IOUERR in error messages. 
C  08 Oct 98 @DynTel-ZK 2-000 Explicit Typing 
C  31 Jul 98  @DynTel-BAG 2-000  Comment statements 
C 
C  8 Apr 98  @DynTel-RJD 2-654  Added an extra field to the BER parameter 
C             record to determine if the record pertains to starting or 
C             running emissions. 
C  Oct-30-97 @DynTel-gjr 2-645 Eliminating option for interactive input. 
C 
C  11 Jun 97 DynTel@MLA 2-622  Altered code to read Mobile6 transition 
C            input formats. 2 digits are allowed for vehicle types. 
C            Added IVTYPE.I for Mobile6. 
C 
C  Input on call: 
C 
C    Constants:  
C    MAXVEH from include file IVTYPE.I 
C 
C    argument list: INERR 
C    common blocks: 
C    /BASEQ4/ MAXCT 
C    /FLAGS2/ NEWFLG 
C    /IOUCOM/ IOUGEN,IOUOUT 
C    /YEARS4/ IYSTART,IYEND 
C    /VVNAMS/ VVHEAVY, VVLDGAS  (IVTYPE.I) 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C    /BASEQ3/ KINKS 
C    /BASEQ4/ NEWPAR,BERNEW,NEWFIT,NEWCT,NEWOML 
C    /FLAGS2/ NEWFLG 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  A50NEW   R    above 50K miles deterioration rate, nonzero for 1981+ LDGV 
C                and 1991+ LDGT, HC/CO cases only (local variable: not used in 
C                HDCGPM) 
C  KINK50   I    kinked at 50K miles ef curve flag: 1 = no, 2 = yes 
C  NEWPOL   I    new emission rate's pollutant code 
C  NEWREG   I    new emission rate's region code 
C  SECTION  I    Dummy argument for CHKLAB. 
C  STDLAB   C    Dummy argument for CHKLAB. 
C 
C  Notes: 
C 
C  Kinked EFs for LDGTs were added for MOBILE5. 
C 
      IMPLICIT NONE  
      INCLUDE 'BASEQ3.I' 
      INCLUDE 'BASEQ4.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'YEARS4.I' 
C 
      INTEGER, EXTERNAL :: JUNIT   
C 
      CHARACTER   INPREC*(*) 
      INTEGER     INERR 
      INTEGER     RC 
      INTEGER     BER 
      LOGICAL     EXTDATA 
      CHARACTER   FILENAME*80 
      CHARACTER   LABEL*19 
      INTEGER     LABTYPE 
      CHARACTER   LASTREC*150 
      CHARACTER   M6REC*150 
      INTEGER     SECTION 
      INTEGER     SOURCE 
      CHARACTER   STDLAB*40 
      INTEGER IFAC, KINK50, NEWREG, NEWPOL 
      REAL A50NEW  
      
ccs   Add declarations of M6DIR and M6FILE
      CHARACTER(LEN=200) M6DIR
      CHARACTER(LEN=280) M6FILE
C 
C     Note: The rules for the federal evaporative test procedure have 
C     been finalized and the parameters of the final rulemakeing have  
C     been incorporated into MOBILE6.  MOBILE6 does not allow users to 
C     enter any parameters for the federal evap test procedure. This 
C     means the values of 3 and 4 for NEWFLG that were valid in MOBILE5 
C     are no longer acceptable in MOBILE6. The values of 1 (default BER) 
C     and 2 (user-supplied BER) are still acceptable. 
C 
C     Note that the user no longer has the option to disable the Clean 
C     Air Act Amendment (CAAA).  Therefore, values 5 (disable 1990 CAAA 
C     requirements with default BER), and 6 (disable 1990 CAAA 
C     requirements with user-supplied BER) are no longer applicable. 
C 
C     THE FEDERAL EVAPORATIVE TEST PROCEDURE IS NOT FULLY IMPLEMENTED 
C     YET IN MOBILE6. THE CURRENT VERSION OF MOBILE6 STILL ALLOWS THE 
C     ENTRY OF EVAP TEST PARAMTERS. 
C 
      RC = 0 
      IF (NEWFLG==1) THEN 
        NEWFLG = 2 
      ELSEIF (NEWFLG==3) THEN 
        NEWFLG = 4 
      ELSEIF (NEWFLG==5) THEN 
        NEWFLG = 6 
      ELSE 
        CALL QUITER(5.0, NEWFLG, 1, INERR) 
        RC = 1 
        GOTO 99 
      ENDIF 
C 
C     Find out whether we are reading from an external data file. 
C                  
      IF (INPREC(20:20)=='@') THEN 
        EXTDATA = .FALSE. 
        SOURCE = IOUGEN 
        M6REC = INPREC 
      ELSE 
C 
C       Read the data from an external data file. 
C 
        EXTDATA = .TRUE. 
        SOURCE = JUNIT() 
C 
C       Find the token number of the record label that must appear 
C       on a record somewhere in the header of the external data file. 
C 
        CALL CHKLAB('BASIC EFS',STDLAB,BER,SECTION) 
C       The next two assignment statements are present solely to suppress 
C       compiler generated warnings. They have no functional importance. 
        LABEL = STDLAB(1:1) 
        LABTYPE = SECTION 
C 
C       Open the external data file. 
C 

ccs        CALL UCCOMP(INPREC(21:), FILENAME) 
        CALL LJCOMP(INPREC(21:), FILENAME)

ccs     All M6 input files live in SMK_M6PATH
        CALL GETENV( "SMK_M6PATH", M6DIR )
        M6FILE = M6DIR( 1:LEN_TRIM( M6DIR ) ) // '/' // FILENAME

        OPEN(SOURCE,FILE=M6FILE,STATUS='OLD',ACTION='READ', 
     *       IOSTAT=RC,ERR=80) 
        WRITE(IOUREP,120) FILENAME 
  120   FORMAT(/'* Reading non-default BASIC EMISSION RATES from ', 
     &          'the following external',/'* data file: ',A) 
C 
C       Read the required header line from the file. 
C 
        CALL NXTREC(SOURCE,IOUREP,M6REC,INERR,RC) 
        IF (RC.NE.0) GOTO 85 
C 
C       Check for a valid labeled input record. 
C 
        CALL UCCOMP(M6REC(1:19),LABEL) 
        CALL CHKLAB(LABEL,STDLAB,LABTYPE,SECTION) 
        IF (LABTYPE.NE.BER) GOTO 90 
C 
      ENDIF 
C 
C     Read the data field from the input record. 
C 
      READ(M6REC,130,ERR=92,END=92,IOSTAT=RC) NEWCT 
  130 FORMAT(T22,BN,I3) 
      IF(NEWCT<1.OR.NEWCT>MAXCT) CALL QUITER(0.,NEWCT,23,INERR) 
C 
C     Find the first data line. 
C 
      CALL NXTREC(SOURCE,IOUREP,M6REC,INERR,RC) 
      IF (RC.NE.0) GOTO 98 
      BACKSPACE(SOURCE) 
C 
C  Read in and process one record at a time. 
C 
      KINKS=1 
      DO 50 IFAC=1,NEWCT 
      READ(SOURCE,*,ERR=98,END=98) NEWOML(IFAC),NEWREG,NEWVEH,NEWPOL 
     *                        ,NEWMYF,NEWMYL,ZMLNEW,DRNEW,A50NEW 
C 
      IF(NEWOML(IFAC).NE.'R'.AND.NEWOML(IFAC).NE.'r'.AND. 
     *     NEWOML(IFAC).NE.'S'.AND.NEWOML(IFAC).NE.'s') 
     *     CALL QUITER(0.,0,71,INERR)   
      IF(NEWREG<1.OR.NEWREG>2) CALL QUITER(0.,NEWREG,62,INERR) 
      IF(NEWVEH<1.OR.NEWVEH>MAXVEH) 
     &     CALL QUITER(0.,NEWVEH,25,INERR) 
      IF(NEWPOL<1.OR.NEWPOL>3) CALL QUITER(0.,NEWPOL,26,INERR) 
      CALL YRTEST(NEWMYF,27,IYSTART,IYEND,INERR) 
      CALL YRTEST(NEWMYL,27,IYSTART,IYEND,INERR) 
      IF(NEWMYL<NEWMYF) CALL QUITER(0.,0,29,INERR) 
      KINK50=1 
      IF( (VVLDGAS(NEWVEH)==1) .AND. (NEWMYL>=1981.AND.NEWMYL<=1993) ) 
     &      CALL QUITER(0.0,IFAC,154,INERR) 
      IF(ZMLNEW<=0.0) CALL QUITER(ZMLNEW,0,30,INERR) 
      IF(DRNEW<0.0) CALL QUITER(DRNEW,0,31,INERR) 
      IF(KINK50==2.AND.A50NEW<0.0) CALL QUITER(A50NEW,0,31,INERR) 
C 
      IF (VVLDGAS(NEWVEH)==1.AND.NEWMYL>=1994) THEN 
        CALL QUITER(0.0,IFAC,74,INERR) 
      ENDIF 
C 
C  If the emission rate is for a heavy duty vehicle class, convert the rate's 
C  units from grams per brakehorsepower to grams per mile. 
C 
      IF(VVHEAVY(NEWVEH)==1.AND.HDCGPM_FLAG==1) CALL HDCGPM 
C 
C  Save screened data in temporary storage areas. 
C 
      NEWPAR(1,IFAC)=NEWREG 
      NEWPAR(2,IFAC)=NEWVEH 
      NEWPAR(3,IFAC)=NEWPOL 
      NEWPAR(4,IFAC)=NEWMYF 
      NEWPAR(5,IFAC)=NEWMYL 
      BERNEW(1,IFAC)=ZMLNEW 
      BERNEW(2,IFAC)=DRNEW 
      IF(KINK50==1) BERNEW(3,IFAC)=-1.0 
      IF(KINK50==2) BERNEW(3,IFAC)=A50NEW 
   50 CONTINUE 
C 
      IF (EXTDATA) CLOSE (SOURCE) 
C 
C  Store new ef equations for later access by BEF. 
C 
      CALL STOERQ(INERR) 
C 
   60 GOTO 99 
C 
C     Error opening the external data file. 
C 
   80 WRITE(IOUERR,180) FILENAME 
  180 FORMAT(/'*** ERROR: unable to open the following external ', 
     &        'data file:',/'*** ',A) 
      INERR = INERR + 1 
      GOTO 99 
C 
C     Error reading the header of the external data file. 
C 
   85 SELECT CASE (RC) 
        CASE (:-1) 
          WRITE(IOUERR,185) FILENAME 
  185     FORMAT(/'*** ERROR: no header record was found in ', 
     &            'the following external data file:',/'***',A) 
          RC = 2 
          INERR = INERR + 1 
        CASE (1:) 
          WRITE(IOUERR,186) FILENAME 
  186     FORMAT(/'*** ERROR reading the header of the ', 
     &            'following external data file:',/'*** ',A) 
        CASE DEFAULT 
      END SELECT 
      GOTO 99 
C 
C     Invalid label on the header record. 
C 
   90 WRITE(IOUERR,190) M6REC 
  190 FORMAT(/'*** ERROR: the following record is not allowed in ', 
     &        'the header of the external data file:', 
     &      /'*** ',A) 
      RC = 3 
      INERR = INERR + 1 
      GOTO 99 
C 
   92 WRITE(IOUERR,192) M6REC 
  192 FORMAT(/'*** ERROR reading the data field of the following ', 
     &        'labeled input record:', 
     &       /'*** ',A) 
      INERR = INERR + 1 
      GOTO 99 
C 
   98 CALL QUITER(0.,0,128,INERR) 
      IF (RC<0) THEN 
        WRITE(IOUERR,198) 
  198   FORMAT('*** (Unexpected end of file.)') 
      ELSE 
        BACKSPACE(SOURCE) 
        READ(SOURCE,199,ERR=99,END=99) LASTREC 
  199   FORMAT(A) 
        WRITE(IOUERR,200) LASTREC 
  200   FORMAT('*** The following record caused the error:', 
     &        /'*** ',A) 
      ENDIF 
C 
   99 RETURN 
      END 
