      SUBROUTINE SMKDRIVER( INPFILE, LDEV, EFLAG )
C 
C  DRIVER calls MOBILE as a subroutine. 
C 
C  Calls MOBILE. 
C 
C  Changes: (Last change first) 
C 
C  12 Jan 02 AIR Task Bug383: Moved creation of M6ERROR.TXT so that 
C            improper program aborts do not leave file uneditable. 
C  16 Jan 02 AIR Task 11: Added INERR to PROCHDR parameter list to 
C            improve error handling. Improved TOTERR error handling. 
C            Fixed unassigned output on unit 2. 
C  10 Jan 01 @EPA-bag bug371 change error message filename 
C  13 Dec 01 AIR Task 29: Added INQUIRE to error message logic. 
C  10 Apr 11 AIR Task 11: Redirected error messages to IOUERR file. 
C  29 Sep 00 AIR Task 07: Corrected syntax error in 150 FORMAT. 
C  02 Oct 00 @EPA-djb Bug230, Removed unused I/O devices. 
C   08 aUG 00 @EPA-bag,  Placed a runtime diagnostic using DATE_AND_TIME 
C             intrinsic subroutine 
C   15 Jun 00 @EPA-bag,  Bug120 Removed unused IOUSSH and SSOPEN variables 
C  11 Apr 00 AIR Task 01: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C            so that Lahey Standard for keyboard/monitor I/O is followed. 
C            Also changed IOUGEN,IOUREP, and IOUERR. Disabled output of runtime. 
C  10 Mar 00 @DynTel-MLA 1-013  database code 
C   2 Aug 99 @DynTel-MLA 2-684  New version of DRIVER for Mobile6. 
C 
C  Local array subscripts: 
C 
C  Local variable / array dictionary: 
C 
C   Name      Type              Description 
C  ------     ----  ------------------------------------------------ 
C  BACK         L   Tells intrinsic IDX to search backw. from end of string 
C  BATMODE      L   TRUE if we are processing a batch file 
C  BATREC       C   An input (data) record from the batch file 
C  BUFFER1      C   Used to compose message strings 
C  BUFFER2      C   Used to compose message strings 
C  CMDOUTPUT    C   Report file name spec'd interactively or in batch file 
C  CMDACTION    C   Open action specified interactively or in batch file 
C  CMIN         C   Character representation of NMIN 
C  CSEC         C   Character representation of NSEC 
C  ELAPSED      I   Processing time elapsed 
C  ENDTIME      I   Time at end of processing 
C  FEXIST       L   TRUE if specified file exists 
C  FOPEN        L   TRUE if specified file is open 
C  IDX          I   Pointer to a position in a character string 
C  INERR        I   Number of errors reported by MOBILE 
C  INPFILE      C   Input file name 
C  INPLINE      C   One (data) line from the input file 
C  LABEL        C   Label field from an input file record 
C  NHOUR        I   Hour value returned by TIMER 
C  NMIN         I   Minute value returned by TIMER 
C  NSEC         I   Seconds value returned by TIMER 
C  NUMDID       I   Number of batch files processed 
C  PASS         I   1 on first pass through file prompt loop, 2 otherwise 
C  PTR1         I   Pointer to a pos. in a character string (for NXTTOK) 
C  PTR2         I   Pointer to a pos. in a character string (for NXTTOK) 
C  RC           I   Return code 
C  ROOT         C   Input file name (with extension removed) 
C  RUNMODE      L   TRUE if input is regular input file (not a batch file) 
C  STARTTIME    I   Time at start of processing 
C 
C  Notes: 
C 
C  29 Jun 98 @Dyntel-BAG 2-661 Output filenames have the extensions 
C            Descriptive output file     : '.txt' 
C            Spreadsheet tab delimited   : '.prn' 
C            Spreadsheet comma delimited : '.csv' 
C            (Spreadsheet sheet output files are opened in OUTPUT.) 
C 
      USE DATABASE 
C 
      IMPLICIT NONE 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'VDATA.I' 
C     
      INTEGER, EXTERNAL :: JUNIT 
C 
      EXTERNAL BD01,     BD01P5,   BD02P5,   BD03,     BD04,      
     &         BD04P5,   BD05,     BD06,     BD07,     BD08P5,    
     &         BD10,     BD11,     BD12,     BD15,     BD16,  
     &         BD17,     BD19,     BD20,     BD23,     BD24,    
     &         BD25,     BD27,     BD28,     BD29,     BD30,      
     &         BD31,     BD33,     BD34,     BD35,     BD36,  
     &         BD37,     BD38,     BD39,     BD40,     BD41,      
     &         BD42,     BD43,     BD44,     BD45,     BD46,      
     &         BD47,     BD48,     BDEV01,   BDSPEED1, BDSPEED2,  
     &         BDUNIT,   BD_AC,    BD_IM1,   BD_IM2,   BD_NGV,    
     &         BD_SFTP,  BD_TIER2, BD_T2PM   
C 
      CHARACTER(LEN=80), INTENT (IN) :: INPFILE    ! Mobile6 input file name
      INTEGER          , INTENT (IN) :: LDEV       ! SMOKE log file unit number
      LOGICAL          , INTENT(OUT) :: EFLAG      ! true: an error has occured
C 
      LOGICAL,PARAMETER         ::  BACK=.TRUE. 
      LOGICAL                   ::  BATMODE 
      CHARACTER*150             ::  BATREC 
      CHARACTER*150             ::  BUFFER1 
      CHARACTER*150             ::  BUFFER2 
      CHARACTER*80              ::  CMDOUTPUT 
      CHARACTER*10              ::  CMDACTION 
      LOGICAL                   ::  FEXIST 
      LOGICAL                   ::  FOPEN 
      INTEGER                   ::  IDX 
      INTEGER                   ::  INERR 
ccs      CHARACTER*80              ::  INPFILE 
      CHARACTER*170             ::  INPLINE 
      CHARACTER*20              ::  LABEL 
      CHARACTER*10              ::  CSTARTTIME 
      CHARACTER*10              ::  CENDTIME 
      CHARACTER*8               ::  CDATE 
      CHARACTER*5               ::  CZONE 
      INTEGER, DIMENSION(8)     ::  ISTARTTIME 
      INTEGER, DIMENSION(8)     ::  IENDTIME 
      INTEGER, DIMENSION(3)     ::  ITOTAL_TIME 
      INTEGER                   ::  NUMDID 
      INTEGER                   ::  PASS 
      INTEGER                   ::  PTR1 
      INTEGER                   ::  PTR2 
      INTEGER                   ::  RC 
      INTEGER                   ::  TOTERR 
      CHARACTER*80              ::  ROOT 
      LOGICAL                   ::  RUNMODE 

ccs   Add declarations of M6DIR and M6FILE
      CHARACTER(LEN=50) M6DIR
      CHARACTER(LEN=80) M6FILE
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
C On the MAC, set the floating point control reg to control  
C rounding precision used for floating point calculations. 
C Uncomment the following lines for MacFortran II version only.  
C MACC Begin system dependent code for Macintosh 
C MAC      IMASK = getfpcontrol() 
C MAC      IMASK = (IMASK.AND.RPCLEAR) 
C MAC      IMASK = (IMASK.OR.RPSINGLE) 
C MAC      CALL setfpcontrol(IMASK) 
C MACC End of system dependent code for Macintosh 
C Suggested change for different computers 
C on the MAC IOUGEN=9 and on the IBM IOUGEN=1 
C on the MAC IOUIN=9 and on the IBM IOUIN=5 
C on the MAC IOUOUT=9 and on the IBM IOUOUT=6 
C on the MAC IOUERR=6 and on the IBM IOUERR=2 
C 
C  IOUGEN  : Device used to read the command file. 
C  IOUIMD  : Device used to read the file containing I/M credits. 
C  IOUIN   : Device used to read interactive input from the screen. 
C  IOUALT  : Device used to read data from external files. 
C  IOUBAT  : Device used to read commands from a batch file. 
C  IOUXIM  : Device used to read the external file containing I/M parameters. 
C  IOUREP  : Device used to write the descriptive output. 
C  IOUERR  : Device used to write error messages. 
C  IOUOUT  : Device used to write to the screen. 
C  IOUNULL : Device used to dump unneeded output. 
C  DBUNIT(7) : Devices used to write database output. 
C 
C See IOUCOM and Block Data 16 for initializations. 
C 
C     Input Units: 
      IOUGEN = 0       ! opened in SMKDRIVER.F 
      IOUXIM = 0       ! opened in GETXIM.F
      IOUIMD = 0       ! opened in GETM5IM.F and RDT12C.F
ccs      IOUIN  =  5 
      IOUALT = 0       ! opened in DBOPTS.F, GET94E.F, GET94P.F, GETATHAPS.F, 
                       !   GETCONV.F, GETDSACT.F, GETFVMT.F, GETHSACT.F, 
                       !   GETMAR.F, GETMPG.F, GETNGVEF.F, GETNGVFR.F, GETNH3BER.F, 
                       !   GETNH3SDR.F, GETPMDDR1.F, GETPMDDR2.F, GETPMDZML.F, 
                       !   GETPMGDR1.F, GETPMGDR2.F, GETPMGZML.F, and GETREG.F  
      IOUBAT = 0        ! opened in SMKDRIVER.F  
C 
C     NOTE: Unit numbers 51-57 are reserved for the database output 
C           tables and are initialized in MODULE DATABASE. 
C     Output units: 
      IOUREP = 0       ! opened in PROCHDR.F
      IOUERR = 0       ! opened in SMKDRIVER.F
      IOUOUT = LDEV    ! SMOKE log file
      IOUNULL = -1 
 
      IOUPM  = 0        ! opened in PROCHDR.F 
      IOUAT  = 0        ! opened in PROCHDR.F 
      IOUSS  = 0        ! opened in PROCHDR.F 
C 
C     Initialize total errors 
C 
      TOTERR = 0 
C 
C     Initialize the name of the program. 
C 
      INCLUDE 'VNAME.I' 
C 
      WRITE(IOUOUT,200) VERSION 
  200 FORMAT (' ',A80) 
C 
C     Get the name of the Mobile6 input file. It defaults to MOBILE6.IN 
C     if the user enters a blank line in response to the first prompt. 
C 
ccs      WRITE(IOUOUT,100) 
ccs  100 FORMAT(' Enter the name of the Mobile6 input file:'/) 
C 
C     Make sure the file exists. (PASS counts the number of passes 
C     through the next DO loop.) 
C 
ccs      PASS = 1 
ccs      DO  ! Forever 
C 
ccs         READ(IOUIN,105) INPLINE 
ccs  105    FORMAT(A) 
C 
C        If the blank line is entered on the first pass through the loop, 
C        use the default input file name of MOBILE6.IN and the default 
C        descriptive output file name of MOBILE6.TXT (generated later, 
C        unless another output file name is specified in the header of 
C        the input file). If a blank line is entered after the first pass 
C        through the loop, stop the program. 
C 
ccs         IF (INPLINE == ' ' .AND. PASS == 1) THEN 
ccs            INPLINE = 'MOBILE6.IN' 
ccs         ENDIF 
C 
ccs         PASS = 2 
C 
ccs         IF (INPLINE == ' ') THEN 
ccs            GOTO 99 
ccs         ELSE 
C 
C           Look for an input file name (and a possible ouptut file name 
C           and open action) on the line entered by the user in response 
C           to the prompt. 
C 
ccs            PTR2 = 0 
ccs            CALL NXTTOK(INPLINE,INPFILE,PTR1,PTR2) 
ccs            CALL NXTTOK(INPLINE,CMDOUTPUT,PTR1,PTR2) 
ccs            CALL NXTTOK(INPLINE,CMDACTION,PTR1,PTR2) 
C 
ccs         ENDIF 
C 
C        If the file name provided by the user contains a period, use the 
C        name exactly as given. If the name provided does not contain a 
C        period, look FIRST for a file with the given name and an extension 
C        of '.IN' and THEN, if that file is not found, look for a file with 
C        the given name and no added extension. 
C 
ccs         CALL UCCOMP(INPFILE,ROOT) 
ccs         CALL LJCOMP(INPFILE,ROOT)

ccs         IDX = INDEX(ROOT,'.',BACK) 
ccs         IF (IDX == 0) THEN 
ccs            INPFILE = ROOT(1:LEN_TRIM(ROOT)) // '.IN' 
ccs         ELSE 
ccs            INPFILE = ROOT 
ccs         ENDIF 
ccs         INQUIRE(FILE=INPFILE,EXIST=FEXIST) 
ccs         IF (.NOT. FEXIST .AND. IDX == 0) THEN 
ccs            INPFILE = ROOT 
ccs            INQUIRE(FILE=INPFILE,EXIST=FEXIST) 
ccs         ENDIF 
ccs         IF (FEXIST)  EXIT 
C 
ccs         WRITE(IOUOUT,110) INPFILE(1:LEN_TRIM(INPFILE)) 
ccs  110    FORMAT(' The file ',A,' does not exist.', 
ccs     &         /' Enter another file name. (Enter a blank line', 
ccs     &          ' to stop the program.)'/) 
C 
ccs      END DO 
C 
C     Open the IOUERR file 
C 
ccs   All M6 input files live in SMK_M6PATH 
      CALL GETENV( "SMK_M6PATH", M6DIR ) 

ccs   Make sure SMK_M6PATH is set
      IF( M6DIR == ' ' ) THEN
          WRITE( *,* ) 'ERROR: Environment variable SMK_M6PATH ' //
     &                 'is not set.'
          EFLAG = .TRUE.
          RETURN
      END IF

      M6FILE = TRIM( M6DIR ) // '/' // 'M6ERROR.TXT' 
       
      IOUERR = JUNIT() 
      OPEN(IOUERR,FILE=M6FILE) 
C 
C     Open the input file for reading. 
C 
ccs      IOUGEN = JUNIT() 
ccs      OPEN(IOUGEN,FILE=INPFILE,STATUS='OLD',ACTION='READ',IOSTAT=RC) 
ccs      IF (RC /= 0) GOTO 80 
C 
C     Read the first line and determine whether it is a batch file. Skip 
C     blank lines and comments. (Whether the file is a batch file or a 
C     regular input file, it will be reopened later. Don't echo comments 
C      here -- i.e. pass IOUNULL to NXTREC -- or they will appear twice.) 
C 
ccs      CALL NXTREC(IOUGEN, IOUNULL, INPLINE, INERR, RC) 
ccs      IF (RC /= 0) GOTO 85 
ccs      CLOSE(IOUGEN) 
C 
C     Look for the batch/run mode header on the first data line of the file. 
C 
ccs      CALL UCCOMP(INPLINE(1:19),LABEL) 
ccs      BATMODE = (LABEL == 'MOBILE6 BATCH FILE') 
ccs      RUNMODE = (LABEL == 'MOBILE6 INPUT FILE') 
ccs      IF (.NOT. BATMODE .AND. .NOT. RUNMODE) GOTO 90 
C 
ccs      IF (BATMODE) THEN 
ccs         WRITE(IOUOUT,115) INPFILE(1:LEN_TRIM(INPFILE)) 
ccs      ELSE IF (RUNMODE) THEN 
ccs         WRITE(IOUOUT,116) INPFILE(1:LEN_TRIM(INPFILE)) 
ccs      ENDIF 
ccs  115 FORMAT(' Batch file name: ',A) 
ccs  116 FORMAT(' Input file name: ',A) 
C 
C     If batch mode, reopen the file as the batch input; skip line 1. 
C     (Comments beginning with '>' will be echoed to the terminal.) 
C 
ccs      DBBATCH = ' ' 
ccs      IF (BATMODE) THEN 
ccs         DBBATCH = INPFILE 
ccs         IOUBAT = JUNIT() 
ccs         OPEN(IOUBAT,FILE=DBBATCH,STATUS='OLD',ACTION='READ',IOSTAT=RC) 
ccs         IF (RC /= 0) GOTO 80 
ccs         CALL NXTREC(IOUBAT, IOUOUT, INPLINE, INERR, RC) 
ccs         IF (RC /= 0) GOTO 85 
ccs      ENDIF 
C 
C     Collect timing statistics (start timer). 
C 
      CALL DATE_AND_TIME(CDATE,CSTARTTIME,CZONE,ISTARTTIME) 
      WRITE(IOUOUT,120) CSTARTTIME(1:2),CSTARTTIME(3:4),CSTARTTIME(5:10) 
  120 FORMAT(/'    Processing start time is ',A2,':',A2,':',A6,'.'/) 
C 
C     Process each input file. 
C 
ccs      NUMDID = 0 
      DBNFILE = 0 
ccs      DBSTRNUM = 0 
ccs      DBLINENUM = 0 
      CALL INILAB() 
ccs      MAINLOOP: DO WHILE (BATMODE .OR. RUNMODE) 
C 
ccs         RUNMODE = .FALSE.  
C 
ccs         IF (BATMODE) THEN 
C 
C           Get the next non-blank, non-comment line from the batch file. 
C           (Comments beginning with '>' will be echoed to the terminal.) 
C 
ccs            CALL NXTREC(IOUBAT, IOUOUT, BATREC, INERR, RC) 
ccs            IF (RC < 0) EXIT 
ccs            IF (RC > 0) THEN 
ccs               INPFILE = DBBATCH   ! For error message code only. 
ccs               GOTO 85 
ccs            ENDIF 
C 
C           Look for the input file name, an output file name, and an open 
C           action for the output file. The last two items are optional. 
C 
ccs            PTR2 = 0 
ccs            CALL NXTTOK(BATREC,INPLINE,PTR1,PTR2) 

ccs            CALL UCCOMP(INPLINE,INPFILE) 
ccs            CALL LJCOMP(INPLINE,INPFILE)

ccs            CALL NXTTOK(BATREC,CMDOUTPUT,PTR1,PTR2) 
ccs            CALL NXTTOK(BATREC,CMDACTION,PTR1,PTR2) 
C 
ccs         END IF 
C 
C        Make sure the file exists; open it. 
C 
         DBNFILE = DBNFILE + 1 
ccs         IDX = INDEX(INPFILE,'.',BACK) 
ccs         IF (IDX == 0) THEN 
ccs            ROOT = INPFILE 
ccs            INPFILE = ROOT(1:LEN_TRIM(ROOT)) // '.IN' 
ccs         ENDIF 
ccs         INQUIRE(FILE=INPFILE,EXIST=FEXIST) 
ccs         IF (.NOT. FEXIST .AND. IDX == 0) THEN 
ccs            INPFILE = ROOT 
ccs            INQUIRE(FILE=INPFILE,EXIST=FEXIST) 
ccs         ENDIF 
ccs         IF (.NOT. FEXIST) THEN 
ccs            WRITE(BUFFER1,125) DBNFILE 
ccs  125       FORMAT(' *** ERROR: Unable to find the following', 
ccs     &             ' input file (file ',I4,'):') 
ccs            CALL LJCOMP(BUFFER1,BUFFER2) 
ccs            WRITE(IOUOUT,126) BUFFER2(1:LEN_TRIM(BUFFER2)), 
ccs     &                 INPFILE(1:LEN_TRIM(INPFILE)) 
ccs  126       FORMAT(' ',A,/' *** ',A/) 
ccs            CYCLE MAINLOOP 
ccs         END IF 
C 
C        In batch mode, display the name of each input file. 
C 
ccs         IF (BATMODE) WRITE(IOUOUT,128) INPFILE(1:LEN_TRIM(INPFILE)) 
ccs  128    FORMAT(' * Input file:  ',A)     
C 
         IOUGEN = JUNIT()
         OPEN(IOUGEN,FILE=INPFILE,STATUS='OLD',ACTION='READ',IOSTAT=RC) 
         IF (RC /= 0) GOTO 80 
         DBINPUT = INPFILE 
C 
C        Read the header of the input file. Set the output flags and 
C        open the output files, using default names if necessary. 
C 
         CALL PROCHDR(INPFILE,CMDOUTPUT,CMDACTION,INERR,RC) 
         IF (RC /= 0) THEN 
            WRITE(BUFFER1,130) RC,INPFILE 
  130       FORMAT(' *** ERROR ',I9,' processing the header of', 
     &             ' input file ',A80) 
            CALL LJCOMP(BUFFER1,BUFFER2) 
            WRITE(IOUOUT,"(' ',A)") BUFFER2(1:LEN_TRIM(BUFFER2)) 
            WRITE(IOUERR,"(' ',A)") BUFFER2(1:LEN_TRIM(BUFFER2)) 
            TOTERR = TOTERR + INERR 
ccs            EXIT MAINLOOP 
         ENDIF 
ccs         NUMDID = NUMDID + 1 
C 
         RC = 0 
         INERR = 0 
         DBNRUN = 0 
         DBNSCEN = 0 
         DO 
C 
            DBNRUN = DBNRUN + 1 
            DBNSCEN = 0 
            CALL MOBILE(INPFILE,DBNRUN,DBNFILE,INERR,RC) 
            TOTERR = TOTERR + INERR 
C 
            IF (RC < 0) EXIT 
C 
            WRITE(IOUOUT,140) DBNRUN,INERR 
  140       FORMAT(/' Run #',I3,', INERR =',I3/) 
            WRITE(IOUERR,140) DBNRUN,INERR 
C 
         END DO 
C 
         INQUIRE(IOUGEN,OPENED=FOPEN) 
         IF (FOPEN) CLOSE(IOUGEN) 
         INQUIRE(IOUREP,OPENED=FOPEN) 
         IF (FOPEN) CLOSE(IOUREP) 
         IF (DBFLAG == 2) THEN 
            DO IDX = 1,7 
               CLOSE(DBUNIT(IDX)) 
            END DO 
            DBFLAG = 1 
         ENDIF 
 
ccs      END DO MAINLOOP 
C 
      WRITE(IOUOUT,145) 
  145 FORMAT(/' DRIVER calls completed.') 
      IF (RC > 0) THEN 
         WRITE(BUFFER1,146) RC 
  146    FORMAT(' Program completion code is ',I6) 
         CALL LJCOMP(BUFFER1,BUFFER2) 
         WRITE(IOUOUT,"(' ',A)") BUFFER2(1:LEN_TRIM(BUFFER2)) 
      ENDIF 
C 
      CALL DATE_AND_TIME(CDATE,CENDTIME,CZONE,IENDTIME) 
      CDATE=CDATE(3:8) 
      CZONE=CZONE(1:5) 
      WRITE(IOUOUT,150) CENDTIME(1:2),CENDTIME(3:4),CENDTIME(5:10) 
  150 FORMAT(/' Processing end time is ',A2,':',A2,':',A6,'.') 
 
      CALL RUNTIMES(IENDTIME,ISTARTTIME,ITOTAL_TIME) 
      WRITE(IOUOUT,155) ITOTAL_TIME(1), ITOTAL_TIME(2), ITOTAL_TIME(3) 
  155 FORMAT(' The total run time is ',I2,' hrs., ',I2, 
     &       ' minutes, and ',I2,' seconds.') 
C 
ccs      IF (BATMODE) THEN 
ccs          WRITE(IOUOUT,170) DBNFILE 
ccs          WRITE(IOUOUT,175) NUMDID 
ccs  170     FORMAT(/' ',I4,' batch file(s) read') 
ccs  175     FORMAT(' ',I4,' batch file(s) processed') 
ccs      ENDIF 
C 
ccs      IF (IOUOUT.EQ.9) THEN 
ccs         WRITE(IOUOUT,180) 
ccs  180    FORMAT(/' Processing complete.', 
ccs     &           ' Hit Return to exit from the program.') 
ccs         READ(IOUIN,185) 
ccs  185    FORMAT() 
ccs      ENDIF 
C 
      GOTO 99 
C 
   80 CONTINUE 
      WRITE(BUFFER1,580) RC,INPFILE 
  580 FORMAT(' *** ERROR ',I9,' opening file ',A80) 
      CALL LJCOMP(BUFFER1, BUFFER2) 
      WRITE(IOUOUT,'(A150)') BUFFER2 
      WRITE(IOUERR,'(A150)') BUFFER2 
      INERR=INERR+1 
      TOTERR=TOTERR+INERR
      GOTO 99 
C 
   85 CONTINUE 
      WRITE(BUFFER1,585) RC,INPFILE 
  585 FORMAT(' *** ERROR ',I9,' reading file ',A80) 
      CALL LJCOMP(BUFFER1, BUFFER2) 
      WRITE(IOUOUT,'(A150)') BUFFER2 
      WRITE(IOUERR,'(A150)') BUFFER2 
      TOTERR=TOTERR+INERR
      INERR=INERR+1 
      GOTO 99 
C 
   90 WRITE(IOUOUT,590) INPFILE 
  590 FORMAT(/' *** The text MOBILE6 INPUT FILE or MOBILE6 BATCH', 
     &        ' FILE must be present in', 
     &       /' *** the label field of the first data record of', 
     &        ' every Mobile6 input file.' 
     &       /' *** (The label field is the first 19 columns of', 
     &        ' the record. The label may', 
     &       /' *** be entered in upper or lower case.) This', 
     &        ' record was not found in the', 
     &       /' *** following input file: ',A80) 
C     GOTO 99 
C 
   99 IF(TOTERR.EQ.0) THEN 
        EFLAG = .FALSE.
        M6FILE = TRIM( M6DIR ) // '/' // 'M6ERROR.TXT'
        INQUIRE(FILE=M6FILE,EXIST=FEXIST)  ! Ensure file exists before trying to delete it 
        IF(FEXIST) CLOSE (IOUERR,STATUS='DELETE') 
      ELSE 
      	EFLAG = .TRUE.
        WRITE(IOUOUT,595) TOTERR 
        INQUIRE(IOUREP,OPENED=FOPEN)
        IF(FOPEN) WRITE(IOUREP,595) TOTERR 
  595 FORMAT(// 
     *       ' ATTENTION: At least ',I3,' error(s) occurred during'/ 
     *       '            this MOBILE6 run. Please review the'/ 
     *       '            messages in the file M6ERROR.TXT,'/ 
     *       '            correct the problems, and try again.') 
      ENDIF 
C 
ccs   Make sure any files that are open get closed
      INQUIRE( IOUGEN, OPENED=FOPEN )
      IF( FOPEN ) CLOSE( IOUGEN )
      INQUIRE( IOUXIM, OPENED=FOPEN )
      IF( FOPEN ) CLOSE( IOUXIM )
      INQUIRE( IOUIMD, OPENED=FOPEN )
      IF( FOPEN ) CLOSE( IOUIMD )
      INQUIRE( IOUALT, OPENED=FOPEN )
      IF( FOPEN ) CLOSE( IOUALT )
      INQUIRE( IOUBAT, OPENED=FOPEN ) 
      IF( FOPEN ) CLOSE( IOUBAT ) 
      INQUIRE( IOUREP, OPENED=FOPEN )
      IF( FOPEN ) CLOSE( IOUREP )
      INQUIRE( IOUERR, OPENED=FOPEN )
      IF( FOPEN ) CLOSE( IOUERR )
      INQUIRE( IOUPM, OPENED=FOPEN )
      IF( FOPEN ) CLOSE( IOUPM )
      INQUIRE( IOUAT, OPENED=FOPEN )
      IF( FOPEN ) CLOSE( IOUAT )
      INQUIRE( IOUSS, OPENED=FOPEN )
      IF( FOPEN ) CLOSE( IOUSS )

      RETURN 
C 
C # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
C 
      CONTAINS 
C 
      SUBROUTINE RUNTIMES(ENDTIME,STARTTIME,RUNTIME) 
 
      INTEGER, DIMENSION(8), INTENT(IN)  ::  ENDTIME 
      INTEGER, DIMENSION(8), INTENT(IN)  ::  STARTTIME 
      INTEGER, DIMENSION(3), INTENT(OUT) ::  RUNTIME 
      INTEGER :: NUMSECS1 
      INTEGER :: TOTALTIME 
      INTEGER :: ENDHOUR 
 
      ENDHOUR=ENDTIME(5) 
      IF(ENDTIME(5) .LT. STARTTIME(5)) ENDHOUR=ENDTIME(5)+24 
 
      TOTALTIME=(ENDHOUR-STARTTIME(5))*3600+ 
     *          (ENDTIME(6)-STARTTIME(6))*60+ 
     *           ENDTIME(7)-STARTTIME(7) 
     * 
      RUNTIME(1) = INT(TOTALTIME/3600) 
      NUMSECS1   = MOD(TOTALTIME,3600) 
      RUNTIME(2) = INT(NUMSECS1/60) 
      RUNTIME(3) = MOD(NUMSECS1,60) 
C 
      RETURN 
      END SUBROUTINE RUNTIMES 
 
      END SUBROUTINE SMKDRIVER 
