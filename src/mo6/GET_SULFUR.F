      SUBROUTINE GET_SULFUR(INPREC,INERR,RC) 
C 
C    GET_SULFUR reads in the user supplied gasoline sulfur 
C    content in ppm (variable S_PPM) and if S_PPM is less than 
C    30 ppm GET_SULFUR also requires the User Supplied Sulfur 
C    Correction Factors (array USSCF(3,4,2,2) ). 
C 
C  Called by DSPTCH. 
C 
C  Calls CHKLAB,  QUITER. 
C 
C 
C  Changes: (Last change first) 
C 
C 10 Apr 11 AIR Task 11: Redirected error messages to IOUERR file. 
C 16 Oct 00 @EPA-bag bug277 Open external data files with READ only 
C  08 Jun 00: Removed type declarations for intrinsic functions. 
C  11 Apr 00 AIR Task 01: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C            so that Lahey Standard for keyboard/monitor I/O is followed. 
C  21 May 99 @DynTel-ddj 2-699, New Subroutine 
C 
C 
C  Input on call: 
C 
C    argument list: INERR 
C    common blocks: 
C    /FLAGS2/ IMFLAG 
C    /IOUCOM/ IOUALT,IOUOUT,IOUGEN 
C 
C  Output on return: 
C 
C    argument list: INERR 
C    common blocks: 
C    /FLAGS1/ TAMFLG 
C    /SULFUR/ S_PPM, USSCF 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C SUL_FLAG  I    Flag indicating source of sulfur correction factors: 
C                    1. Database regression coeficients. 
C                    2. User supplied sulfur correction factors. 
C S_PPM     R    Gasoline sulfur content (ppm). 
C 
C 
C  Notes: 
C 
      IMPLICIT NONE 
c 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'SULFUR.I' 
C  
      INTEGER, EXTERNAL :: JUNIT    
C 
      INTEGER     INERR, RC, SOURCE, I, NCHAR 
C 
      CHARACTER   FILENAME*80 
      CHARACTER   INPREC*(*) 
      CHARACTER   LASTREC*150 
       
ccs   Add declarations of M6DIR and M6FILE 
      CHARACTER(LEN=50) M6DIR 
      CHARACTER(LEN=80) M6FILE 
C 
C     Find out whether we are reading from an external data file. 
C 
      READ(INPREC(21:),*,IOSTAT=RC,END = 10,ERR=97) 
     &        S_PPM, FILENAME 
C 
C       Open the external data file. 
C 
  10  IF ( RC.EQ.0 ) THEN 
         SUL_FLAG = 2 
         SOURCE = JUNIT() 
 
ccs      All M6 input files live in SMK_M6PATH 
         CALL GETENV( "SMK_M6PATH", M6DIR ) 
         M6FILE = M6DIR( 1:LEN_TRIM( M6DIR ) ) // '/' // FILENAME 
         
         OPEN(SOURCE,FILE=M6FILE,STATUS='OLD',ACTION='READ', 
     *        IOSTAT=RC,ERR=80) 
C 
         READ(SOURCE,*,ERR=98,END=98,IOSTAT=RC) USSCF 
         CLOSE(SOURCE) 
         WRITE(IOUREP,300) S_PPM 
  300    FORMAT(10X,'***** User Supplied Sulfur Correction Factors used' 
     &      ,/, 10X,'***** in place of regression coefficients.' 
     &      ,//,10x,'User supplied gasoline sulfur content = ',F5.1, 
     &              ' ppm.',/) 
      ELSE 
         WRITE(IOUREP,301) S_PPM 
  301    FORMAT(/,10x,'User supplied gasoline sulfur content = ',F5.1, 
     &              ' ppm.',/) 
         IF (S_PPM.GE.30) THEN 
C 
C     Check ppm input to verify it is within upper limit of 600 ppm. 
C 
            RC = 0 
            IF(S_PPM.GT.600) THEN 
              WRITE(IOUREP,302) S_PPM 
  302         FORMAT(/,10X,'Error occurred in reading user supplied', 
     &                     ' gasoline sulfur content.'/, 
     &                 10x,'Sulfur content of ',F6.1, 
     &                     ' is greater than upper allowed', 
     &                     ' bound of 600 ppm.',/, 
     &                 10X,'Sulfur content set equal to 600 ppm',/) 
               S_PPM = 600. 
            END IF 
            SUL_FLAG = 1 
C 
C     Input sulfur level less than 30 ppm.  Check whether external 
C     file containing User Supplied Sulfur Correction Factors (USSCF) 
C     has been provided in input. 
C 
         ELSE 
            RC = 0 
            WRITE(IOUREP,303) S_PPM 
  303       FORMAT(/,10x,'User supplied gasoline sulfur content of ', 
     &              F5.1,' ppm requires an',/,10x, 
     &                   ' external file but none was provided.',/,10x, 
     &                   'Sulfur content set to 30 ppm.',/) 
            S_PPM = 30 
         END IF 
      END IF 
C 
C     Error opening the external data file. 
C 
      IF (RC.NE.0) THEN 
   80    WRITE(IOUERR,180) FILENAME 
  180    FORMAT(/'*** ERROR: unable to open the following external ', 
     &        'data file:',/'*** ',A) 
   97    RC = MOD(RC,256) 
         IF (RC.EQ.71) THEN 
           WRITE(IOUERR,181) 
  181      FORMAT('*** (The file name or path could not be found.)') 
           DO I = LEN(FILENAME),1,-1 
             NCHAR = I 
             IF (FILENAME(I:I).NE.' ') EXIT 
           END DO 
           IF (NCHAR.GT.50) THEN 
             WRITE(IOUERR,182) 
  182        FORMAT('*** (The file name must not be longer than 50 ', 
     &        'characters.)') 
             INERR = INERR + 1 
           END IF 
         ELSE 
   98      IF (RC.NE.0) THEN 
              CALL QUITER(0.,0,177,INERR) 
              IF (RC.LT.0) THEN 
                 WRITE(IOUERR,198) 
  198            FORMAT('*** (Unexpected end of file.)') 
              ELSE 
                 BACKSPACE(SOURCE) 
                 READ(SOURCE,199,ERR=99,END=99) LASTREC 
  199            FORMAT(A) 
                 WRITE(IOUERR,200) LASTREC 
  200            FORMAT('*** The following record caused the error:', 
     &              /'*** ',A) 
              END IF 
           END IF 
         ENDIF 
      END IF 
C 
   99 RETURN 
      END 
