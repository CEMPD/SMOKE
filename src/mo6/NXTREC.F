      SUBROUTINE NXTREC(SOURCE, REPORT, M6REC, INERR, RC) 
C 
C  This subroutine returns the next data record from the specified file. 
C  Comments are discarded or echoed to the output file, as appropriate. 
C 
C  Called by GETATP, FNDLAB. 
C 
C  Changes: (Last change first) 
C 
C  11 Apr 00 AIR Task 01: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C            so that Lahey Standard for keyboard/monitor I/O is followed. 
C  24 Aug 99 @DynTel-MLA 2-684  Don't echo trailing blanks on comments. 
C            Don't echo comments at all if the unit number for the report 
C            file is invalid (less than 0). Add a leading blank to comments 
C            that are echoed to the terminal. 
C  25 Nov 98 @DynTel-Zk 2-000 Explicit Typing 
C  Input on call: 
C 
C    SOURCE  - The I/O unit to read for input records. This may be the 
C              standard M6 input file or an external data file. 
C    REPORT  - The I/O unit connected to the output file. 
C 
C  Output on return: 
C 
C    M6REC   - The next data record from the SOURCE input file. 
C    RETCODE - The return code. It will have a value greater than 0 if an 
C              error was reported by a FORTRAN READ or WRITE operation, a 
C              value less than 0 if EOF was encountered before a valid  
C              data record was found, and 0 otherwise. See the NOTES for  
C              suggestions out coding a call to this routine from another  
C              subroutine. 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C 
C  Notes: 
C 
C  NXTREC was added for MOBILE6. 
C 
C  The call to this subroutine should be followed in the calling program  
C  by a CASE statment branching on RETCODE, as follows. Note that the 
C  CASE (:-1) handles all cases less than zero and the CASE (1:) handles 
C  all cases greater than zero. 
C      CALL NXTREC(SOURCE, REPORT, INPREC, INERR, RC) 
C      SELECT CASE (RETCODE) 
C        CASE (:-1) 
C          GOTO aaaa 
C        CASE (1:) 
C          GOTO bbbb 
C        CASE DEFAULT 
C      END SELECT 
C  This will have the same effect as a FORTRAN READ statment of the 
C  form READ(IOUNIT,FMT,IOSTAT=IOS,END=aaaa,ERR=bbbb)... bbbb RC=IOS. 
C 
      IMPLICIT NONE 
      INCLUDE 'IOUCOM.I' 
C 
      CHARACTER*(*)  M6REC 
      INTEGER        IOS 
      INTEGER        LT 
      INTEGER        RC 
      INTEGER        REPORT 
      INTEGER        SOURCE 
      INTEGER        INERR 
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
  100 RC = 0 
C 
      M6REC = ' ' 
   20 READ(SOURCE,120,IOSTAT=IOS,ERR=80,END=85) M6REC 
  120 FORMAT(A) 
  121 FORMAT(' ',A) 
C 
      IF (M6REC.EQ.' '.OR.M6REC(1:1).EQ.'*') THEN 
        GOTO 20 
      ENDIF 
C 
      IF (M6REC(1:1).EQ.'>') THEN 
        IF (REPORT.GE.0) THEN 
          M6REC(1:1) = '*' 
          LT = LEN_TRIM(M6REC) 
          IF (REPORT == IOUOUT) THEN 
            WRITE(REPORT,121,IOSTAT=IOS,ERR=80) M6REC(1:LT) 
          ELSE 
            WRITE(REPORT,120,IOSTAT=IOS,ERR=80) M6REC(1:LT) 
          ENDIF 
        ENDIF 
        GOTO 20 
      ENDIF 
C 
      GOTO 90 
C 
   80 RC = IOS 
      INERR = INERR + 1 
      GOTO 90 
C 
   85 RC = -1 
C 
   90 RETURN 
      END 
