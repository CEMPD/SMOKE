      SUBROUTINE FNDLAB(M6REC,STDLAB,LABTYPE,SECTION,INERR,RC) 
C 
C  This subroutine returns the next labeled input record from the source 
C  file. 
C 
C  Called by MOBILE. 
C 
C  Changes: (Last change first) 
C 
C  28 Jul 99 @DynTel-MLA 2-000 Removed the USE LABELS statement. The 
C            module subroutine CHKLAB was removed from that module 
C            and converted to an ordinary external procedure. 
C  11 May 99 @DynTel-MLA 2-674  Modified calls to CHKLAB because the 
C            calling sequence for that subroutine was changed. 
C   3 Mar 99 @DynTel-MLA 2-659  Changed IOUREP to IOUERR in error message. 
C  02 Oct 1998 @DynTel-ZK 2-000 explicit Typing 
C  Input on call: 
C 
C    /IOUCOM/  IOUGEN,IOUREP 
C 
C  Output on return: 
C 
C    LABTYPE - The token number for the label on the current input rec. 
C    M6REC   - The next data record from the input file. 
C    SECTION - The input section where the label belongs (1=global data 
C              only, 2=global or scenario, 3=scenario section only) 
C    STDLAB  - The full label text (unabbreviated, uppercase, left justified) 
C    RC      - The return code. It will have a value greater than 0 if an 
C              error was reported by a FORTRAN READ or WRITE operation, a 
C              value less than 0 if EOF was encountered before a valid  
C              data record was found, and 0 otherwise. 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C 
C  Notes: 
C 
C  NXTLAB was added for MOBILE6. 
C 
      IMPLICIT NONE 
      INCLUDE 'IOUCOM.I' 
C 
      INTEGER        INERR 
      CHARACTER      LABEL*19 
      INTEGER        LABTYPE 
      CHARACTER      M6REC*(*) 
      INTEGER        RC 
      INTEGER        SECTION 
      CHARACTER      STDLAB*40 
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
   10 RC = 0 
C 
   20 CALL NXTREC(IOUGEN,IOUREP,M6REC,INERR,RC) 
C 
      SELECT CASE (RC) 
        CASE (:-1) 
          M6REC = 'END OF RUN         : 000000000000000000' 
        CASE (1:) 
          INERR = INERR + 1 
          GOTO 90 
        CASE DEFAULT 
      END SELECT 
C 
C     Check the label type. If the label type is unrecognized, print 
C     out an error message and continue with the next input record. 
C 
      CALL UCCOMP(M6REC(1:19),LABEL) 
   25 CALL CHKLAB(LABEL,STDLAB,LABTYPE,SECTION) 
C 
      IF (LABTYPE.EQ.0) THEN 
        WRITE(IOUERR,130) M6REC 
  130   FORMAT('*** ERROR: The following record does not contain', 
     &         ' a valid input label:',/A) 
        INERR = INERR + 1 
        GOTO 20 
      ENDIF 
C 
      IF (RC.LT.0) THEN 
        WRITE(IOUERR,140) 
  140   FORMAT('*** WARNING: A missing END OF RUN record was', 
     &         ' supplied.') 
      ENDIF 
C 
   90 RETURN 
      END 
