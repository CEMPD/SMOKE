      SUBROUTINE OUTETP 
C 
C  OUTETP prints the user supplied Evap Test Procedure parameters. 
C 
C  Called by OUTHD3. 
C 
C  Changes: (Last change first) 
C 
C  06 Oct 00 @EPA-BG  BUG155 removed NEWETP flag previously used for user input 
C            of ETP phase in ; removed code asociated with NEWETP and user supplied 
C            phase in 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C  03 Dec 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    common blocks: 
C    /ETPPAR/ ETP,IETPYR,ETPPIV 
C    /IOUCOM/ IOUREP 
C 
C  Local variable / array dictionary: 
C 
C  ETPNAM(3)  -  ETPNAM ( KEYEQ ) 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  ETPNAM  C    names of the ETP groups: Problem Free, Purge Failure, 
C               and Pressure Failure 
C 
C  Notes: 
C 
C  OUTETP was added for MOBILE5v1. 
C 
      IMPLICIT NONE      
      INCLUDE 'ETPPAR.I' 
      INCLUDE 'IOUCOM.I' 
C 
      INTEGER KEYEQ, IETPY, IETP 
      CHARACTER(LEN=16), DIMENSION(3) :: ETPNAM=(/ 
     *            'Problem Free    ', 
     *            'Purge Failure   ', 
     *            'Pressure Failure'/) 
C 
      WRITE(IOUREP,200) (IETPYR(IETPY),IETPY=1,NETPYR) 
  200 FORMAT('0','Evaporative Test Procedure Phase-in ', 
     *           'Years and Percentages:'/' ','Model Year:',4I5) 
      WRITE(IOUREP,205) (ETPPIV(IETPY),IETPY=1,NETPYR) 
  205 FORMAT(' ','Percentage:',4F5.0) 
C 
      DO KEYEQ=1,3 
        WRITE(IOUREP,220) ETPNAM(KEYEQ),(ETP(IETP,KEYEQ),IETP=1,7) 
  220   FORMAT(' ',A16,3X,7(F4.0,5X)) 
      ENDDO 
C 
      RETURN 
      END 
