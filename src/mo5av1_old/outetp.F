C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTETP
C
C  OUTETP prints the user supplied Evap Test Procedure parameters.
C
C  Called by OUTHD3 and OUTHD4.
C
C  Input on call:
C
C    common blocks:
C    /ETPPAR/ ETP,NEWETP,IETPYR,ETPPIV
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
C
      CHARACTER*16 ETPNAM(3)
C
      COMMON /ETPPAR/ ETP(7,3),NEWETP,IETPYR(4),ETPPIV(4),NETPYR
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
C
      DATA ETPNAM/'Problem Free    ',
     *            'Purge Failure   ',
     *            'Pressure Failure'/
C
      WRITE(IOUREP,200) (IETPYR(IETPY),IETPY=1,NETPYR)
  200 FORMAT('0','Evaporative Test Procedure Phase-in ',
     *           'Years and Percentages:'/' ','Model Year:',4I5)
      WRITE(IOUREP,205) (ETPPIV(IETPY),IETPY=1,NETPYR)
  205 FORMAT(' ','Percentage:',4F5.0)
C
      IF(NEWETP.EQ.2) THEN
        WRITE(IOUREP,210)
  210   FORMAT('0',18X,'Evaportive Test Procedure Reductions (%)'/
     *  '0',39X,'Diurnals',21X,'Losses'/
     *  ' ',26X,35('-'),1X,17('-')/
     *  ' ',17X,'Hot Soak   Full   Mult.2-3',
     *  '  Mult.4+  Partial  Running  Resting'/
     *  ' ',16X,7(' ',8('-')))
C
        DO 10 KEYEQ=1,3
        WRITE(IOUREP,220) ETPNAM(KEYEQ),(ETP(IETP,KEYEQ),IETP=1,7)
  220   FORMAT(' ',A16,3X,7(F4.0,5X))
   10   CONTINUE
      ENDIF
C
      RETURN
      END
