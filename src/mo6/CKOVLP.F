      SUBROUTINE CKOVLP(INERR,RC) 
C 
C  CKOVLP checks for overlaps between I/M programs. In Mobile6, it is 
C  illegal to have two I/M programs in operation in the same calendar year 
C  that cover the same model year for a given vehicle type. 
C 
C  Called by TESTIM 
C 
C  Changes: (Last change first) 
C 
C 31 Aug 00 AIR Task 08: Modified code to independently check for exhaust 
C                        and evap I/M program overlaps. 
C 10 Feb 99 @DynTel-MLA 2-663  This is a new routine for Mobile6. 
C 
C 
C  Input on call: 
C 
C    common blocks: 
C    /CHKLST/ IMCHECK,MAXIMPGM 
C    /FLAGS2/ IMFLAG,EVIMFLAG 
C    /IMPAR1/ ICYIM,MODYR 
C    /IMPAR2/ ILDT,ITEST 
C    /OIUCOM/ IOUERR 
C    /IVTYPE/ MAXVEH 
C    /VVNAMS/ VVGAS     (IVTYPE.I) 
C    /VVTYPS/ VTGASBUS  (IVTYPE.I) 
C 
C  Output on return: 
C    argument list: INERR, RC 
C 
C    Function: None. (Error messages are written to the error file.) 
C 
C 
C  Notes: 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'CHKLST.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'IMPAR1.I' 
      INCLUDE 'IMPAR2.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
C 
      INTEGER INERR 
      INTEGER RC 
C 
      INTEGER IV 
      INTEGER IVTAM 
      LOGICAL OVERLAP 
      INTEGER PGM1 
      INTEGER PGM2 
C 
C     If IMFLAG is less than or equal to two, no more than one exh I/M 
C     program is defined and there is no possibility of overlap. 
C 
      IF (IMFLAG.LE.2) GOTO 10 
C 
C     The following code checks for exhaust I/M overlap only. 
C 
      DO PGM1 = 1, MAXIMPGM - 1 
C 
C       If no data was entered for this program, skip it. 
C 
        IF (IMCHECK(1,PGM1).EQ.1) CYCLE 
C 
C       Since exhaust and evap I/M programs can appear to overlap, do not 
C       include evap programs in this loop. 
C 
        IF(ITEST(PGM1).GE.20) CYCLE        
C 
C       Check this program against the remaining programs. 
C 
        DO PGM2 = PGM1+1, MAXIMPGM 
C 
C         If no data was entered for this program, skip it. 
C 
          IF (IMCHECK(1,PGM2).EQ.1) CYCLE 
C 
C         Since exhaust and evap I/M programs can appear to overlap, do not 
C         include evap programs in this loop. 
C 
          IF(ITEST(PGM2).GE.20) CYCLE        
C 
C         In order to conflict, the programs must be in operation during 
C         the same calendar year, the model year coverage must overlap, 
C         and they must cover the at least one vehicle type in common. 
C 
          IF (ICYIM(1,PGM1).LE.ICYIM(1,PGM2)) THEN 
C           Program 1 starts first. Overlap occurs whenever program one's 
C           end date is greater than or equal to program two's start date. 
            OVERLAP = (ICYIM(2,PGM1).GE.ICYIM(1,PGM2)) 
          ELSE 
C           Program 2 starts first. Overlap occurs whenever program two's 
C           end date is greater than or equal to program one's start date. 
            OVERLAP = (ICYIM(2,PGM2).GE.ICYIM(1,PGM1)) 
          ENDIF 
C 
          IF (.NOT.OVERLAP) CYCLE 
C 
          IF (MODYR(1,PGM1).LE.MODYR(1,PGM2)) THEN 
C           Program 1 covers earlier model years. Overlap occurs whenever 
C           the last model year covered by program one is greater than 
C           the first model year covered by program two. 
            OVERLAP = (MODYR(2,PGM1).GE.MODYR(1,PGM2)) 
          ELSE 
C           Program 2 covers earlier model years. Overlap occurs whenever 
C           the last model year covered by program two is greater than 
C           the first model year covered by program one. 
            OVERLAP = (MODYR(2,PGM2).GE.MODYR(1,PGM1)) 
          ENDIF 
C 
          IF (.NOT.OVERLAP) CYCLE 
C 
C         Now check to see whether the programs have any vehicle types 
C         in common. 
C 
          OVERLAP = .FALSE. 
          DO IV = 1, MAXVEH 
C 
            IF(VVGAS(IV).EQ.0) CYCLE 
            IVTAM = IV 
            IF (IV.EQ.VTGASBUS) THEN 
              IVTAM = 14 
            ENDIF 
C 
            IF (ILDT(IVTAM,PGM1).EQ.2.AND. 
     *          ILDT(IVTAM,PGM2).EQ.2) THEN 
              OVERLAP = .TRUE. 
              EXIT 
            ENDIF 
C 
          END DO 
C         End of IV loop. 
C 
          IF (OVERLAP) THEN 
            WRITE(IOUERR,120) PGM1, PGM2 
  120       FORMAT('*** ERROR: Exhaust I/M programs ',I1,' and ', 
     *             I1,' overlap. They operate during', 
     *            /12X,'the same calendar year and they both cover ', 
     *             'a common model year', 
     *            /12X,'range for at least one vehicle type.') 
            RC = 1 
            INERR = INERR + 1 
          END IF 
C         
        END DO 
C       End of loop over PGM2. 
      END DO 
C     End of loop over PGM1. 
C 
C     If EVIMFLAG is less than or equal to two, no more than one evap I/M 
C     program is defined and there is no possibility of overlap. 
C 
   10 IF (EVIMFLAG.LE.2) GOTO 99 
C 
C     The following code checks for evap I/M overlap only.  
C 
      DO PGM1 = 1, MAXIMPGM - 1 
C 
C       If no data was entered for this program, skip it. 
C 
        IF (IMCHECK(1,PGM1).EQ.1) CYCLE 
C 
C       Since exhaust and evap I/M programs can appear to overlap, do not 
C       include exhaust programs in this loop. 
C 
        IF(ITEST(PGM1).LT.20) CYCLE        
C 
C       Check this program against the remaining programs. 
C 
        DO PGM2 = PGM1+1, MAXIMPGM 
C 
C         If no data was entered for this program, skip it. 
C 
          IF (IMCHECK(1,PGM2).EQ.1) CYCLE 
C 
C         Since exhaust and evap I/M programs can appear to overlap, do not 
C         include exhaust programs in this loop. 
C 
          IF(ITEST(PGM2).LT.20) CYCLE        
C 
C         In order to conflict, the programs must be in operation during 
C         the same calendar year, the model year coverage must overlap, 
C         and they must cover the at least one vehicle type in common. 
C 
          IF (ICYIM(1,PGM1).LE.ICYIM(1,PGM2)) THEN 
C           Program 1 starts first. Overlap occurs whenever program one's 
C           end date is greater than or equal to program two's start date. 
            OVERLAP = (ICYIM(2,PGM1).GE.ICYIM(1,PGM2)) 
          ELSE 
C           Program 2 starts first. Overlap occurs whenever program two's 
C           end date is greater than or equal to program one's start date. 
            OVERLAP = (ICYIM(2,PGM2).GE.ICYIM(1,PGM1)) 
          ENDIF 
C 
          IF (.NOT.OVERLAP) CYCLE 
C 
          IF (MODYR(1,PGM1).LE.MODYR(1,PGM2)) THEN 
C           Program 1 covers earlier model years. Overlap occurs whenever 
C           the last model year covered by program one is greater than 
C           the first model year covered by program two. 
            OVERLAP = (MODYR(2,PGM1).GE.MODYR(1,PGM2)) 
          ELSE 
C           Program 2 covers earlier model years. Overlap occurs whenever 
C           the last model year covered by program two is greater than 
C           the first model year covered by program one. 
            OVERLAP = (MODYR(2,PGM2).GE.MODYR(1,PGM1)) 
            ENDIF 
C 
          IF (.NOT.OVERLAP) CYCLE 
C 
C         Now check to see whether the programs have any vehicle types 
C         in common. 
C 
          OVERLAP = .FALSE. 
          DO IV = 1, MAXVEH 
C 
            IF(VVGAS(IV).EQ.0) CYCLE 
            IVTAM = IV 
            IF (IV.EQ.VTGASBUS) THEN 
              IVTAM = 14 
            ENDIF 
C 
            IF (ILDT(IVTAM,PGM1).EQ.2.AND. 
     *          ILDT(IVTAM,PGM2).EQ.2) THEN 
              OVERLAP = .TRUE. 
              EXIT 
            ENDIF 
C 
          END DO 
C         End of IV loop. 
C 
          IF (OVERLAP) THEN 
            WRITE(IOUERR,122) PGM1, PGM2 
  122       FORMAT('*** ERROR: Evap I/M programs ',I1,' and ', 
     *             I1,' overlap. They operate during', 
     *            /12X,'the same calendar year and they both cover ', 
     *             'a common model year', 
     *            /12X,'range for at least one vehicle type.') 
            RC = 1 
            INERR = INERR + 1 
          END IF 
C         
        END DO 
C       End of loop over PGM2. 
      END DO 
C     End of loop over PGM1. 
C 
   99 CONTINUE 
      RETURN 
      END 
