      INTEGER FUNCTION FINDIM(ICY,MY,IV,IEXEV) 
C 
C  FINDIM returns the I/M program number for a given calendar year, model year, 
C  vehicle type case. The routine differentiates between exhaust and evap 
C  I/M requests 
C 
C  Called by IMPROG 
C 
C  Changes: (Last change first) 
C 
C 03 Dec 01 @EPA-bag        Fix I/M EXEMPTION AGE conditional for exhaust 
C                           and evaporative I&M 
C 30 Oct 00 @EPA-bag bug284 Fix I/M EXEMPTION AGE conditional from G.E. to G.T.; 
C           exemption age of 1 should exempt all vehicles age 1 or older( emissions  
C           should be the same as in the no I/M case.) 
C 15 Sep 00 AIR Task 08: Added switch to recognize need for either exhaust 
C                        or evap I/M program. Also fixed bug where I/M program 
C                        was selected 1 year too late. 
C 10 Aug 99 @DynTel-bg  2-663  Replaced constants with descriptive parameters. 
C  2 Jun 99 @DynTel-MLA 2-000  Account for grace periods and exemption ages. 
C 10 Feb 99 @DynTel-MLA 2-663  This is a new routine for Mobile6. 
C 
C  Input on call: 
C 
C    common blocks: 
C    /CHKLST/ IMCHECK,MAXIMPGM 
C    /IMPAR1/ EXEMPTAGE,ICYIM,IMGRACE,MODYR 
C    /IMPAR2/ ILDT 
C    /IVTYPE/ MAXVEH 
C    /VVNAMS/ VVGAS     (IVTYPE.I) 
C    /VVTYPS/ VTGASBUS  (IVTYPE.I) 
C 
C  Output on return: 
C    argument list: 
C    function FINDIM 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------ 
C  IMPGM    I    I/M program number (loop control variable) 
C  IVTAM    I    Veh. type (mapped into range of vehs. subject to I/M) 
C  IEXEV    I    Exhaust I/M = 1, Evap I/M = 2 
C 
C  Notes: 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'CHKLST.I' 
      INCLUDE 'IMPAR1.I' 
      INCLUDE 'IMPAR2.I' 
      INCLUDE 'IVTYPE.I' 
C 
      INTEGER, INTENT(IN) :: ICY 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: MY 
      INTEGER, INTENT(IN) :: IEXEV 
C 
      INTEGER IMPGM 
      INTEGER IVTAM 
      INTEGER, SAVE :: EXHAUSTIM=1,EVAPIM=2 
      INTEGER, SAVE :: BEGIN_DATE=1,END_DATE=2 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
      FINDIM = 0 
C 
      IF(VVGAS(IV).EQ.0) GOTO 99 
      IVTAM = IV 
      IF (IV.EQ.VTGASBUS) THEN 
        IVTAM = 14 
      ENDIF 
C 
      DO IMPGM = 1, MAXIMPGM 
C 
C       If no data was entered for this program, skip it. 
C 
        IF (IMCHECK(1,IMPGM).EQ.1) CYCLE 
C 
C       Do not allow exhaust requests to select evap programs and visa-versa 
C 
        IF(IEXEV==EXHAUSTIM.AND.ITEST(IMPGM)<20) THEN 
          IF (ILDT(IVTAM,IMPGM).EQ.2.AND. 
     *      ICY.GE.ICYIM(BEGIN_DATE,IMPGM).AND. 
     *      ICY.LE.ICYIM(END_DATE,IMPGM).AND. 
C 
     *      MY.GE.MODYR(BEGIN_DATE,IMPGM).AND. 
     *      MY.LE.MODYR(END_DATE,IMPGM).AND. 
C 
     *      ICY.GE.(MY+IMGRACE(IMPGM)).AND. 
C 
     *      ICY.LE.(MY+EXEMPTAGE(IMPGM)) ) THEN 
C 
            FINDIM = IMPGM 
            GOTO 99 
          ENDIF 
        END IF 
C 
        IF(IEXEV==EVAPIM.AND.ITEST(IMPGM)>=20) THEN 
          IF (ILDT(IVTAM,IMPGM).EQ.2.AND. 
     *      ICY.GE.ICYIM(BEGIN_DATE,IMPGM).AND. 
     *      ICY.LE.ICYIM(END_DATE,IMPGM).AND. 
C 
     *      MY.GE.MODYR(BEGIN_DATE,IMPGM).AND. 
     *      MY.LE.MODYR(END_DATE,IMPGM).AND. 
C 
     *      ICY.GT.(MY+IMGRACE(IMPGM)).AND. 
C 
     *      ICY.LE.(MY+EXEMPTAGE(IMPGM)) ) THEN 
C 
            FINDIM = IMPGM 
            GOTO 99 
          END IF 
        ENDIF 
C 
      END DO 
C     End of loop over IMPGM 
C 
   99 CONTINUE 
      RETURN 
      END 
