      SUBROUTINE SCFCO() 
C 
C  SFCCO Calculate the Off-Cycle coefficient and level and 
C  the speed correction factor from emission versus 
C  speed equation by emission level 
C 
C 02 Oct 00 @EPA-djb Bug230, Removed Open/Write statement. 
C 17 Aug 99 @DynTel-ddj 2-698  Streamlined OCLEVEL, BEFST, SCFF7, SCFF, 
C           BEF30 and BEFSP calculations. 
C  1 March 1999 @DynTel-HQ 2-698 calculate SFC and Off-cycle coefficient 
C 
C  CALLED BY REINIT 
C 
C     common blocks: 
C     /SPEED9/ OCECO, OCY, OCX, OCLEVEL, EMSFC 
C                 
C 
      IMPLICIT NONE 
      INCLUDE 'SPEED9.I' 
C 
      INTEGER I,J,K  
      REAL SP, SPST, SP30, SP65, SPIN, A1, A2, 
     * BEFST(3,4), BEF30, BEF65, BEFSP, MBEF30, BEFFR 
C 
C  initialize the different speed 
C 
      SPST = 19.6 
      SP30 = 30.5 
      SP65 = 65.0 
      SPIN =  5.0 
C 
C  Here we calculate interpolation factors OCEF and OCLEVEL that are 
C  stored in common block SPEED9.  OCEF and OCLEVEL are used in BEFSPAD. 
C 
      DO I = 1, 4 
         A1=(OCY(2,I)-OCY(1,I))/(OCX(2,I)-OCX(1,I)) 
         OCEF(1,I)=(OCY(2,I)-OCX(2,I)*A1)/(1-A1) 
         OCEF(2,I)=A1/(1-A1) 
         A2=(OCY(3,I)-OCY(2,I))/(OCX(3,I)-OCX(2,I)) 
         OCEF(3,I)=(OCY(3,I)-OCX(3,I)*A2)/(1-A2) 
         OCEF(4,I)=A2/(1-A2) 
C 
C Determine the off cycle level 
C 
         DO J = 1, 3 
            OCLEVEL(J,I) = OCX(J,I) - OCY(J,I) 
         END DO 
      END DO 
C 
C  Calculate the Speed Correction Factors 
C 
C Freeway Road 
C Emission for speed 19.6 by g/mi 
      DO I=1,4 
        DO J=1,3 
          BEFST(J,I)=EMSCF(1,J,3,I)/SPST+EMSCF(2,J,3,I) 
C SPEED 7.1 
          SCFF7(J,I)=(EMSCF(1,J,2,I)/7.1+EMSCF(2,J,2,I))/BEFST(J,I) 
C SPEED 10 
          SCFF(J,3,I)=(EMSCF(1,J,2,I)/10.0+EMSCF(2,J,2,I))/BEFST(J,I) 
        END DO 
      END DO 
C  SPEED 15 20 25 30 
      SP=10.0 
      DO K=4,7 
        SP=SP+SPIN 
        DO I=1,4 
          DO J=1,3 
            SCFF(J,K,I)=(EMSCF(1,J,3,I)/SP+EMSCF(2,J,3,I))/BEFST(J,I) 
          END DO 
        END Do 
      END DO 
C  SPEED 35 40 45 50 55 60 65 
      SP=30.0 
      DO K=8,14 
        SP=SP+SPIN 
        DO I=1,4 
          DO J=1,3 
            BEF30=EMSCF(1,J,3,I)/SP30 + EMSCF(2,J,3,I) 
            BEF65=EMSCF(1,J,4,I)+EMSCF(2,J,4,I)*SP65 
            BEFSP=EMSCF(1,J,4,I)+EMSCF(2,J,4,I)*SP 
            IF(BEF30.LE.BEF65) THEN 
              SCFF(J,K,I)=MAX(BEF30,BEFSP)/BEFST(J,I) 
            ELSE 
              MBEF30=BEF30*SP30/SP 
              SCFF(J,K,I)=MAX(MBEF30,BEFSP)/BEFST(J,I) 
            END IF 
          END DO 
        END DO 
      END DO 
C 
C Arterial/Collector 
C 
C Speed 7.1 
C 
      SP=7.1 
      DO I=1,4 
        DO J=1,3 
          BEFSP=EMSCF(1,J,1,I)/SP + EMSCF(2,J,1,I) 
          BEFFR=SCFF7(J,I)*BEFST(J,I) 
          SCFAC7(J,I)=MAX(BEFSP,BEFFR)/BEFST(J,I) 
        END DO 
      END DO 
C 
C Speed 10 15 20 25 30 
C 
      SP=5.0 
      DO K=3,7 
        SP=SP+SPIN 
        DO I=1,4 
          DO J=1,3 
            BEFSP=EMSCF(1,J,1,I)/SP + EMSCF(2,J,1,I) 
            BEFFR=SCFF(J,K,I)*BEFST(J,I) 
            SCFAC(J,K,I)=MAX(BEFSP,BEFFR)/BEFST(J,I) 
          END DO 
        END DO 
      END DO 
C 
C  SPEED 35 40 45 50 55 60 65 
      SP=30.0 
      DO K=8,14 
        SP=SP+SPIN 
        DO I=1,4 
          DO J=1,3 
            MBEF30=(EMSCF(1,J,1,I)+EMSCF(2,J,1,I)*SP30)/SP 
            BEFFR=SCFF(J,K,I)*BEFST(J,I) 
            SCFAC(J,K,I)=MAX(MBEF30,BEFFR)/BEFST(J,I) 
          END DO 
        END DO 
      END DO 
C 
      RETURN 
      END 
