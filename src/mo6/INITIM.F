      SUBROUTINE INITIM() 
C 
C  Calculates arrays used in the computation of I/M benefits for 1981 to 
C  1993 model year light duty gasoline vehicles and trucks. INITIM is 
C  called before the evaluation of the basic emission factors for each 
C  scenario. 
C 
C  Called by CHKINP. 
C 
C  Calls 
C 
C  Changes : 
C  03 Dec 01 @EPA-djb Bug342 Cap IM240 IDR rates. 
C  03 Dec 01 @EPA-BAG        Remove AGEIM, VMTIM, and associated variables/ see 
C                            also IMPROG and FINDIM 
C  17 Oct 01 @EPA-BAG bug331 Biennial/ annual I&M program discrepancy 
C  02 Nov 00 @EPA-djb Bug299, initialized AGEIM array. 
C  15 Sep 00 AIR Task 08: Expanded I/M arrays to MAXIMPGM possible programs. 
C            Added code to process exhaust and evap I/M programs separately. 
C  07 Jul 00 @EPA-elg    All ASM tests now have final IM240 cutpoints as a basis. 
C  04 Aug 99 @DynTel-bg  2-663 Added pollutant, emissions, and year 
C            partition type parameters; comments. 
C  23 Nov 98 @DynTel-RJD 2-663 Created INITIM to calculate arrays needed 
C            for the evaluation of the 81-93 model year light duty 
C            gasoline vehicle and truck I/M benefits. 
C 
C  Input on call: 
C 
C    argument list: 
C    common blocks: 
C    /IMPAR1/ IMGRACE, EXEMPTAGE 
C    /IMPAR4/ CUTHC, CUTCO, CUTNO 
C    /IMPAR6/ IFREQ 
C    /LDGIMC/ HCCUT,COCUT,NOCUT 
C    /MAXIMA/ MAXYRS 
C    /MYRSAV/ AMAR 
C    /RSDINF/ RSDDIRTY,RSDCLEAN 
C 
C  Output on return: 
C 
C    common blocks: 
C    /IMPAR1/ OLDFRAC, NEWFRAC 
C    /IMPAR2/ ITEST 
C    /LDGIMC/ IDR240 
C    /LDGSIM/ RSDCVR1 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  EPART    R    Fraction of vehicles that RSD does not cover because of 
C                an insufficient number of RSD tests 
C  IAIM     I    The integer part age of the vehicle in years on the 
C                date of its previous IM test. 
C  IDR240   R    Identification rate for high emitters as a result of an 
C                IM240 test. 
C  ITN      I    Counter for the number of RSD tests. 
C  IMTT     I    I/M test type nunmber 
C  IMVG     I    I/M vehicle group number 
C  IY       I    IY is reset to 1 immediately following any I/M test. It 
C                is then incremented by 1 on each of the anniversaries 
C                of that test. It counts the number of years between tests. 
C  KY       I    Counter for the number of periodic IM tests taken by a 
C                vehicle 
C  LAMBDA   R    The total number of RSD tests during a given year 
C                expressed as a fraction of the total number of vehicles 
C                in the fleet and weighted by the mileage accumulation 
C                rate. 
C  NFCTRL   I    Factorial: NFCTRL(X)=X! 
C 
      IMPLICIT NONE 
      INCLUDE 'CHKLST.I' 
      INCLUDE 'EVAL.I' 
      INCLUDE 'IMPAR1.I' 
      INCLUDE 'IMPAR2.I' 
      INCLUDE 'IMPAR4.I' 
      INCLUDE 'IMPAR6.I' 
      INCLUDE 'LDGIMC.I' 
      INCLUDE 'LDGSIM.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'MYRSAV.I' 
      INCLUDE 'RSDINF.I' 
C 
      REAL    LAMBDA,EPART 
      REAL, DIMENSION(26,5) :: VMTIM 
      INTEGER JDX,IV,IAIM,IMPGM,IMTT,ITN,NFCTRL,IMVG 
      INTEGER, SAVE :: RUN=1,START=2 
      INTEGER, SAVE :: HC=1,CO=2,NOX=3 
C 
      VMTIM = 0.0 
C 
C Set the final cutpoints for ASM tests (See page 57 
C 2nd paragraph of M6.IM.001) 
C Only peform this for exhaust I/M programs 
C 
      DO IMPGM=1,MAXIMPGM 
        IF(ITEST(IMPGM).GE.20) CYCLE 
        DO IMVG=1,4 
          DO JDX=1,MAXYRS 
            IMTT = ITEST(IMPGM) 
            IF (IMTT.GE.6 .AND. IMTT.LE.11) THEN 
              CUTHC(JDX,IMVG,IMPGM) = 0.8 
              CUTCO(JDX,IMVG,IMPGM) = 15.0 
              CUTNO(JDX,IMVG,IMPGM) = 2.0 
            END IF 
          END DO 
        END DO 
      END DO 
C 
C Calculate the high emitter identification rates IDR for IM240 
C tests from the IM cutpoints (See pages 25-27 for running emissions 
C and page 54 for start emissions of M6.IM.001) 
C Only peform this for exhaust I/M programs 
C 
      DO IMPGM=1,MAXIMPGM 
        IF(ITEST(IMPGM).GE.20) CYCLE 
        DO IMVG=1,4 
          DO JDX=1,MAXYRS 
C 
C           Cap HC and CO cutpoints. 
            IF(CUTHC(JDX,IMVG,IMPGM).GT.5.0) THEN 
              CUTHC(JDX,IMVG,IMPGM)=5.0 
            ENDIF 
            IF(CUTCO(JDX,IMVG,IMPGM).GT.100.0) THEN 
              CUTCO(JDX,IMVG,IMPGM)=100.0 
            ENDIF 
C 
C           IDR for HC running emissions 
            IDR240(JDX,RUN,HC,IMVG,IMPGM) 
     *          =1.1451 - 0.1365*LOG(CUTHC(JDX,IMVG,IMPGM)) 
     *                  - 0.1069*LOG(CUTCO(JDX,IMVG,IMPGM)) 
            IF(IDR240(JDX,RUN,HC,IMVG,IMPGM).LT.0.0) THEN 
              IDR240(JDX,RUN,HC,IMVG,IMPGM)=0.0             !Set Minimum 
            ENDIF 
C 
C           IDR for CO running emissions 
            IDR240(JDX,RUN,CO,IMVG,IMPGM) 
     *          = 1.1880 - 0.1073*LOG(CUTHC(JDX,IMVG,IMPGM)) 
     *                   - 0.1298*LOG(CUTCO(JDX,IMVG,IMPGM)) 
            IF(IDR240(JDX,RUN,CO,IMVG,IMPGM).LT.0.0) THEN 
              IDR240(JDX,RUN,CO,IMVG,IMPGM)=0.0             !Set Minimum 
            ENDIF 
C 
C           IDR for NOx running emissions 
            IDR240(JDX,RUN,NOX,IMVG,IMPGM) 
     *          = 0.5453 + 0.7568*CUTNO(JDX,IMVG,IMPGM) 
     *                   - 0.3687*CUTNO(JDX,IMVG,IMPGM)**2 
     *                   + 0.0406*CUTNO(JDX,IMVG,IMPGM)**3 
            IF(CUTNO(JDX,IMVG,IMPGM).GT.4.5) THEN 
              IDR240(JDX,RUN,NOX,IMVG,IMPGM)=0.0             !Set Minimum 
            ENDIF 
C 
C           IDR for HC start emissions 
            IDR240(JDX,START,HC,IMVG,IMPGM) 
     *          = 0.9814 - 0.1590*LOG(CUTHC(JDX,IMVG,IMPGM)) 
     *                   - 0.1409*LOG(CUTCO(JDX,IMVG,IMPGM)) 
            IF(IDR240(JDX,START,HC,IMVG,IMPGM).LT.0.0) THEN 
              IDR240(JDX,START,HC,IMVG,IMPGM)=0.0           !Set Minimum 
            ENDIF 
C 
C           IDR for CO start emissions 
            IDR240(JDX,START,CO,IMVG,IMPGM) 
     *          = 1.1460 - 0.1593*LOG(CUTHC(JDX,IMVG,IMPGM)) 
     *                   - 0.1707*LOG(CUTCO(JDX,IMVG,IMPGM)) 
            IF(IDR240(JDX,START,CO,IMVG,IMPGM).LT.0.0) THEN 
              IDR240(JDX,START,CO,IMVG,IMPGM)=0.0           !Set Minimum 
            ENDIF 
C 
C           IDR for NOx start emissions 
            IDR240(JDX,START,NOX,IMVG,IMPGM) = 0.0 
C 
          END DO 
        END DO 
      END DO 
C 
C 
C Calculate the October 1 accummulated vehicle miles traveled. 
C AMAR(JDX,IV) 
C 
      DO IV=1,5 
        VMTIM(1,IV)=10.0*AMAR(1,IV) 
        DO JDX=2,MAXYRS 
          VMTIM(JDX,IV)=VMTIM(JDX-1,IV)+10.0*AMAR(JDX,IV) 
        END DO 
      END DO 
C 
C Calculate the coverage for RSD clean and dirty screening option 1. 
C (See page 39 of M6.IM.001) 
C Only peform this for exhaust I/M programs 
C 
      DO IMPGM=1,MAXIMPGM 
        IF(ITEST(IMPGM).GE.20) CYCLE 
        IF (RSDDIRTY(IMPGM).EQ.1.OR.RSDCLEAN(IMPGM).EQ.1) THEN 
          DO IV=1,5 
            DO IAIM=0,MAXYRS-1 
              EPART  = 0.0 
              LAMBDA = 250.0*RSDTEST(IMPGM)*AMAR(IAIM+1,IV) 
     *               / (TNVHCL(IMPGM)*VMTIM(25,IV)) 
              NFCTRL = 1 
              DO ITN = 1,NTF(IMPGM)-1 
                NFCTRL=NFCTRL*ITN 
                IF (NTF(IMPGM).GT.1) 
     *          EPART=EPART+((LAMBDA**ITN)*EXP(-LAMBDA)/NFCTRL) 
              END DO 
C 
C             Subtract out the zero test probability along 
C             with the other 1 thru NTF(IMPGM)-1 test probabilities 
              RSDCVR1(IAIM+1,IV,IMPGM) = 1.- EXP(-LAMBDA)-EPART 
            END DO 
          END DO 
        END IF 
      END DO 
C 
      OLDFRAC=(REAL(MEVAL)+2.0)/12.0 
      NEWFRAC=1.0-OLDFRAC 
C 
      RETURN 
      END 
