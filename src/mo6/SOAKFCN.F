      SUBROUTINE SOAKFCN(IP) 
C 
C  SOAKFCN creates an array to store the California soak function 
C  and then corrects it for EPA use. This array is used in the SOAKCF  
C  function where it is weighted by the fractions of vehicles using 
C  different catalyst types. 
C 
C  Called by HRLOOP 
C 
C  Calls 
C 
C  Changes: (Last change first) 
C 
C  28 may 99 @DynTel -bg 2-650 Added code to write out results for a check with 
C            EPA documents. Output code is commented out. 
C   3 May 99 @DynTel-bag  2-650 MINUTES now an array, SOK_TIME, rather than 
C            determining each value for each ST_INDX and IP; rearanged 
C            assignment of soak time interval parameters and EPA correction 
C            factor. 
C  20 Oct 98 @DynTel-bag 2-000 Modified (CALSOAK and renamed to SOAKFCN) to 
C            eliminate redundant calculations over hour, vehicle type, ... 
C            loops in HRLOOP and reduce the MOBILE6 run time. Moved TMINS, 
C            calculation parameters, CSA, CSB, CSC, and RATIO into 
C            block data. 
C 
C  Input on call or from calls: 
C 
C     parameter list: IP 
C     common blocks: 
C     /SOAK/ CSA, CSB, CSC, RATIO, TMINS, SOK_TIME 
C 
C     Output on return: 
C     common blocks: 
C     /SOAK/ SOAKFUN 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------ 
C  IP        I     pollutant index  (1=HC,2=CO,3=NOx) 
C  ST_INDX   I     soak times index (1,69) 
C  IECT      I     emission control technology index 
C                                 (1=no cat., 2=cat., 3= elec. cat.)  
C  CN        I     denotes curves 1 and 2 in CARB soak function 
C  EPACF     R     EPA correction factor for the CARB soak function 
C  FLT_MINS  R     local, real value for TMINS(IP,IECT) used to determine 
C                  soak time range 
C 
C  Notes. 
C   
C  18 Mar 98 @DynTel-RJD 6-650 Created this subroutine for STARTS.      
C   4 Aug 98 @DynTel-MLA 6-650 Merged this subroutine into Mobile. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'SOAK.I' 
C 
C      REAL EPACF(69) 
C      REAL CAFCN(69) 
C      INTEGER CN(69) 
C 
      INTEGER CN 
      INTEGER IECT 
      INTEGER IP 
      INTEGER ST_INDX 
      REAL    EPACF 
      REAL    FLT_LMIT 
C 
C 
CCCC        WRITE(6,1113) IP,ratio(ip) 
CCCC1113    FORMAT('SOAKFCN  *****IP= ',i2,2x, 'r=',1p,E10.3 ) 
C 
      DO IECT = 1,2 
C 
C       There Two soak time ranges for the California soak function. 
C       The California soak function has a steeper rise for soak times 
C       less than 120 minutes. The soak time which determines the 
C       beginning of the more gradual increase in the California soak 
C       function with increasing soak time is pollutant and catalyst 
C       dependent. So an array of soak time upper range lower limits, 
C       UPR_RANG_LO_LMT, which are (pollutant and catalyst dependent 
C       and thus) related to the California soak function parameters, CSA, 
C       CSB, and CSC is defined in BD41. (see M6.ste.003 table 4 and BD41; 
C       also see m6.ste.003 June 30, 1998 revision) 
C 
C       from BD41: 
C       DATA UPR_RANG_LO_LMT/ 
C 
C       IECT    IP = 1(HC)   2(CO)   3(NOX) 
C 
C        1            53,     120,     120, 
C        2            90,     117,      62 / 
C 
C 
        FLT_LMIT = FLOAT(UPR_RANG_LO_LMT(IP,IECT)) 
C 
        DO ST_INDX = 1,69 
C         Determine the EPA correction factor to California soak function 
C         for only the early soak time ranges. The EPA correction factor 
C         for the latter/second soak time ranges is 1.0, i.e., no correction. 
C         Default to the second range of soak time intervals, i.e., no 
C         correction to the CA soak function : 
C         IF( SOK_TIME(ST_INDX) .GE. UPR_RANG_LO_LMT ) THEN 
C 
C         Default adjustment to California soak function is no adjustment, 
C         i.e., soak time is in the upper soak time range 
C 
CCCC          CN(ST_INDX) = 2 
CCCC          EPACF(ST_INDX) = 1.0 
C 
          CN = 2 
          EPACF = 1.0 
C 
          IF( SOK_TIME(ST_INDX) .LT. FLT_LMIT ) THEN 
C 
CCCC            CN(ST_INDX) = 1 
C 
            CN = 1 
C 
C           Determine adjustment to California soak function 
C 
            IF( SOK_TIME(ST_INDX) .LE. 10. ) THEN 
C 
CCCC              EPACF(ST_INDX) = RATIO(IP)+(1.0-RATIO(IP))* 
C 
              EPACF = RATIO(IP)+(1.0-RATIO(IP))* 
     &                ((10.-SOK_TIME(ST_INDX) )/10.) 
C 
            ELSE IF( SOK_TIME(ST_INDX) .GT. 10.) THEN 
C 
CCCC              EPACF(ST_INDX) = RATIO(IP)+(1.0-RATIO(IP))* 
C 
              EPACF = RATIO(IP)+(1.0-RATIO(IP))* 
     &             ((SOK_TIME(ST_INDX)-10.)/(FLT_LMIT-10.)) 
C 
            END IF 
C 
          END IF 
C 
C         Determine the reduced California soak function from the above EPA 
C         correction factor, EPACF, and the CA soak function parameters. 
C 
CCCC          cafcn(st_indx) = 
CCCC     &             CSA(CN(ST_INDX),IP,IECT) + 
CCCC     &             CSB(CN(ST_INDX),IP,IECT)*SOK_TIME(ST_INDX) + 
CCCC     &             CSC(CN(ST_INDX),IP,IECT)*SOK_TIME(ST_INDX)* 
CCCC     &             SOK_TIME(ST_INDX) 
CCCC          SOAKFUN(ST_INDX,IECT) = EPACF(ST_INDX) * cafcn(st_indx) 
C 
          SOAKFUN(ST_INDX,IECT) = EPACF * 
     &           ( CSA(CN,IP,IECT) + 
     &             CSB(CN,IP,IECT)*SOK_TIME(ST_INDX) + 
     &             CSC(CN,IP,IECT)*SOK_TIME(ST_INDX)* 
     &             SOK_TIME(ST_INDX) ) 
C 
        END DO 
C       end loop over soak time index 
C 
CCCC        WRITE(6,*) '**********IECT= ',IECT 
CCCC        WRITE(6,1887) ' ST_INDX :',(ST_INDX,ST_INDX =1,69) 
CCCC        WRITE(6,1870) 'SOK_TIME :',(SOK_TIME(ST_INDX),ST_INDX =1,69) 
CCCC        WRITE(6,1887) 'CURVE # :',(CN(ST_INDX),ST_INDX =1,69) 
CCCC        WRITE(6,1871) 'CSA :',(CSA(CN(ST_INDX),IP,IECT),ST_INDX=1,69) 
CCCC        WRITE(6,1871) 'CSB :',(CSB(CN(ST_INDX),IP,IECT),ST_INDX=1,69) 
CCCC        WRITE(6,1871) 'CSC :',(CSC(CN(ST_INDX),IP,IECT),ST_INDX=1,69) 
CCCC        WRITE(6,1871) '  CAFCN :',(cafcn(ST_INDX),ST_INDX=1,69) 
CCCC        WRITE(6,1871) '  EPACF :',(EPACF(ST_INDX),ST_INDX =1,69) 
CCCC        WRITE(6,1871) 'SOAKFCN :',(SOAKFUN(ST_INDX,iect),ST_INDX=1,69) 
CCCC        WRITE(6,*) ' ' 
C 
      END DO 
C     end loop over catalyst type 
C 
CCCC1870    FORMAT(1X,A12,69(F8.3,3X)) 
CCCC1871    FORMAT(1X,A12,69(1P,E10.3,1X)) 
CCCC1887    FORMAT(2X,A11,3X,69(I3,8X)) 
C 
      RETURN 
      END        
