      RECURSIVE REAL FUNCTION GLL_WGT(MY,ACTUAL_MY,HOT_DIURN_RUN,AGE) 
C 
C GLL_WGT calculates the HOTSOAK Gross Leakers weighting fraction 
C       It is written so that it can be used for HOTSOAK Diurnal and 
C       HOTSOAK Running Loss 
C 
C CALLED by HOTSOK 
C 
C Calls GLL_WGT recursively. 
C 
C Changes Last First 
C 
C  13 Dec 00 @EPA-bag t-x6 MY and ACTUAL_MY added to parameter list for CAA switch 
C 18 Nov 1999 @Dyntel-nh 1-014 New function 
C 
C Input Arguments 
C 
C MY INTEGER Model Year 
C HOT_DIURN_RUN INTEGER. Steering index 1=HOTSOAK, 2=DIURNAL, 3=RUNNINGLOSS 
C AGE   INTEGER Veicle age 
C 
C OUTPUT ARGUMENTS 
C 
C FUNCTION VALUE Weighting 
C 
      USE LEAKYDAT, ONLY : ILIQUID_LEAKS,GLL_COEF,PHASE_IN 
      IMPLICIT NONE 
C 
C Input Arguments 
C 
      INTEGER, INTENT(IN) :: MY        ! Model Year 
      INTEGER, INTENT(IN) :: ACTUAL_MY 
      INTEGER HOT_DIURN_RUN !=2 DIURNAL =3 RUN 
      INTEGER, INTENT(IN) :: AGE 
      INTEGER J  ! Phase IN INDEX 
      INTEGER K! DIURNAL=1, RUNNINGLOSS=3, HOTSOAK =NOT APPLICABLE 
      REAL WGTDIURN 
      REAL WGTRUN 
      REAL HOTS_GLL 
      REAL A,B,C,D ! Coefficients for calculating weighting 
C 
      REAL APRE 
      REAL BPRE 
      REAL CPRE 
      REAL DPRE 
C 
      REAL GLLPRE 
      REAL APOST 
      REAL BPOST 
      REAL CPOST 
      REAL DPOST 
C 
      REAL GLLPOST 
      INTEGER IDIURNAL 
      INTEGER IRUNLOSS 
      INTEGER IHOTSOAK 
      PARAMETER (IHOTSOAK=1,IDIURNAL=2,IRUNLOSS=3) 
      INTEGER, PARAMETER :: INO=0 
C 
      GLL_WGT=0.0 
      IF(ILIQUID_LEAKS == INO ) GOTO 99 
C 
      SELECT CASE (HOT_DIURN_RUN) 
          CASE (IHOTSOAK) 
            WGTDIURN = GLL_WGT(MY,ACTUAL_MY,IDIURNAL,AGE) 
            WGTRUN   = GLL_WGT(MY,ACTUAL_MY,IRUNLOSS,AGE) 
            HOTS_GLL = (WGTDIURN+WGTRUN)-(WGTDIURN*WGTRUN) 
            GLL_WGT=HOTS_GLL 
            RETURN 
          CASE (IDIURNAL,IRUNLOSS) 
 
            K=2*(HOT_DIURN_RUN-2)+1 
 
            SELECT CASE (MY) 
                CASE (:1995) ! Use One coefficient set 1995 and earlier 
                    A=GLL_COEF(1,K) 
                    B=GLL_COEF(2,K) 
                    C=GLL_COEF(3,K) 
                    D=GLL_COEF(4,K) 
C 
                    GLL_WGT=A/(1+B*EXP((C*AGE)+(D*AGE*AGE))) 
C 
                CASE (1999:) ! Another Coefficient set 1999 and later. 
                    A=GLL_COEF(1,K+1) 
                    B=GLL_COEF(2,K+1) 
                    C=GLL_COEF(3,K+1) !This C is half that of MY95-, for durability 
                    D=GLL_COEF(4,K+1) 
C 
                    GLL_WGT=A/(1+B*EXP((C*AGE)+(D*AGE*AGE))) 
C 
                CASE (1996:1998)! Transition period uses both 
C 
                    APRE=GLL_COEF(1,K) 
                    BPRE=GLL_COEF(2,K) 
                    CPRE=GLL_COEF(3,K) 
                    DPRE=GLL_COEF(4,K) 
C 
                    GLLPRE=APRE/(1+BPRE*EXP((CPRE*AGE) 
     &              +(DPRE*AGE*AGE))) 
C 
C 
                    APOST=GLL_COEF(1,K+1) 
                    BPOST=GLL_COEF(2,K+1) 
                    CPOST=GLL_COEF(3,K+1)  !This C is half that of MY95-, for durability 
                    DPOST=GLL_COEF(4,K+1) 
C 
                    GLLPOST=APOST/(1+BPOST*EXP((CPOST*AGE) 
     &              +(DPOST*AGE*AGE))) 
C 
                    J=ACTUAL_MY-1995 
                    GLL_WGT=PHASE_IN(J)*GLLPRE + 
     &                           (1-PHASE_IN(J))*GLLPOST 
 
             END SELECT 
C                 
      END SELECT 
      IF (GLL_WGT<0.0) THEN 
        GLL_WGT=0.0 
      ELSEIF (GLL_WGT>1.0) THEN 
        GLL_WGT=1.0 
      END IF 
C 
99    RETURN 
      END 
