C ****************************************************************************
C *
C * File Version: @(#)$Id$
C *
C * COPYRIGHT (C) 1996, Daewon Byun
C *
C *
C * Pathname: $Source$
C * Last updated: $Date$ 
C *
C
****************************************************************************/

**********************************************************************
*
*                        SUBROUTINE PBLPWR
*
***********************************************************************
      SUBROUTINE PBLPWR(UL,TL,QL,PL,ZL,ZF,NZMAX, HGT,
     &                  PBL, AMOL, ZRUF, TSFC, PSFC,
     &                  UHGT,THGT,QHGT,PHGT,RHGT)
C
C    SUBROUTINE TO INTERPOLATE U,T,Q,P,RH AT GIVEN HEIGHT FROM
C    LAYERED MM4 MET DATA
C
C    THIS VERSION USES PBL SIMILARITY PROPOSED BY BYUN (1991)
C    TO GET AROUND THE LIMIT OF USING SURFACE LAYER SIMILARITY.
C     --- BY DAEWON BYUN, FEB. 24, 1991
C    
C    DOCUMENTATION MODIFIED  BY SUNNY LO FOR MCIP OUTPUT
C
C**********************************************************************C
C
C    INPUT:
C---- FROM MM4 MET DATA FILE                                  MCIP
C
C           LAYERED WIND SPEED        (M/S) :  UL   = MET_DOT_3D : UWIND (M/S)
C                                                                  VWIND (M/S)
C           LAYERED TEMPERATURE       (K)   :  TL   = MET_CRO_3D : TA   (K)
C           LAYERED MIXING RATIO      (G/G) :  QL   = MET_CRO_3D : QV   (G/G)
C           LAYER PRESSURE            (MB)  :  PL   = MET_CRO_3D : PRES (PA)
C           LAYER HEIGHT              (M)   :  ZL   = GRID_CRO_3D: ZH   (M)
C           LAYER FACE HEIGHT         (M)   :  ZF   = GRID_CRO_3D: ZF   (M)
C
C           NUMBER OF LAYERS          -     :  NZMAX 
C
C           BOUNDARY LAYER HEIGHT     (M)   :  PBL  = MET_CRO_2D : PBL  (M)
C           MONIN-OBUKHOV LENGTH      (M)   :  AMOL = MET_CRO_2D : MOL  (M)
C           ROUGHNESS LENGTH          (M)   :  ZRUF = GRID_CRO_2D: ZZERO(M)
C
C
C           SURFACE AIR TEMPERATURE   (K)   :  TSFC = MET_CRO_2D : TGD  (K)
C           SURFACE PRESSURE          (MB)  :  PSFC = GRID_CRO_2D: P*0  (PA)
C                                                            Plus calculation
C---- HEIGHT AT WHICH INTERPOLATION IS MADE
C           INTERPOLATION HEIGHT      (M)   :  HGT   
C
C
C OUTPUT:
C           INTERPOLATED WIND SPEED   (M/S) :  UHGT
C           INTERPOLATED TEMPERATURE  (K)   :  THGT
C           INTERPOLATED MIXING RATIO (G/G) :  QHGT
C           INTERPOLATED PRESSURE     (MB)  :  PHGT
C           INTERPOLATED RELATIVE HUMIDITY (0-1.0):  RHGT
C
C**********************************************************************C
C
C     PARAMETER (KMAX=15)
      DIMENSION UL(NZMAX),TL(NZMAX),QL(NZMAX),
     &          PL(NZMAX),ZL(NZMAX),ZF(0:NZMAX)
      DIMENSION DZF(NZMAX), THTVL(NZMAX), DTHTVL(NZMAX), DQL(NZMAX)
C
C DATA FOR CONSTANTS - SURFACE LAYER SIMILARITY THEORY
C
C*OLD
      DATA BETAM/4.7/, BETAH/6.3513514/, PRO/0.74/, ALAMDAO /0.07/
C* NEW HOGSTROM (1988)
C      DATA BETAM/6.0/, BETAH/8.2105263/, PRO/0.95/, ALAMDAO /0.3/
C*OLD
      DATA GAMMAM/15.0/, GAMMAH/9.0/, VKAR /0.4/, RICR /0.20/
C* NEW HOGSTROM (1988)
C      DATA GAMMAM/19.3/, GAMMAH/11.6/, VKAR /0.4/, RICR /0.15/
      DATA ETAS /0.1/, DELTAM /0.01/, DELTAH /0.01/, ITMAX /20/
      DATA SMALL /0.0001/, AN1 /0.3/, AN2/0.3/, AN3 /0.8/, CH/0.8/
      DATA AM1 /3.0/, AM2 /0.0/, AM3/0.0/, ZAU /1.0/, ZAT /0.074/
      DATA B4S /13.4/, CU /0.0/, CV /2.0/, ALPHAG /45./
C
      DATA UMIN /0.1/, ZMIN /0.01/, VKAR /0.4/, 
     &     RDRY /287.04/, GAMMA /.0098/, G /9.8/
     &     USTMIN/0.0005/, QMIN/.001/, RBMAX/0.25/
C sl &     USTMIN/0.05/, QMIN/.001/, RBMAX/0.25/
C
C
C***
C FUNCTION DEFINITION
C
C*** FUNCTION VIRTUAL POTENTIAL TEMPERATURE
Corig FVPTMP(T,P,Q)= T*(1.+0.61*Q) * (1000.0/P)**.286
      FVPTMP(T,P,Q)= T*(1000.0/P)**.286
C
C*** FUNCTIONS FOR STABLE AND UNSTABLE CASES
C
      FX(ZETA) = CVMGP(0.0,ABS(1.-GAMMAM*ZETA)**.25,ZETA)
      FY(ZETA) = CVMGP(0.0,ABS(1.-GAMMAH*ZETA)**.5,ZETA)
C
      FPHIM(ZETA)=CVMGP(1.0+BETAM*ZETA,ABS(1.-GAMMAM*ZETA)**(-.25),ZETA)
      FPHIH(ZETA)=CVMGP( PRO*(1.0+BETAH*ZETA),
     &                   ABS(1.-GAMMAH*ZETA)**(-.5), ZETA)
C
      FPSIM(ZETA,X) = CVMGP(-BETAM*ZETA,
     &  2.*ALOG(X+1.)+ALOG(1.+X*X)-2.*ATAN(X), ZETA)
      FPSIH(ZETA,Y) = CVMGP(-BETAH*ZETA, 2.*ALOG(Y+1.), ZETA)
C
      FGAMU(ZETA,X) = CVMGP(-.5*BETAM*ZETA*ZETA,
     & -( 4.*X**3/3.-X**4 +
     &  (X**4-1.)*(ALOG(1.+X*X)+2.*ALOG(1.+X)-2.*ATAN(X)) )/GAMMAM
     & ,ZETA)
C
C FUNCTIONS FOR SURFACE POWER
      FPS(AMU,ETAS,ETAO,XS,XO) = FPHIM(ETAS*AMU)/
     & (ALOG(ETAS/ETAO)-FPSIM(ETAS*AMU,XS)+FPSIM(ETAO*AMU,XO))
      FRS(AMU,ETAS,ETAO,YS,YO) = FPHIH(ETAS*AMU)/
     & (PRO*(ALOG(ETAS/ETAO)-FPSIH(ETAS*AMU,YS)+FPSIH(ETAO*AMU,YO)))
C
C FUNCTIONS FOR COEFFICIENTS FOR PBL PROFILES
C
      FAU(AMU,ETAO,P,XA,XO)=-((P+1.)/(P*(1.-ETAO)**P))*
     &       ( ETAO*ALOG(ETAO)/(1.-ETAO) + 1. - FPSIM(AMU,XA)
     &     + (FGAMU(AMU,XA) - FGAMU(AMU*ETAO,XO))/(AMU*(1.-ETAO)) )
      FAV(AMU,ETAO,ALAMDA,Q)=((Q+1.)/(Q*(1.-ETAO)**(Q+1)))*VKAR/ALAMDA
      FAT(AMU,ETAO,AKAPA,R)= (1./(R*(1-ETAO)**(R-1)))*
     &                       ( VKAR*AKAPA/PRO-FPHIH(AMU)/PRO-1. )
C
C FUNCTIONS FOR THE ESTIMATION OF POWERS
C
      FHM(AMU,ETAS,ETAO,XA,XS,XO) =
     &  ABS( ALOG(ETAS/ETAO) - FPSIM(AMU*ETAS,XS) + FPSIM(AMU*ETAO,XO))
     &/ ABS( ETAO*ALOG(ETAO)/(1.-ETAO) + 1. - FPSIM(AMU,XA)
     &+ (FGAMU(AMU,XA) - FGAMU(AMU*ETAO,XO))/(AMU*(1.-ETAO)) )
      FHH(AMU,ETAS,ETAO,YS,YO,AKAPA) =
     &  ABS( ALOG(ETAS/ETAO) - FPSIH(AMU*ETAS,YS) + FPSIH(AMU*ETAO,YO))
     &/ (ABS( VKAR*AKAPA/PRO-FPHIH(AMU)/PRO-1. )*(1.-ETAO))
C
C FUNCTIONS FOR NORMALIZED PBL PROFILES
C     UPBLN = VKAR*UPBL/USTAR
      UPBLN(AMU,ETA,ETAO,P,AU,X,XO) =
     &                FPSIM(AMU*ETAO,XO)+ ALOG(ETA/ETAO)
     &              - FPSIM(AMU*ETA,X) + AU*(ETA-ETAO)**P
C     VPBLN = VKAR*VPBL/USTAR
      VPBLN(AMU,ETA,ETAO,Q,AV) = -AV*(ETA-ETAO)**Q
      FTAUX(Q,AN) = 1. - AN*(Q+1.-AN**Q)/Q
      FTAUY(AMU,ETA,ETAO,P,AN,X,XO,XA,ALAMDA) =
     & ALAMDA*( ETA*(ALOG(ETA)-1.) - ETAO*(ALOG(ETAO)-1.)   +
     &          (ETAO*ALOG(ETAO)-ETAO+1.)*AN*(P+1.-AN**P)/P +
     &          AN*(AN**P-1.)*(1.-ETAO)*FPSIM(AMU,XA)/P   -
     &          (FGAMU(ETA*AMU,X) - FGAMU(ETAO*AMU,XO))/AMU +
     & AN*(P+1.-AN**P)*(FGAMU(AMU,XA)-FGAMU(AMU*ETAO,XO))/(AMU*P) )/VKAR
C     DTPBLN = VKAR*(TPBL-THTETO)/(PRO*THSTAR)
      DTPBLN(AMU,ETA,ETAO,R,AT,Y,YO) =
     &                 FPSIH(AMU*ETAO,YO) + ALOG(ETA/ETAO)
     &               - FPSIH(AMU*ETA,Y) + AT*(ETA-ETAO)**R
C
C**********************************************************************C
C
      ROG=RDRY/G
C
C COMPUTE SURFACE VIRTUAL POTENTIAL TEMPERATURE, SPECIFIC HUMIDITY
C
C*      ES     = 6.1078*EXP( 5384.21*(1./273.15 - 1./TL(1)) )
C* REPLACED MIN(TL(1),TSFC) TO ESTIMATE SURFACE VAPOR PRESSURE
C*DWB 7/24/89

      ESAT = 6.1078*EXP( 5384.21*(1./273.15 - 1./TL(1)) )
      QSAT = 0.622*ESAT/(PL(1)-ESAT)
      RHSFC = QL(1)/QSAT
Cdon      TES = AMIN1(TSFC,TL(1))
Cdon      ES     = 6.1078*EXP( 5384.21*(1./273.15 - 1./TES) )
      ES     = 6.1078*EXP( 5384.21*(1./273.15 - 1./TSFC) )
C*
Cdon      QSFC   = 0.622*ES/(PSFC-ES)
      QSFC   = RHSFC * 0.622*ES/(PSFC-ES)
      VPTSFC = FVPTMP(TSFC,PSFC,QSFC)
C
C COMPUTE VIRTUAL POTENTIAL TEMPERATURE FOR WHOLE LAYERS
C
      DO 100 K=1,NZMAX
      THTVL(K) = FVPTMP(TL(K),PL(K),QL(K))
      DTHTVL(K)= THTVL(K) - VPTSFC
      DQL(K)   = QL(K) - QSFC
  100 CONTINUE

C     write(6,*) ' ptemp computed after do 100'
C
C COMPUTE BULK RICHARDSON NUMBER RB
C
C!!!! RB IS NOT USED FOR THE CURRENT VERSION SINCE PBL/AMOL !!!!
C!!!! IS AVAILABLE (2/23/91)                                !!!!
C
C  RB   = G * Z * DTH / ( THO * U**2 )
C      RB = G*ZL(1)*(THTVL(1) - VPTSFC) / (VPTSFC*UL(1)*UL(1))
C      RB = AMIN1(RB,RBMAX)
C
C DETERMINE LAYER INDEX FOR PBL
C
      DO 200 K=1,NZMAX
      IF(PBL.LT.ZL(K)) THEN
       KPBL=K
       GOTO 210
      END IF
  200 CONTINUE
      KPBL=NZMAX
  210 CONTINUE
C
C DETERMINE LAYER INDEX FOR HGT
C
      DO 250 K=1,NZMAX
      IF(HGT.LT.ZL(K)) THEN
       KHGT=K
       GOTO 260
      END IF
  250 CONTINUE
      KHGT=NZMAX
  260 CONTINUE
C
C COMPUTE DZF
C
      DO 300 K=1,NZMAX
      DZF(K) = ZF(K) - ZF(K-1)
  300 CONTINUE
C
C FIND POWERS P AND R
C
C*** FOR VARIABLE VALUE OF P ( P,R NOT DEPENDENT ON AMU)
      AKAPA= 0.0 ! WE ASSUME STRATIFICATION ALOFT IS NOT IMPORTANT
      ETAS = 0.1
      ETAO = ZRUF/PBL
      AMU  = PBL/AMOL
      XO   = FX(ETAO*AMU)
      XS   = FX(ETAS*AMU)
      XA   = FX(1.0 *AMU)
      YO   = FY(ETAO*AMU)
      YS   = FY(ETAS*AMU)
      YA   = FY(1.0 *AMU)
      ANS = (ETAS-ETAO)/(1.-ETAO)
c     print *, 'Passing Check Point 10'
      ALOGNS = ALOG(ANS)
* start to change - keep P1, R1 C
*     write(6,*) 'compute p and r'
*kay      IF(ABS(AMOL).GT.1000.) THEN
C     write(6,*) 'neutral case'
C ! TREAT AS NEUTRAL
                     P1 = 1.43
                     R1 = P1
*kay  ELSE
C ! DIABATIC CONDITION
C
*      PS   = FPS(AMU,ETAS,ETAO,XS,XO)
*      RS   = FRS(AMU,ETAS,ETAO,YS,YO)
C
C** ESTIMATE P AND R USING THE ASYMPTOTIC MAGNITUDE COMPARISONS
C
C** USE ITERATION FOR P STARTING WITH INITIAL P = PS
C
*     P1 = PS
*     FHM1 = FHM(AMU,ETAS,ETAO,XA,XS,XO)
*     DO 520 IT=1,ITMAX
*     WRITE(*,*) 'P1 = ',P1
C*      PNEW=ALOG(ABS(P1*DELTAM*FHM1/(1.+P1)))/ALOGNS
C* NEWTON-RAPSON METHOD
*     GP = ALOG((P1+1.)/P1) + P1*ALOGNS - ALOG(DELTAM*FHM1)
*     G1P= ALOGNS - 1./(P1*(1.+P1))
*     PNEW = P1 - GP/G1P
*     IF(ABS((PNEW-P1)/PNEW).LT.SMALL*0.01) THEN
*       P1  = PNEW
*ITP = IT
*       GOTO 521
*     ELSE
*       P1 = PNEW
*     END IF
* 520 CONTINUE
C
C NOT CONVERGED WITHIN ITMAX
*     WRITE(10,*) '*************************************************'
*     WRITE(10,*) ' CONVERGENCE FAILED FOR POWER P '
*     WRITE(10,*) '*************************************************'
*     STOP
* 521 CONTINUE
C     WRITE(6,*) 'P1 IS DETERMINED = ',P1
C
C
C** USE ITERATION FOR R STARTING WITH INITIAL R = RS
C
*     R1 = RS
*     FHH1 = FHH(AMU,ETAS,ETAO,YS,YO,AKAPA)
*     DO 540 IT=1,ITMAX
C     WRITE(*,*) 'R1 = ',R1
C*      RNEW = ALOG(R1*DELTAH*FHH1)/ALOGNS
C* USE NEWTON-RAPSON METHOD
*     GR   = -ALOG(R1) + R1*ALOGNS - ALOG(DELTAH*FHH1)
*     G1R  = -1./R1 + ALOGNS
*     RNEW = R1 - GR/G1R
*     IF(ABS((RNEW-R1)/RNEW).LT.SMALL*0.01) THEN
*       R1  = RNEW
*       ITR = IT
*       GOTO 541
*     ELSE
*       R1 = RNEW
*     END IF
* 540 CONTINUE
C
C NOT CONVERGED WITHIN ITMAX
*     WRITE(10,*) '*************************************************'
*     WRITE(10,*) ' CONVERGENCE FAILED FOR POWER R '
*     WRITE(10,*) '*************************************************'
*     STOP
* 541 CONTINUE
C     WRITE(6,*) ' R1 IS DETERMINED  ',R1
*     END IF
C end here
      P = P1
      R = R1
C DIABATIC P,R ALSO DETERMINED
      AU   = FAU(AMU,ETAO,P,XA,XO)
      AT   = FAT(AMU,ETAO,AKAPA,R)
C
C COMPUTE USTHAT
C
      IF(KPBL.LE.1) THEN
         K   = 1
         ETA = ZL(K)/PBL
c     print *, 'Passing Check Point 11'
         X   = FX(ETA*AMU)
         Y   = FY(ETA*AMU)
	 USTAR  = VKAR*UL(K)/UPBLN(AMU,ETA,ETAO,P,AU,X,XO)
         THSTAR = VKAR*DTHTVL(K)/DTPBLN(AMU,ETA,ETAO,R,AT,Y,YO)
         QSTAR  = VKAR*DQL(K)/DTPBLN(AMU,ETA,ETAO,R,AT,Y,YO)
C
         USTHAT = AMAX1(USTAR,USTMIN)
         THSTHAT= THSTAR
         QSTHAT = QSTAR
c     print *, 'Passing Check Point 12'
      ELSE
        USDZSUM = 0.0
        THDZSUM = 0.0
        QDZSUM  = 0.0
        DZSUM   = 0.0

        DO 400 K=1,KPBL-1
         ETA = ZL(K)/PBL
         X   = FX(ETA*AMU)
         Y   = FY(ETA*AMU)
c     print *, 'Passing Check Point 13'
C
C UST, THST, QST
	 USTAR  = VKAR*UL(K)/UPBLN(AMU,ETA,ETAO,P,AU,X,XO)
         THSTAR = VKAR*DTHTVL(K)/DTPBLN(AMU,ETA,ETAO,R,AT,Y,YO)
         QSTAR  = VKAR*DQL(K)/DTPBLN(AMU,ETA,ETAO,R,AT,Y,YO)
c     print *, 'Passing Check Point 14'
C
         USDZSUM = USDZSUM + DZF(K)*USTAR
         THDZSUM = THDZSUM + DZF(K)*THSTAR
         QDZSUM  = QDZSUM  + DZF(K)*QSTAR
         DZSUM   = DZSUM   + DZF(K)
  400 CONTINUE
C
        USTHAT  = AMAX1((USDZSUM/DZSUM),USTMIN)
        THSTHAT = THDZSUM/DZSUM
        QSTHAT  = QDZSUM/DZSUM
c     print *, 'Passing Check Point 15'
C
      END IF
C
C      write(6,*) ' usthat,thsthat,qsthat =',usthat,thsthat,qsthat
C
C NOW, WE HAVE GOOD ESTIMATION OF USTHAT, THSTHAT, QSTHAT
C DEPENDING ON THE RELATIVE POSITION OF PBL AND HGT, ESTIMATE
C METEOROLOGICAL PARAMETERS USING THE PROFILE FUNCTION AND LINEAR
C INTERPOLATION
C
      IF(HGT.LE.PBL) THEN
C
C HGT IS IN PBL LAYER: USE PARAMETERIZED PBL PROFILE
C
      ETA = HGT/PBL
      X   = FX(ETA*AMU)
      Y   = FY(ETA*AMU)
      UHGT = USTHAT*UPBLN(AMU,ETA,ETAO,P,AU,X,XO)/VKAR
      THVHGT=VPTSFC+ THSTHAT*DTPBLN(AMU,ETA,ETAO,R,AT,Y,YO)/VKAR
      QHGT  =QSFC  + QSTHAT *DTPBLN(AMU,ETA,ETAO,R,AT,Y,YO)/VKAR
c     print *, 'Passing Check Point 16'
C
      UHGT  = AMAX1( UHGT, UMIN )
      QHGT  = AMAX1( QHGT, QMIN )
C
      ELSE IF( (PBL.LT.HGT) .AND. (HGT.LT.ZL(KPBL)) ) THEN
C
C HGT IS ABOVE THE BOUNDARY LAYER BUT BELOW ZL(KPBL):
C   USE PARAMETERIZED PROFILE TILL ZSL AND INTERPOLATE WITH
C   VALUES AT PBL AND AT ZL(KPBL)
C
C COMPUTE ETA, X, AND Y AT PBL
      ETA = 1.0
      X   = FX(ETA*AMU)
      Y   = FY(ETA*AMU)
      UPBL = USTHAT*UPBLN(AMU,ETA,ETAO,P,AU,X,XO)/VKAR
      THVPBL=VPTSFC+ THSTHAT*DTPBLN(AMU,ETA,ETAO,R,AT,Y,YO)/VKAR
      QPBL  =QSFC  + QSTHAT *DTPBLN(AMU,ETA,ETAO,R,AT,Y,YO)/VKAR
c     print *, 'Passing Check Point 17'
C
C INTERPOLATE
C       write(6,*) 'interpolating from zbar to hgt'
        CORA   = (ZL(KPBL)-HGT) / (ZL(KPBL) -PBL)
        CORB   = 1. - CORA
C
        UHGT   = AMAX1( (UL(KPBL)*CORB + UPBL*CORA), UMIN )
        THVHGT = THTVL(KPBL)*CORB + THVPBL*CORA
        QHGT   = AMAX1( (QL(KPBL)*CORB + QPBL*CORA), QMIN )
c     print *, 'Passing Check Point 18'
C
      ELSE
C
C HGT IS LOCATED ABOVE THE MIXED LAYER:
C   USE SIMPLE LINEAR INTERPOLATION USING ADJACENT LEVELS
C
C       write(6,*) 'interpolating using zl(khgt-1) and (khgt) for hgt'
        CORA   = (ZL(KHGT)-HGT) / (ZL(KHGT) -ZL(KHGT-1))
        CORB   = 1. - CORA
C
        UHGT   = AMAX1( (UL(KHGT)*CORB + UL(KHGT-1)*CORA), UMIN )
        THVHGT = THTVL(KHGT)*CORB + THTVL(KHGT-1)*CORA
        QHGT   = AMAX1( (QL(KHGT)*CORB+QL(KHGT-1)*CORA), QMIN )
c     print *, 'Passing Check Point 19'
      END IF
C
C**********************************************************************C
C
C CONVERT VIRTUAL POTENTIAL TEMPERATURE TO NORMAL TEMPERATURE
C
      PHGT  = PL(KHGT)*EXP( (ZL(KHGT)-HGT)/(ROG*THTVL(KHGT)) )
Cold  THGT  = (THVHGT/(1.+0.61*QHGT))*(PHGT/1000.)**.286
      THGT  = (THVHGT)*(PHGT/1000.)**.286
C
C COMPUTE RELATIVE HUMIDITY
C
      ESHGT = 6.1078*EXP(5384.21*(1./273.15-1./THGT))
      QSHGT = 0.622*ESHGT/(PHGT-ESHGT)
      EHGT  = PHGT*QHGT/(QHGT+0.622)
      RHGT  = EHGT/ESHGT
c     print *, 'Passing Check Point 20'
C
C EVERYTHING HAS BEEN COMPUTED
C
      RETURN
      END
