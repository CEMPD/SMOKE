      SUBROUTINE BIGIN2(IV,MY,JDX,IDX,IP) 
C 
C  BIGIN2 initializes the vehicle class specific parameters and flags used 
C  subsequently by BIGCFX in the calculation of emission factor correction 
C  factors and offsets. 
C 
C  Called by HRLOOP. 
C 
C  Calls ISPPTR and BIGMET. 
C 
C  Changes: 
C 
C  14 Jun 00 AIR Task 03: Removed BIGIN2 from BIGCFX, since the variables 
C            and arrays are not hourly dependent, as is BIGCFX. Moved 
C            initialization of FCC, FCN, FHC AND FHN to here from BIGIN1. 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C   2 Feb 00 @DynTel-ddj 1-032 Added data value of NEYR to end of 
C            NCATYR to handle vehicle 28 (LDDT34). 
C  10 Mar 00 @DynTel-bg  1-011 added an additional year to NCATYR for second 
C            light diesel category 
C  15 Apr 99 @DynTel-ddj 2-696 Removed initialization of TCF(IB). 
C            Initialization is performed in MAINTCF.  Also modified 
C            initialization of CHFRAC. 
C  22 Oct. 98 @DynTel-HXQ 2-680 delete the last data in NCATYR. 
C  14 Aug 98 @DynTel-zk 2-000  Explicit typing 
C  23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C            certification naming conventions. 
C  20 Apr 98 @DynTel-RJD 2-654 Added an extra dimension to the SALHCF array 
C            to distinguish between starts and running emissions. 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. 
C 
C  Input on call: 
C 
C    constant:  
C    MAXVEH from include file IVTYPE.I 
C    NEYR 
C 
C    argument list: IV,MY,JDX,IDX,IP 
C    common blocks: 
C    /CUMCOM/ CUMMIL 
C    /IVPCOM/ IVPTRL,IVPTRA,IVPTRB 
C    /SPEED4/ LB1STS 
C    /TAMEQ2/ F50K 
C    /VVTYPS/ VTLDDV, VTURBAN   (IVTYPE.I) 
C    /VVNAMS/ MAXVEH (IVTYPE.I) 
C 
C  Output on return: 
C 
C    argument list: 
C    common blocks: 
C    /OPMOD2/ MYSTEP 
C    /OPMOD3/ BFMILE,KINKBF,BFGT50 
C    /SCENE2/ FCC,FCN,FHC,FHN 
C    /SPEED4/ IVB,IVA 
C    /SPEED6/ SALHCF,SCFIDL 
C    /TEMPC3/ CHFRAC 
C 
C  Local array subscripts: 
C 
C  NCATYR(8)       -  NCATYR ( IV ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  NCATYR   I    last non-catalyst MY for the given IVC 
C  LIP      I    Local loop variable for IP 
C 
C  Notes: 
C 
C  SALHCF is reinitialized to 0.0, starting in Version 02. 
C  BIGIN2 was modified for MOBILE4.1. 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'CUMCOM.I' 
      INCLUDE 'IVPCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'OPMOD2.I' 
      INCLUDE 'OPMOD3.I' 
      INCLUDE 'PARAME.I' 
      INCLUDE 'SCENE2.I' 
      INCLUDE 'SPEED4.I' 
      INCLUDE 'SPEED6.I' 
      INCLUDE 'TAMEQ2.I' 
      INCLUDE 'TEMPC3.I' 
C 
C     Declare parameter list 
C 
      INTEGER, INTENT(IN)  :: IV 
      INTEGER, INTENT(IN)  :: MY 
      INTEGER, INTENT(IN)  :: JDX 
      INTEGER, INTENT(IN)  :: IDX 
      INTEGER, INTENT(IN)  :: IP 
C 
C     Declare local variables 
C 
      INTEGER ISR 
      INTEGER LIP 
      INTEGER IVLM 
      INTEGER, DIMENSION(MAXVEH) :: NCATYR=(/1974,1974,1974, 
     *  1978,1978,1983,1983,1983,1983,1983,1983,1983,1983, 
     *  NEYR,NEYR,NEYR,NEYR,NEYR,NEYR,NEYR,NEYR,NEYR,NEYR,NEYR, 
     *  1983,NEYR,NEYR,NEYR/) 
C 
C  Assign operating mode fractions. 
C 
      FCN=0.206 
      FHC=0.273 
      FCC=0.206 
      FHN=(FCC-FCN)+FHC 
C 
C  If needed, calculate the methane offset to composite HC. 
C  The values calculated by BIGMET are needed later in METHANE. 
C 
      IF(IP.EQ.1) CALL BIGMET(IV,MY,IDX,JDX) 
C 
C  Set up the pointers used here and later in BIGSC1 
C   
      IVLM=IVPTRL(IV) 
      IVA=IVPTRA(IV) 
      IVB=IVPTRB(IV) 
C 
      IF(IVLM.EQ.0) RETURN 
C 
      BFMILE=CUMMIL(JDX,IV)*.0001 
C 
C  Select bag coefficients based on vehicle class and model year 
C 
      MYSTEP=1 
      IF(MY.LT.1981.OR.IVLM.GT.5) GOTO 50 
      IF(MY.LT.1984.AND.IVLM.GT.1) GOTO 50 
C 
C  Use yearly stepped bag coefficients for 1981+ LDGV and 1984+ LDGT. 
C  1981+ LDGV class has bag fraction coefficients curve kink at 50K miles. 
C 
C  1981+ LDGV and 1984+ LDGT have kinks at 50+K miles. 
C 
      MYSTEP=2 
      KINKBF=1 
      IF(IVLM.LE.5.AND.BFMILE.GT.F50K) KINKBF=2 
      IF(KINKBF.EQ.2) BFGT50=BFMILE-F50K 
C 
C  Assign cold/hot fractions by model year group. 
C 
  50  IF(MY.GT.NCATYR(IV)) THEN 
         CHFRAC(1)=FCN 
         CHFRAC(3)=FHN 
      ELSE 
         CHFRAC(1)=FCC 
         CHFRAC(3)=FHC 
      END IF 
      CHFRAC(2)=1.-CHFRAC(1)-CHFRAC(3) 
C 
C  Reinitialize SALHCF for the given IDX x IV case to 0.0.  This is done 
C  so that the prior scenario's SALHCF does not show up in the current 
C  scenario's (a possibility with the TF=0.0 bypass of the BIGSAL call). 
C  Also, initialize Idle Speed Correction Factor Adjustments to 1.0. 
C 
      IF(IDLYES) THEN 
        SCFIDL=0.0 
      ELSE 
        DO LIP=1,3 
          DO ISR=1,2 
            SALHCF(IDX,LIP,IV,ISR)=0.0 
          END DO 
        END DO 
      ENDIF 
C 
      RETURN 
      END 
