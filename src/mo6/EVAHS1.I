C 
C Changes : (Last change first) 
C 
C  27 Feb 01 @EPA-BAg bug314 removed HSIV5 and RLV5 variables. 
C  20 Jul 00 @EPA-elg  Added HSIV5 and RLV5 variables. 
C  10 Mar 00 @DynTel-mla 1-013  Added hot soak fraction array, HSRED 
C 10 Nov 99 @Dyntel nh 1-014 added HSEQ_81 AND HSEQ_86 declarations 
C               AND TO common /EVAHS1/ 
C  15 Jan 99 @DynTel-bag 2-664  Added hot soak fraction array, HSFRAC(IMIN), 
C            rather than calculating them in HOTSOK 
C 24 Sep 98 - DynTel@MLA 2-670  Added HSCOEF, the coefficients for the 
C             hourly hot soak distribution function. 
C  3 Apr 97 - DynTel@MLA  Expanded array EFHS to accomodate the new GASBUS  
C             vehicle type. 
C 15 Nov 96 DynTel@MLA Expanded array EFHS to accomodate the new Heavy Duty 
C           Gasoline and LDGT vehicle categories 
C 
C 
C  The data is initialized in BD11. 
C 
      REAL HSEQ, EFHS, HSCOEF, HSFRACT, HS_HRLY_TRP_FRCN, 
     *     HS_TRPS_PR_DAY,HSEQ_81,HSEQ_86 
C 
      REAL,DIMENSION(2),PARAMETER  ::  HSRED = (/0.739,0.714/) 
C 
      COMMON /EVAHS1/ HSEQ(3,2,9), EFHS(3,2,5,15), HSCOEF(4,14,2), 
     *                HSFRACT(60), HS_HRLY_TRP_FRCN(24,2), 
     *                HS_TRPS_PR_DAY(2,2),HSEQ_81(3,2,9),HSEQ_86(3,2,9) 
