C 
C  Changes : (Last change first) 
C 
C 19 Feb 01 AIR Task 11: Added IEFPM, particules. 
C 29 Aug 00 @EPA-elg Removed PCTREM 
C 27 Jul 00 @EPA-djb Added PCTREM and removed HOUREF. 
C 15 Nov 99 @DynTel-BG 1-006 Moved RNLOSSTEM out of common to RNLOSS parameter list 
C 12 Nov 99 @DynTel-BG  1-006  Removed hourly hot soak array, HSHOUR, daily 
C           diurnal results, DIURNL24, hourly resting loss results, and 
C           hourly diurnal results 
C  5 May 1999 @DynTel-HQ 2-664 Add RNLOSSTEM variable 
C  3 Mar 99 @DynTel-bag 2-664  Increased hourly arrays to 24 hourly intervals; 
C           Renamed/deleted diurnal, resting, and running loss variables 
C 30 Oct 98 @DynTel-bag 2-677  Added declaration of DIURNL24(LEAKCAT,IDX), 
C           the full day diurnals (FDD) dependent upon three vapor leak 
C           categories (pressure and purge tests pass and fail rates) 
C           and liquid leaks, DIURNL_VL(IH,IV) and DIURNL_LL(IH,IV) hourly diurnal 
C           totals for vapor and liquid leaks, respectively and the daily 
C           total, DIURNL24TOT(IV). Added hourly resting loss variables, 
C           RSTLSS(IH,IV)and RSTLSSLL(IH,IV). 
C 24 Sep 98 @DynTel-MLA 2-670  Added declaration of HSHOUR(IH,IVGAS), the 
C           hourly hot soak emission factors. 
C  5 Aug 98 @DynTel-MLA 2-650  Added declaration of HOUREF(IP,IH,IV,ISR), 
C           the hourly start and running exhaust emnission factors. 
C 10 Mar 98 DynTel@RJD 2-622 Added integer labels for emission factors to 
C           facilitate the efficient calculation of emission factors for 
C           composite vehicles in EFCVT. 
C 14 Apr 97 DynTel@MLA  Expanded arrays EFEVAP, EFEXH, EFFTP, REFL_GM_MI,  
C           and EFRUNL to accomodate the new BUS vehicle types. 
C  3 Dec 96 DynTel@MLA Expanded arrays EFEVAP, EFEXH, EFFTP, REFL_GM_MI, and  
C           EFRUNL to accomodate the new Heavy Duty Gasoline, Heavy Duty 
C           Diesel, and LDGT vehicle categories. 
C 
C The data is initialized in BD27. 
C 
      INTEGER IEFFTP, IEFEXH, IEFEVAP, IEFREFL, IEFRUNL, IEFPM 
      REAL  EFFTP, EFEXH, EFEVAP, REFL_GM_MI, EFRUNL 
      COMMON /RESUL1/ EFFTP(3,29),EFEXH(29),EFEVAP(29), 
     *                REFL_GM_MI(29),EFRUNL(29), IEFFTP, 
     *                IEFEXH, IEFEVAP, IEFREFL,IEFRUNL, IEFPM 
