C 
C  Changes: (Last change first) 
C 
C 15 Nov 99 @DynTel-BG 1-006 Renamed EFRSTL to G_MI_EVP 
C 10 Mar 98 DynTel@RJD 2-622 Added integer labels for emission factors to 
C           facilitate the efficient calculation of emission factors for 
C           composite vehicles in EFCVT. 
C 14 Apr 97 - DynTel@MLA  Expanded arrays EFIDLE, G_MI_EVP, RLGGAL, and  
C           RSTGPH to accomodate the new BUS vehicle types. 
C  3 Dec 96 DynTel@MLA Expanded arrays EFIDLE, G_MI_EVP, RLGGAL, and 
C           RSTGPH to accomodate the new Heavy Duty Gasoline, Heavy 
C           Duty Diesel, and LDGT vehicle categories. 
C 
C 
C 
C The data is initialized in BD27. 
C 
      INTEGER IRLGGAL, IRSTGPH, IGMIEVP, IEFIDLE 
      REAL RLGGAL, RSTGPH, G_MI_EVP, EFIDLE 
      COMMON /RESUL2/ RLGGAL(29),RSTGPH(29),G_MI_EVP(6,29), 
     *               EFIDLE(3,29),IRLGGAL,IRSTGPH,IGMIEVP,IEFIDLE 
