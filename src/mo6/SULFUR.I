C 
C  Changes: (Last change first) 
C 
C  10 Jul 00 @EPA-ddj Bug184 Added GLFUEL_P. 
C  22 Jun 00 @EPS-ddj bug157 Increased USSCF array IV dimension 
C                     USSCF(4,2,2,4,2)-->USSCF(4,2,2,4,3) 
C    4 May 00 @EPA-ddj, E-7  Added GLS_PPM. 
C   21 May 99 @DynTel-ddj, 2-699  NEW Include FILE. 
C 
C 
C  VARIABLE Dictionary 
C 
C  ARRAYS 
C    S_COEF(IV, IP, ISR, ISTD, INH)  Regression coefficients for sulfur 
C                                    correction factor. 
C    USSCF(ISTD, ISR, INH, IP, IV)   User supplied sulfur correction 
C                                    factor. 
C    S_FUEL(LEVEL,FUEL_P,ICY/MY,IT2) Fuel sulfur level matrix. 
C    LONGTERM(IP)                    Ratio of long-term to short-term 
C                                    sulfur sensitivity for LEV and 
C                                    cleaner vehicles. 
C 
C  CONSTANTS 
C    IH_S_PPM                        Gasoline sulfur content (ppm), 
C                                    in-house use only. 
C    S_PPM                           Gasoline sulfur content (ppm). 
C    SUL_FLAG                        Flag indicating source of sulfur 
C                                    correction factors: 
C                                      1,  use database regression 
C                                          coefficients. 
C                                      2,  user supplied sulfur 
C                                          correction factors (USSCF). 
C    FUEL_P                          Index indicating Fuel Program: 
C                                      1,  Eastern Conventional Gasoline 
C                                      2,  RFG 
C                                      3,  Western Conventional Gasoline 
C                                      4,  Alternate user supplied 
C 
C 
C     The data is initialized in BD05 
C 
      INTEGER GLFUEL_P 
      INTEGER SUL_FLAG 
C 
      REAL IH_S_PPM 
      REAL GLS_PPM 
      REAL S_COEF 
      REAL S_PPM 
      REAL USSCF 
      REAL S_FUEL 
      REAL LONGTERM 
C 
      COMMON /S_DATA/ IH_S_PPM,     S_PPM,         USSCF(4,2,2,4,3), 
     &                SUL_FLAG,     S_COEF(4,2,4,2,2), 
     &                LONGTERM(3),  S_FUEL(2,4,16,2), GLS_PPM, 
     &                GLFUEL_P 
