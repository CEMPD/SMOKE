C Data is initialized in BD_IM, REINIT and INIFLG 
C 
C 15 Sep 00 AIR Task 08 Increased I/M-related array sizes from 5 to 7 to 
C           allow for EVAP I/M 
C 16 Aug 99 @DynTel-bg  2-663 Renamed SIMFLG and SIMONLY to IMFLG_NP and 
C           IMONLY_NP, respectively. These are for the 'non-periodic' COIM 
C           and OBD programs. 
C 28 Oct 98 @DynTel-RJD 2-663 LDGSIM contains data used the 1981-1993 light 
C           duty gasoline vehicle and truck selective IM benefits calculations 
C 
      REAL COIM,RSDEFF,RSDCUT,RSDCVR1,RSDCVR2,RSDNOX,HCLOSS, 
     *     COLOSS,NOXLOSS,RSDFAIL,NVHCL,RSDTEST,TNVHCL,RSDCSE, 
     *     IMONLY_NP 
      INTEGER IMFLG_NP,NTF,IMPN_NP 
      COMMON /LDGSIM/ IMFLG_NP(7),NTF(7),RSDEFF(15,3), 
     *RSDCUT(7),RSDCVR1(25,5,7),RSDCVR2(25,7), 
     *RSDFAIL(25,7),NVHCL(25,7),COIM(7),RSDTEST(7),TNVHCL(7), 
     *RSDNOX(7),HCLOSS(7),COLOSS(7),NOXLOSS(7),IMPN_NP(2,25), 
     *RSDCSE(4,3,2),IMONLY_NP(7) 
