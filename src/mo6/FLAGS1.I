C 
C     The data is initialized in BD17. 
C 
C  Changes: (Last change first) 
C 
C 27 Nov 00 @EPA-bag bug#301 Removed FVMTFLG, HVMTFLG, SVMTFLG 
C  23 Mar 1999 @DynTel-HXQ 2-698 Add flag for Speed Correction 
C  16 Mar 99 @DynTel-BG 2-664 Removed RNLFLG. 
C  29 Sep 98 @DynTel-HQ 2-678  Add CONVFLG for conversion factor 
C  25 Jun 98 @DynTel-MLA 2-659  Changed the PROMPT flag to USR94P, since  
C            its only purpose now is to indicate whether the user has  
C            entered new implementation or EF data for 94+ LDG vehicles. 
C 
C             
C 
C 
      INTEGER TAMFLG,VMFLAG,SPCOFLG,SCFFLG,OXYFLG,DSFLAG,DTGFLG, 
     *  USR94P,GLOXFLAG,GLDSFLAG,GLVMFLG,GLDTGFLG,CONVFLG,SCFUSER 
      COMMON /FLAGS1/ TAMFLG,VMFLAG,SPCOFLG,SCFFLG,OXYFLG,DSFLAG, 
     *                DTGFLG,GLOXFLAG,GLDSFLAG,GLDTGFLG,GLVMFLG, 
     *                USR94P,CONVFLG,SCFUSER 
