      SUBROUTINE REINIT 
      USE USERDAT 
C 
C  REINIT re-initializes common block variables and arrays to their 
C  compile time values.  These variables and arrays initial values are 
C  zeroes and ones and easily coded.  They are not BDSAVEd and then 
C  RESTORed, but they do need to be reset at the outset of each call to 
C  MOBILE. 
C 
C  Called by MOBILE. 
C 
C  Calls INITAP 
C 
C  Changes: (Last change first) 
C 
C  10 May 01 AIR Task 11/5: added initialization of IDLEPMFLG 
C  15 May 01 @EPA-BAG Bug120, Removed SSFLAG and SSONLY. 
C  06 Oct 00 @EPA-BAG BUG155 Remove OBDPC 
C  06 Oct 00 @EPA-BG  BUG155 removed NEWETP flag previously used for user input 
C            of ETP phase in 
C  26 Sep 00 @EPA-EG      Added reinitialization of HD conversion factor arrays. 
C  15 Sep 00 AIR Task 08: Added initialization of GLLRDC, EVFBEN, EVSAW. 
C            Changed I/M initializations to depend on MAXIMPGM. Changed 
C            MIL to WMIL for easier variable tracking. 
C  31-Jul-00 @EPA bug185  changed WAVPRE81 and WAV81PLUS from 0.2 to 0.05 
C  08 Jun 00 AIR Task 03: Removed nonblock DO constructs. 
C  15 Jun 00 @EPA-ddj Bug149 Changed initialization of MOFFMTH from 1.0 
C            to 0.0. 
C  16 May 00 @EPA      Removed initialization of NODESC 
C  01 May 00 @EPA      Changed OBDIM to 0.99 from 0.9 according to decision 
C                      made in TRPD group meeting 
C  01 May 00 @EPA-elg  Added TAMBAS reinitialization 
C  25 Apr 00 @EPA-elg  Added MOFFMTH Array and NGV Reinitialization. 
C  10 Apr 00 AIR Task 02: Added initialization of EXPEXH and NODESC. 
C  13 Apr 00 @Dyntel-BG 1-014 Removed call to INIGLL 
C  16 Nov 99 @Dyntel-nh 1-014 Initialize Gross Leaker status 
C  01 Nov 99 @Dyntel-nh 1-004 Initialize DEFEAT status. 
C  28 Sep 99 @DynTel-ddj for EPA  Added ISTD index to OFFMTH. 
C  24 Aug 99 @DynTel-MLA 2-684  Initialize the database flags. 
C   5 Aug 99 @DynTel-ddj 2-699  Removed initializes of arrays BYTAM. 
C  21 Jul 99 @Dyntel=NH 2-660 Call to reinitialize 1981-1993 user input 
C  23 Mar 99 @Dyntel-HXQ 2-698  Initialize speed correction related 
C            coefficient 
C  22 Feb 99 @DynTel-ddj 2-696 Added hourly index to BYTAM. 
C  10 Mar 99 @Dyntel-RJD 2-694  Reinitialized new optional user-input 
C            for OBD. 
C  2  Mar 99 @DynTel-BAG 2-664  Expanded TLVMTU array for hourly 
C            running losses (weekly and hourly indices) 
C  16 Feb 99 @Dyntel-RJD 2-663  Reinitialized new optional user-input 
C            for 81-93 model year LDGV/T I/M programs. 
C  23 Dec 1998 @DynTel-ZK 2-000 Explicit Typing 
C  20 Nov 98 @DynTel-MLA 2-663  Added initialization for TTCDAT. 
C   6 Nov 98 @DynTel-MLA 2-663  A new dimension was added to ICYIM to 
C            allow the array to hold program end dates as well as start 
C            dates. All references to ICYIM in existing code were 
C            updated accordingly. 
C 
C  Input on call: 
C 
C    Constants:  
C    MAXVEH from include file IVTYPE.I 
C    NEYR 
C 
C    common blocks: 
C    /IMPAR5/ DCNTSV 
C    /MAXIMA/ MAXPOL, MAXYRS 
C    /PPEFFD/ PPEFF1 
C    /CHKLST/ MAXIMPGM 
C 
C  Output on return: 
C 
C    common blocks: 
C    /ATPAR1/ LAPSY 
C    /ATPAR2/ ATPPGM, ATPFQT, CRATP, DISTYP 
C    /BASEQ6/ NUMERU 
C    /BYMYC2/ BYBEF4, BYFER 
C    /BYMYC3/ BYEVAP, BYRUNL, BYREFL, BYRSTL 
C    /ETPPAR/ NETPYR 
C    /EVAL/   MEVAL 
C    /EVPDAT/ EVPAER 
C    /IMPAR1/ ICYIM, CRIM, WAVPRE81, WAV81PLUS, EXEMPTAGE, IMGRACE 
C    /IMPAR2/ ILDT 
C    /IMPAR4/ CUTHC, CUTCO, CUTNO 
C    /IMPAR5/ DISCNT 
C    /IMPAR6/ INTYP 
C    /IM240P/ DSIZE 
C    /LDGIMC/ WAVRDC 
C    /LDGOBD/ WMIL, OBDF(3), OBDIM, OBDREP 
C    /LDGSIM/ COIM 
C    /MYRSAV/ NEWCUM 
C    /NGVDAT2/ NGV_NZML, NGV_NDET, NGV_HIGH, 
C              NGV_HDZML, NGV_HDDET, NGV_PC 
C    /OFFSET/ OFFMTH, MOFFMTH 
C    /OPCNTL/ OPTABL 
C    /QUITXQ/ N1QUIT 
C    /PPEFFD/ PPEFF, PPDES 
C    /REGION/ INITPR 
C    /RUNLS1/ TLVMTU 
C    /TAMOUT/ TAMBAG,THS,TDU,TCC,TAMBAS 
C    /TAMOU2/ TOB,GCONLY 
C    /TAMPB1/ PBBTR 
C    /TTCINF/ TTCDAT 
C    /IDLCOM/ IDLEPMFLG 
C 
C 
C  Notes: 
C 
C  REINIT was added for Version 01's conversion of MOBILE4 to a 
C  subroutine.  REINIT was modified in Version 02 to reinitialize the 
C  by model year emission factor components arrays. 
C  REINIT was modified in MOBILE4.1v6 to zero out PBBTR. 
C  REINIT was modified in MOBILE5.3 to initialize AER in INITAP. 
C  24-August-1994 @CSC-pme request 2-446. Third dimension of OFFMTH 
C            array was  removed. 
C  May-1995 @ CSC-zp Req. 2-495 "Add MOB5a-H to MOBILE5b" 
C            Re-initialiize default values of the pressure/purge test 
C            effectiveness and descriptive record. 
C  1 Feb 96 @DynTel-MLA 2-611  Reinitialize DISCNT to default values. 
C  7 May 96 @DynTel-MLA 2-620 Reinitialize IREAD to permit reading new 
C           I/M credit files in second and subsequent runs in batch 
C           mode. 
C  5 Jun 96 @DynTel-MLA 2-620 Removed initialization of IREAD. Default 
C           I/M credits can now be carried over from one run to another. 
C 25 Feb 97 MLA  Added reinitialization of SSFLAG and OPTABL (from 
C           OPCNTL.I) for Mobile6 transition code. Added IVTYPE.I for 
C           Mobile6. 
C 02 Mar 98 @DynTel-rjd 2-204 Reintialized MEVAL for lead phase out 
C           calculations. 
C 20 Apr 98 @DynTel-BAG 2-661 Added SSONLY and LEVNAM initialization to 
C           .FALSE. and to '  ', respectively 
C 21 Apr 98 @DynTel-RJD 2-654 Added an extra dimension to the OFFMTH 
C           and BYBEF4 arrays to distinguish between start and running 
C           emissions. 
C 20 Apr 98 @DynTel-BAG 2-000 Removed BYIMCR and BYIMRE initializations  
C 30 Jun 98 @DynTel-MLA 2-000 Removed initialization for LEVNAM. This    
C           variable has been renamed as P94AR an is now initialized 
C           in INIFLG. 
C 23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C           certification naming conventions. 
C 
      USE DATABASE 
C 
      IMPLICIT NONE 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'ATPAR1.I' 
      INCLUDE 'ATPAR2.I' 
      INCLUDE 'BASEQ6.I' 
      INCLUDE 'BYMYC2.I' 
      INCLUDE 'BYMYC3.I' 
      INCLUDE 'ETPPAR.I' 
      INCLUDE 'EVAL.I' 
      INCLUDE 'EVPDAT.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'HDCCOM.I' 
      INCLUDE 'IMPAR1.I' 
      INCLUDE 'IMPAR2.I' 
      INCLUDE 'IMPAR4.I' 
      INCLUDE 'IMPAR5.I' 
      INCLUDE 'IMPAR6.I' 
      INCLUDE 'IM240P.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'LDGIMC.I' 
      INCLUDE 'LDGOBD.I' 
      INCLUDE 'LDGSIM.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'MYRSAV.I' 
      INCLUDE 'NGV.I' 
      INCLUDE 'OFFSET.I' 
      INCLUDE 'OPCNTL.I' 
      INCLUDE 'PARAME.I' 
      INCLUDE 'QUITXQ.I' 
      INCLUDE 'PPEFFD.I' 
      INCLUDE 'REGION.I' 
      INCLUDE 'RUNLS1.I' 
      INCLUDE 'TAMOUT.I' 
      INCLUDE 'TAMOU2.I' 
      INCLUDE 'TAMPB1.I' 
      INCLUDE 'TTCINF.I' 
      INCLUDE 'CHKLST.I' 
      INCLUDE 'EVIM.I' 
      INCLUDE 'IDLE.I' 
C 
      INTEGER IMPGM 
      INTEGER I 
      INTEGER J 
      INTEGER K 
C 
C  /ATPAR1/: Signify no ATP program by setting ATP start year to 2051. 
C  /ATPAR2/: Reset all ATP parameters. 
C 
      LAPSY=NEYR 
C 
      CALL DEALLOC_8193 !Deallocate 1981-1993 User input 
      CALL INIDEFEAT 
C 
      ATPPGM =   1 
      ATPFQT =   1 
      CRATP  = 100. 
      DISTYP =   1 
C 
C  /BASEQ6/: Re-initializing NUMERU to zeroes makes inaccessible any user EFs 
C            entered during previous calls. 
C 
      NUMERU=0 
C 
C   DATABASE (MODULE) 
C     Initialize a checklist of the items entered by the user at the 
C     run level. This data will be used to form the database output. 
C 
      DBRUNFLAGS = 'N' 
      DBRUNFILES = ' ' 
C 
C  Reset the TTC credits to the default values for all I/M programs. 
C 
C  /TTCINF/ 
C 
      DO IMPGM=1,MAXIMPGM 
        DO I=1,18 
          DO J=1,3 
            DO K=1,24 
               TTCDAT(K,J,I,IMPGM)=TTCDEF(K,J,I) 
            END DO 
          END DO 
        END DO 
      END DO 
C 
C  Zero out prior call's last scenario's by model year ef component results. 
C  Upcoming run may bypass scenario reinit in EF/HCCALX via REG/VMTMIX check. 
C 
C  /BYMYC2/,/TAMOUT/: 
C 
      BYBEF4 = 0.0 
      BYFER  = 0.0 
C 
      TAMBAG = 0.0 
      TAMBAS = 0.0 
C 
C  /BYMYC3/,/TAMOUT/,/TAMOU2/: 
C 
      BYEVAP = 0.0 
      BYRSTL = 0.0 
C 
      BYRUNL = 0.0 
      BYREFL = 0.0 
      GCONLY = 0.0 
      TCC = 0.0 
      TDU = 0.0 
      TOB = 0.0 
      THS = 0.0 
C 
C  /EVAL/ MEVAL 
C 
C  MEVAL is the integer month index. It is used in lead phase out calculations 
C  before it is reset to the calander month for emission factor calculations.  
C 
      MEVAL = 1 
C 
C  /ETPPAR/: Reset number of phase-in years. 
C 
      NETPYR=4 
C 
C  /IM240P/: Initialize DSIZE to zero. 
C 
      DSIZE = 0.0 
C 
C  /MYRSAV/: Insure CUMMIL is calculated. 
C 
      NEWCUM=1 
C 
C  /OFFSET/: Zero out methane offset, in case subsequent call asks for Total HC. 
C 
      OFFMTH  = 0.0 
      MOFFMTH = 0.0 
C 
C  /OPCNTL/ 
C 
      OPTABL = 1 
C 
      EXPEXH = .FALSE. 
C 
C  /IDLCOM/ 
C 
      IDLEPMFLG = .FALSE. 
C 
C  /QUITXQ/ Reset QUITER 1st call switch. 
C 
      N1QUIT=1 
C 
C  /PPEFFD/ 
C    Initialize the default values of the pressure/purge  
C    descriptive record 
C 
      PPDES = ' ' 
C 
      PPEFF = PPEFF1 
C 
C  /REGION/: Turn on "print run header" switch. 
C 
      INITPR=1 
C 
C  /RUNLS1/: Insure that TLVMTU can be reset to default, if necessary 
C 
      TLVMTU(1,1,1)=-1 
C 
C  /EVPDAT/,/IMPAR1/,/IMPAR2/,/IMPAR5/,/IMPAR6/: Remaining clean up items. 
C  Signify no I/M program by setting ICYIM (start year) to 2051. 
C  Reset I/M parameters used in FAIL routine. 
C  ICYIM, CRIM, ILDT, INTYP were expanded to two dimensions for 
C  IM240 MOBILE5 expansion program Jan 92. 
C  DISCNT expanded to two dimensions, reinitialized from DCNTSV. 
C 
      EVPAER = 1.0 
      DISCNT = DCNTSV 
C 
      CALL INITAP 
C 
C  Reset default values of IM parameters that the user can optionally 
C  override. Compliance rate is set at 85% or the non-compliance is 
C  set to 15%(see p.4 M6.IM.001). Cost waivered percent of vehicles 
C  have a default value of 20%(see p.4 M6.IM.001). 
C 
C  /IMPAR1/ 
      EXEMPTAGE = MAXYRS 
      IMGRACE   = 1 
      CRIM      = 85.0 
      WAVPRE81  = 0.05 
      WAV81PLUS = 0.05 
      EVFBEN    = 1.0 
C  /LDGIMC/ 
      WAVRDC    = 0.80 
C  /LDGSIM/ 
      COIM      = 16.0 
C  /IMPAR2/ 
      ILDT      = 1 
C  /IMPAR4/ 
      CUTHC     =  1.2 
      CUTCO     = 20.0 
      CUTNO     =  3.0 
C 
C     INTYP (test only or test & repair) is a required I/M input, so 
C     run level reinitialization is not required 
C  /IMPAR6/ 
      INTYP(1)  = 1 
      INTYP(2)  = 1 
C 
C  /LDGOBD/ 
      WMIL     = 0.85   !see p.13 of M6.EXH.009 
      OBDF(1) = 0.9    !0 < mileage < 36,001       (see p.14 of M6.EXH.009 and BD_IM) 
      OBDF(2) = 0.1    !36,000 < mileage < 80,001  (see p.14 of M6.EXH.009) 
      OBDF(3) = 0.0    !80,000 < mileage           (see p.14 of M6.EXH.009) 
      OBDIM   = 0.99   !OBD I&M effectiveness (see p.15&16 of M6.EXH.009) 
      OBDREP  = 1.5    !emissions after OBD induced repairs = 1.5*CRT050k (see p.14 of M6.EXH.009) 
C 
C  /EVIM/ 
      GLLRDC = 0.0     !gross liquid leaker reductions 
C 
C  /TAMPB1/: Zero out misfueling categories tampering rates save area. 
C 
      PBBTR = 0.0 
C 
C Determine SCFF and SCFAC, OCECO, OCLEVEL 
C For Speed Correction 
C 
      SCFUSER = 1 
      CALL SCFCO() 
C 
C     NGV EF reinitialization. 
C 
      NGVFLG1 = 0 
      NGVFLG2 = 0 
      NGV_START = 1994 
C 
      NGVPHIN = 0.0 
C 
      NGVZML  = -1.0 
      NGVDET  = -1.0 
      NGVHIGH = -1.0 
C 
      HDCFAC_FLG = 0 
      HDCFAC2 = HDCFAC 
C 
      RETURN 
      END 
