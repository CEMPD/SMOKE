C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE REINIT
C
C  REINIT re-initializes common block variables and arrays to their compile
C  time values.  These variables and arrays initial values are zeroes and ones
C  and easily coded.  They are not BDSAVEd and then RESTORed, but they do need
C  to be reset at the outset of each call to MOBILE.
C
C  Called by MOBILE.
C
C  Calls INITAP
C
C  Input on call:
C
C    constant: NEYR
C
C    common blocks:
C    /IMPAR5/ DCNTSV
C    /MAXIMA/ MAXVEH,MAXPOL,MAXYRS
C    /PPEFFD/ PPEFF1
C
C  Output on return:
C
C    common blocks:
C    /ATPAR1/ LAPSY
C    /ATPAR2/ ATPPGM,ATPFQT,CRATP,DISTYP
C    /BASEQ6/ NUMERU
C    /BYMYC2/ BYBEF4,BYTAM,BYFER
C    /BYMYC3/ BYEVAP,BYRUNL,BYREFL,BYRSTL
C    /BYMYC5/ BYIMCR,BYIMRE
C    /ETPPAR/ NEWETP,NETPYR
C    /EVPDAT/ EVPAER
C    /FUELS2/ IMTIER
C    /IMPAR1/ ICYIM,CRIM
C    /IMPAR2/ ILDT
C    /IMPAR5/ DISCNT
C    /IMPAR6/ INTYP
C    /IM240P/ DSIZE
C    /MYRSAV/ NEWCUM
C    /OFFSET/ OFFMTH
C    /QUITXQ/ N1QUIT
C    /PPEFFD/ PPEFF,PPDES,RBHPPF
C    /REGION/ INITPR
C    /RUNLS1/ TLVMTU
C    /TAMOUT/ TAMBAG,THS,TDU,TCC
C    /TAMOU2/ TOB,GCONLY
C    /TAMPB1/ PBBTR
C
C
C  Notes:
C
C  REINIT was added for Version 01's conversion of MOBILE4 to a subroutine.
C  REINIT was modified in Version 02 to reinitialize the by model year
C  emission factor components arrays.
C  REINIT was modified in MOBILE4.1v6 to zero out PBBTR.
C  REINIT was modified in MOBILE5.3 to initialize AER in INITAP.
C  REINIT was modified in MOBILE5v13 to initialize IMTIER.
C  24-August-1994 @CSC-pme request 2-446. Third dimension of OFFMTH array
C                          was  removed.
C
C  May-1995 @ CSC-zp Req. 2-495 "Add MOB5a-H to MOBILE5b"
C  Re-initialiize default values of the pressure/purge  test effectiveness
C  and  descriptive record
C 
C  1 Feb 96 @DynTel-MLA 2-611  Reinitialize DISCNT to default values.
C
C  7 May 96 @DynTel-MLA 2-620 Reinitialize IREAD to permit reading new
C      I/M credit files in second and subsequent runs in batch mode.
C  5 Jun 96 @DynTel-MLA 2-620 Removed initialization of IREAD. Default
C      I/M credits can now be carried over from one run to another.
C
      INCLUDE 'ATPAR1.I'
      INCLUDE 'ATPAR2.I'
      INCLUDE 'BASEQ6.I'
      INCLUDE 'BYMYC2.I'
      INCLUDE 'BYMYC3.I'
      INCLUDE 'BYMYC5.I'
      INCLUDE 'ETPPAR.I'
      INCLUDE 'EVPDAT.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'IMPAR1.I'
      INCLUDE 'IMPAR2.I'
      INCLUDE 'IMPAR5.I'
      INCLUDE 'IMPAR6.I'
      INCLUDE 'IM240P.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'MYRSAV.I'
      INCLUDE 'OFFSET.I'
      INCLUDE 'PARAME.I'
      INCLUDE 'QUITXQ.I'
      INCLUDE 'PPEFFD.I'
      INCLUDE 'REGION.I'
      INCLUDE 'RUNLS1.I'
      INCLUDE 'TAMOUT.I'
      INCLUDE 'TAMOU2.I'
      INCLUDE 'TAMPB1.I'
C
C  /ATPAR1/: Signify no ATP program by setting ATP start year to 2051.
C  /ATPAR2/: Reset all ATP parameters.
C
      LAPSY=NEYR
C
      ATPPGM=1
      ATPFQT=1
      CRATP=100.
      DO 5 IV=1,MAXVEH
      DISTYP(IV)=1
   5  CONTINUE
C
C  /BASEQ6/: Re-initializing NUMERU to zeroes makes inaccessible any user EFs
C            entered during previous calls.
C
      DO 10 IR=1,MAXREG
      DO 10 IV=1,MAXVEH
      DO 10 IP=1,MAXPOL
      NUMERU(IP,IV,IR)=0
   10 CONTINUE
C
C  Zero out prior call's last scenario's by model year ef component results.
C  Upcoming run may bypass scenario reinit in EF/HCCALX via REG/VMTMIX check.
C
C  /BYMYC2/,/TAMOUT/:
C
      IVTAM=4
      DO 30 IV=1,MAXVEH
      DO 30 JDX=1,MAXYRS
      DO 30 IP=1,MAXPOL
      BYBEF4(IP,JDX,IV)=0.0
      BYFER(IP,JDX,IV)=0.0
   30 CONTINUE
C
      DO 40 IV=1,IVTAM
      DO 40 JDX=1,MAXYRS
      DO 40 IP=1,MAXPOL
      BYTAM(IP,JDX,IV)=0.0
      DO 35 IB=1,3
      TAMBAG(IB,IP,JDX,IV)=0.0
   35 CONTINUE
   40 CONTINUE
C
C  /BYMYC3/,/TAMOUT/,/TAMOU2/:
C
      DO 50 IV=1,MAXVEH
      DO 50 JDX=1,MAXYRS
      BYEVAP(JDX,IV)=0.0
      BYRSTL(JDX,IV)=0.0
   50 CONTINUE
C
      DO 60 IV=1,IVTAM
      DO 60 JDX=1,MAXYRS
      BYRUNL(JDX,IV)=0.0
      BYREFL(JDX,IV)=0.0
      GCONLY(JDX,IV)=0.0
      TCC(JDX,IV)=0.0
      TDU(JDX,IV)=0.0
      TOB(JDX,IV)=0.0
      DO 55 IFDS=1,2
      THS(IFDS,JDX,IV)=0.0
   55 CONTINUE
   60 CONTINUE
C
C  /BYMYC5/:
C
C  Note: REDSUM is covered for MOBILE call reinitialization by EF/HCCALX code.
C
      DO 70 IVIM=1,4
      DO 70 JDX=1,MAXYRS
      DO 70 IPIM=1,3
      BYIMCR(IPIM,JDX,IVIM)=0.0
      BYIMRE(IPIM,JDX,IVIM)=0.0
   70 CONTINUE
C
C
C  /FLAGS2/ IMTIER
C
C   IMTIER flag signals the application or absence of the
C          Special I/M for Tier 1 Vehicles
C          (1994-1999+ LDGV & LDGT1/2)
C
      IMTIER = 0
C
C  /ETPPAR/: Reset new ETP data flag and number of phase-in years.
C
      NEWETP=1
      NETPYR=4
C
C  /IM240P/: Initialize DSIZE to zero.
C
      DO 75 IV=1,IVTAM
      DO 75 JDX=1,MAXYRS
      DSIZE(JDX,IV)=0.0
  75  CONTINUE
C
C  /MYRSAV/: Insure CUMMIL is calculated.
C
      NEWCUM=1
C
C  /OFFSET/: Zero out methane offset, in case subsequent call asks for Total HC.
C
      DO 80 JDX=1,MAXYRS
C
        DO 80 IV=1,MAXVEH
          OFFMTH(JDX,IV)=0.0
C
   80 CONTINUE
C
C  /QUITXQ: Reset QUITER 1st call switch.
C
      N1QUIT=1
C
C  /PPEFFD/
C    Initialize the default values of the pressure/purge 
C    test effectiveness and descriptive record
C
      PPDES = ' '
      RBHPPF=1
C
      DO 82 K1=1,2
C
        DO 82 J1=1,25
          PPEFF(J1,K1)=PPEFF1(J1,K1)
C
  82  CONTINUE
C
C  /REGION/: Turn on "print run header" switch.
C
      INITPR=1
C
C  /RUNLS1/: Insure that TLVMTU can be reset to default, if necessary
C
      TLVMTU(1)=-1
C
C  /EVPDAT/,/IMPAR1/,/IMPAR2/,/IMPAR5/,/IMPAR6/: Remaining clean up items.
C  Signify no I/M program by setting ICYIM (start year) to 2051.
C  Reset I/M parameters used in FAIL routine.
C  ICYIM, CRIM, ILDT, INTYP were expanded to two dimensions for
C  IM240 MOBILE5 expansion program Jan 92.
C  DISCNT expanded to two dimensions, reinitialized from DCNTSV.
C
      DO 85 IAY=1,2
C
        DO 85 IHG=1,2
          EVPAER(IHG,IAY)=1.0
C
   85 CONTINUE
C
      DO 87 IPIM = 1, 3
C
        DO 87 IP = 1, 3
          DISCNT(IP, IPIM) = DCNTSV(IP, IPIM)
C
  87  CONTINUE
C
      CALL INITAP
      ICYIM(1)=NEYR
      ICYIM(2)=NEYR
      CRIM(1)=0.0
      CRIM(2)=0.0
C
      DO 90 IV=1,IVTAM
        ILDT(IV,1)=1
        ILDT(IV,2)=1
   90 CONTINUE
C
      INTYP(1)=1
      INTYP(2)=1
C
C  /TAMPB1/: Zero out misfueling categories tampering rates save area.
C
      DO 95 IVTAM=1,4
C
        DO 95 IDX=1,MAXYRS
C
          DO 95 IAY=1,2
C
            DO 95 IDPB=1,3
              PBBTR(IDPB,IAY,IDX,IVTAM)=0.0
C
  95  CONTINUE
C         
      RETURN
      END
