C Version @(#)$Id$ $Source$ $Date$ 

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
C    common blocks:
C    /MAXIMA/ MAXVEH,MAXPOL,MAXYRS
C
C  Output on return:
C
C    common blocks:
C    /ATPAR1/ LAPSY
C    /ATPAR2/ ATPPGM,ATPFQT,CRATP,DISTYP
C    /BASEQ6/ NUMERU
C    /BASEQ7/ IPMOD
C    /BYMYC2/ BYBEF4,BYTAM,BYFER
C    /BYMYC3/ BYEVAP,BYRUNL,BYREFL,BYRSTL
C    /BYMYC5/ BYIMCR,BYIMRE
C    /ETPPAR/ NEWETP,NETPYR
C    /EVPDAT/ EVPAER
C    /FUELS2/ IMTIER
C    /IMPAR1/ ICYIM,CRIM
C    /IMPAR2/ ILDT
C    /IMPAR6/ INTYP
C    /IM240P/ DSIZE
C    /MYRSAV/ NEWCUM
C    /OFFSET/ OFFMTH
C    /QUITXQ/ N1QUIT
C    /REGION/ INITPR
C    /RUNLS1/ TLVMTU
C    /TAMID2/ IEK
C    /TAMID3/ OFFIDL
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
C
C
      REAL JULMYR
      INTEGER ATPPGM,ATPFQT,DISTYP
C
      COMMON /ATPAR1/ LAPSY,LAP1ST,LAPLST,LVTFLG(4)
      COMMON /ATPAR2/ ATPPGM,ATPFQT,CRATP,DISTYP(8)
      COMMON /BASEQ6/ MYGERU(12,2,3,8,2),MAXERU,NUMERU(3,8,2),KEYER,IGER
      COMMON /BASEQ7/ IPMOD(3)
      COMMON /BYMYC2/ BYBEF4(3,25,8),BYTAM(3,25,4),BYFER(3,25,8)
      COMMON /BYMYC3/BYEVAP(25,8),BYRUNL(25,4),BYREFL(25,4),BYRSTL(25,8)
      COMMON /BYMYC5/ BYIMCR(3,25,4),BYIMRE(3,25,4),REDSUM(3,4)
      COMMON /ETPPAR/ ETP(7,3),NEWETP,IETPYR(4),ETPPIV(4),NETPYR
      COMMON /EVPDAT/ EVPAER(3,2),EVPTGS(3,2,4),CCEMI(7,4)
      COMMON /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
      COMMON /IMPAR1/ ICYIM(2),ISTRIN(2),MODYR(2,2),WAIVER(2,2),CRIM(2)
      COMMON /IMPAR2/ ILDT(4,2),ITEST(2)
      COMMON /IMPAR6/ IFREQ(2),INTYP(2),IFILET(2)
      COMMON /IM240P/ DSIZE(25,4)
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /MYRSAV/ AMAR(25,8),JULMYR(25,8),NEWCUM
      COMMON /OFFSET/ OFFCO(25,3),OFFMTH(25,8,2)
      COMMON /QUITXQ/ N1QUIT
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
      COMMON /RUNLS1/ TLEMI(6,3,4,4,4,4),TLVMT(6),TLVMTU(6),MAXTL
      COMMON /TAMID2/ MYIDLE(3,3),IEK
      COMMON /TAMID3/ OFFIDL(3,25,4)
      COMMON /TAMOUT/ TAMBAG(3,3,25,4),THS(2,25,4),TDU(25,4),TCC(25,4)
      COMMON /TAMOU2/ TOB(25,4),GCONLY(25,4)
      COMMON /TAMPB1/ LPOD,PBBTR(3,2,25,4)
C
C  /ATPAR1/: Signify no ATP program by setting ATP start year to 2020.
C  /ATPAR2/: Reset all ATP parameters.
C
      LAPSY=2020
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
C  /BASEQ7/: IPMOD must be reset, else idle calulation cases may be skipped due
C            to prior calls' user supplied EFs.
C
      DO 20 IP=1,MAXPOL
      IPMOD(IP)=0
   20 CONTINUE
C
C  Zero out prior call's last scenario's by model year ef component results.
C  Upcoming run may bypass scenario reinit in EF/HCCALX via REG/VMTMIX check.
C
C  /BYMYC2/,/TAMID3/,/TAMOUT/:
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
      OFFIDL(IP,JDX,IV)=0.0
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
      DO 80 IV=1,MAXVEH
      DO 80 ICAP=1,2
      OFFMTH(JDX,IV,ICAP)=0.0
   80 CONTINUE
C
C  /QUITXQ: Reset QUITER 1st call switch.
C
      N1QUIT=1
C
C  /REGION/: Turn on "print run header" switch.
C
      INITPR=1
C
C  /RUNLS1/: Insure that TLVMTU can be reset to default, if necessary
C
      TLVMTU(1)=-1
C
C  /EVPDAT/,/IMPAR1/,/IMPAR2/,/IMPAR6/,/TAMID2/:  Remaining
C  clean up items.
C  Signify no I/M program by setting ICYIM (start year) to 2020.
C  Reset I/M parameters used in FAIL routine.
C  ICYIM, CRIM, ILDT, INTYP were expanded to two dimensions for
C  IM240 MOBILE5 expansion program Jan 92.
C
      DO 85 IAY=1,2
      DO 85 IHG=1,2
        EVPAER(IHG,IAY)=1.0
   85 CONTINUE
      CALL INITAP
      ICYIM(1)=2020
      ICYIM(2)=2020
      CRIM(1)=0.0
      CRIM(2)=0.0
      DO 90 IV=1,IVTAM
      ILDT(IV,1)=1
      ILDT(IV,2)=1
   90 CONTINUE
      IEK=1
      INTYP(1)=1
      INTYP(2)=1
C
C  /TAMPB1/: Zero out misfueling categories tampering rates save area.
C
      DO 95 IVTAM=1,4
      DO 95 IDX=1,MAXYRS
      DO 95 IAY=1,2
      DO 95 IDPB=1,3
      PBBTR(IDPB,IAY,IDX,IVTAM)=0.0
  95  CONTINUE
C
      RETURN
      END
