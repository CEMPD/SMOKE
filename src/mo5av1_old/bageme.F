C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE BAGEME
C
C  BAGEME multiplies the effects group sizes times the corresponding
C  weighted emission impact rates to get the tampering bag emission effects
C  and idle tampering offsets.
C
C  Called by TAMPER.
C
C  Input on call:
C
C    common blocks:
C    /EGSCAL/ EGS
C    /EMECAL/ EMI
C    /FLAGS4/ IDLFLG
C    /LOOKUP/ IVTAM
C    /MAXIMA/ MAXPOL
C    /MYCODE/ IDX
C    /TAMID1/ EMIDL
C
C  Output on return:
C
C    common blocks:
C    /LOOKUP/ IPG
C    /TAMID2/ IEK
C    /TAMID3/ OFFIDL
C    /TAMOUT/ TAMBAG
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  EMSUM    R    tampering emission effects sum for given IVTAM, IDX, IP & IB
C  EISUM    R    idle tampering emission effects sum for given IDX & IP
C  Notes:
C
C  None.
C
C
      INTEGER PRTFLG,HCFLAG,COLDFG
C
      COMMON /EGSCAL/ AER(11,11,2,2),TGS(11,6,4),EGS(7,2)
      COMMON /EMECAL/ EM3W(3,3,3),EMOX(3,3,3),EMEGR(4,3,3),EMI(7,3,3)
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /LOOKUP/ IVTAM,IQG,IPG,JPGD,IHG,IGCSF
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /MYCODE/ MY,IDX,JDX,LDXSY,LMYRVT,IAY,IMDXSY,IMKINK
      COMMON /TAMID1/ EMI3W(3,3,3),EMIOX(3,3,3),EMIDL(3,4,3)
      COMMON /TAMID2/ MYIDLE(3,3),IEK
      COMMON /TAMID3/ OFFIDL(3,25,4)
      COMMON /TAMOUT/ TAMBAG(3,3,25,4),THS(2,25,4),TDU(25,4),TCC(25,4)
C
      DO 40 IP=1,MAXPOL
C
C  Effects group sizes are identical for HC and CO, but the emissions
C  impacts are not.
C
      IPG=1
      IF(IP.EQ.3) IPG=2
C
      DO 20 IB=1,3
      EMSUM=0.0
C
      DO 10 IEG=1,7
      EMSUM=EMSUM+EGS(IEG,IPG)*EMI(IEG,IB,IP)
   10 CONTINUE
C
      TAMBAG(IB,IP,IDX,IVTAM)=EMSUM
   20 CONTINUE
C
C  If idle ef requested, do idle tampering offset for current IP & IDX case.
C
      IF(IDLFLG.EQ.1) GOTO 40
      EISUM=0.0
C
      DO 30 IEG=1,4
      IEK=IEG-1
      IF(IEG.EQ.1) IEK=1
      IGTIDL=ITAMPT(5)
      EISUM=EISUM+EGS(IEG,IPG)*EMIDL(IGTIDL,IEG,IP)
   30 CONTINUE
C
      OFFIDL(IP,IDX,IVTAM)=EISUM
   40 CONTINUE
C
      RETURN
      END
