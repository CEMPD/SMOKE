      SUBROUTINE BAGEME(IDX,ITTT) 
C 
C  BAGEME multiplies the effects group sizes times the corresponding 
C  weighted emission impact rates to get the tampering bag emission effects. 
C 
C  Called by TAMPER. 
C 
C  Changes : 
C 
C  04 Dec 00 @EPA-bag Task X6, removed MY, JDX, and IDX from MYCODE.I; 
C            pass IDX through parameter list 
C  24 Jan 00 EPA-elg  Added ITTT to subroutine parameter list. 
C 
C  Input on call: 
C 
C    common blocks: 
C    /EGSCAL/ EGS 
C    /EMECAL/ EMI 
C    /LOOKUP/ IVTAM 
C    /MAXIMA/ MAXPOL 
C 
C  Output on return: 
C 
C    common blocks: 
C    /LOOKUP/ IPG 
C    /TAMOUT/ TAMBAG 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  EMSUM    R    tampering emission effects sum for given IVTAM, IDX, IP & IB 
C 
C  Notes: 
C 
C  None. 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'EGSCAL.I' 
      INCLUDE 'EMECAL.I' 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'LOOKUP.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'MYCODE.I' 
      INCLUDE 'TAMOUT.I' 
C 
      INTEGER, INTENT(IN) :: IDX 
      INTEGER IP, IEG, IB, ITTT 
      REAL EMSUM        
      DO 30 IP=1,MAXPOL 
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
      IF(ITTT.EQ.1) THEN 
        TAMBAG(IB,IP,IDX,IVTAM)=EMSUM 
      ELSEIF(ITTT.EQ.2) THEN 
        TAMBAS(IB,IP,IDX,IVTAM)=EMSUM 
      ENDIF 
C 
   20 CONTINUE 
   30 CONTINUE 
C 
      RETURN 
      END 
