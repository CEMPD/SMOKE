      SUBROUTINE SCFTAB(IV,IP) 
C 
C SCFTAB calculate the speed correction factor table for 
C HDGV, LDDV, LDDT, HDDV, BUS. 
C 
C  Changes: (Last change first) 
C  
C  13 Jan 00 @DynTel-ddj  1-013  Placed AVSPR and AVSPL into common 
C            block (SPEED9.I). 
C  17 Aug 99 @DynTel-ddj  2-698  Changed CALL to BIGSC1 to a function. 
C  18 May 1999 @DynTel-HQ 2-698 delete the code to call scfsub 
C   1 Mar 1999 @DynTel-HQ 2-698 calculate SCF table 
C                          using MOBILE5 method. 
C 
C    Called By EFCALX, HCCALX 
C 
C    Calls BIGSC1 
C 
C    argument list: IV IP 
C    common blocks: 
C    /IVPCOM/    IVPTRL,IVPTRA,IVPTRB 
C    /SPEED4/    IVB 
C    /SPEED6/    SALHCF 
C    /SPEED9/    M5SCF 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'SPEED6.I' 
      INCLUDE 'SPEED4.I' 
      INCLUDE 'IVPCOM.I' 
      INCLUDE 'SPEED9.I' 
C 
      INTEGER IP,IV,I 
      REAL BIGSC1 
      REAL SP 
C 
      IVB=IVPTRB(IV) 
C     
C calculate the SCF for local road at 12.9 
C 
      SCFL = BIGSC1(IP,IV,AVSPL) 
C 
C calculate the SCF for RAMP road at 34.9 
C 
      SCFR = BIGSC1(IP,IV,AVSPR) 
C 
C calculate the SCF for Freeway, Arterials and Collectors 
C 
      M5SCF(1,IP) = BIGSC1(IP,IV,2.5) 
C 
      SP=0.0 
      DO I=2,14 
         SP=SP+5.0 
         M5SCF(I,IP) = BIGSC1(IP,IV,SP) 
      END DO 
C 
100   RETURN 
      END 
