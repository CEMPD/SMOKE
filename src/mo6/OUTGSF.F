      SUBROUTINE OUTGSF(ICY) 
C 
C  OUTGSF prints the user supplied diesel sales fractions (1.0-GSFRAC) 
C  for LDV and LDT. 
C 
C  Called by OUTHD3. 
C 
C  Changes: (Last change first) 
C 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C  03 Dec 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: ICY 
C    common blocks: 
C    /IOUCOM/ IOUREP 
C    /MAXIMA/ MAXYRS 
C    /REGISF/ USRGSF 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  NAMVEH   C*4  the two classes for which user inputs diesel fraction 
C  USRDSL   R    used to convert USRGSF to diesel fraction for output 
C  USC      C*4  underscore strings for table and column headings 
C 
C  Notes: 
C 
C  Echos input for DSFLAG=2 case.  DSFLAG was added for MOBILE4.1 vD06 
C 
      IMPLICIT NONE 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'REGISF.I' 
      INTEGER, INTENT(IN) :: ICY 
      INTEGER MY, IG, MYL,  MYF, IVTEMP, MAXIDX, JDX, I, IDX 
      REAL, DIMENSION(10) :: USRDSL 
      CHARACTER(LEN=4) :: USC4 = '____' 
      CHARACTER(LEN=4), DIMENSION(2) :: NAMVEH=(/'LDDV','LDDT'/) 
C 
      WRITE(IOUREP,200) 
  200 FORMAT('0', 
     *'Replacement Diesel Sales Fractions Input by User:') 
C 
      MY=ICY-MAXYRS+1 
      DO 30 IG=1,3 
      MYF=10*(IG-1)+1 
      MYL=MIN0(MYF+9,MAXYRS) 
      WRITE(IOUREP,210) (I,I=MY,MY+MYL-MYF) 
  210 FORMAT('0',7X,10(2X,I4,1X)) 
      WRITE(IOUREP,220) (USC4,I=MYF,MYL) 
  220 FORMAT('+',7X,10(2X,A4,1X)) 
      DO 20 IVTEMP=1,2 
      MAXIDX=MYL-MYF+1 
      DO 10 IDX=1,MAXIDX 
      JDX=MAXYRS-IDX+1 
      USRDSL(IDX)=1.0-USRGSF(JDX-MYF+1,IVTEMP) 
   10 CONTINUE 
      WRITE(IOUREP,230) NAMVEH(IVTEMP),(USRDSL(IDX),IDX=1,MAXIDX) 
  230 FORMAT(' ',A4,':',2X,10F7.3) 
   20 CONTINUE 
      MY=MY+10 
   30 CONTINUE 
C 
      RETURN 
      END 
