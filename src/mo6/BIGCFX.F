      SUBROUTINE BIGCFX(IDX,MY,IV,IP,IH) 
C 
C  BIGCFX computes correction factors applied to composite emission 
C  factors, e.g., speed, optional (humidity, trailer towing, a/c), 
C  temperature, by-bag RVP, operating mode, methane and tampering. 
C 
C  Called by EFCALX. 
C 
C  Calls BIGIN1 and BIGTCF. 
C 
C  Changes : (Last change first) 
C 
C  27 Jul 01 AIR Task 15-21: Removed unused IDX from call to BIGTCF. 
C  14 Jun 00 AIR Task 03: Removed calls of BIGIN2, BIGMTH and BIGTAM 
C            and moved them to BIGIN2. These subroutines are not hourly 
C            dependent. 
C  07 Jun 00 AIR Task 03: Integrated the MOBILE5b code into M6. 
C  10 Apr 00 AIR Task 02: Substituted DBSELPOL for PRTFLG 
C  13 Jan 00 @DynTel-ddj 1-013  Removed call to BIGALH. 
C  10 Nov 99 @DynTel-ddj 1-003  Added initialization of OMTCF for 
C            diesel vehicles. 
C  28 Sep 99 @DynTel-ddj for EPA  Removed ISR from call to BIGMET. 
C  1 Mar  99 @DynTel-HQ 2-698 move the calculation of scf to efcalx, 
C                             hccalx and hrloop 
C  21 Apr 99 @DynTel-bg  0-632 Temporarily commented out call to BIGMET 
C   1 Apr 99 @DynTel-DDJ 2-696 Added hourly loop. 
C  30 Jul 98 @DynTel-zk 2-000 Explicit typing 
C 
C  30 Jul 98 @Dyntel -RJD 2-660 Disabled temperature corrections for 
C            vehicles with separated start and running emissions. 
C      
C 
C  Input on call: 
C 
C    Constants:  
C    argument list: IDX,JDX,MY,IV,IP,IH 
C    common blocks: 
C    /IVTYPE/ VM6M5B 
C    /MYRCAL/ TF 
C 
C  Output on return: 
C 
C    common blocks: 
C    /TEMPC4/ LOWCO 
C 
C  Notes: 
C 
C  BIGCFX was modified for MOBILE4 In-House Version 02 to skip 
C  calculation of correction factors when because travel fraction for 
C  the IV x IDX case is 0, the CFs will not be used.  SALHCF for such 
C  cases will print 0.0 in the BYMY table. 
C 
C  10 July 98 @Dyntel -MMS 2-666 inseted VVSTARTS(IV) to seperate 
C            computation of start and running modes 
C  15 Apr 98 @DynTel-RJD 2-654 Added ISR loop variable to distinguish 
C            between start and running emissions. BIGSAL, BIGMET and 
C            BIGOMT are called inside ISR loops. 
C  Oct-21-97 DynTel@gjr 2-645 Removed ICY from call to BIGTCF. 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified 
C            the code to use the symbolic vehicle type names and vehicle 
C            type vectors that were introduced in Mobile6. 
C 
C 
      USE DATABASE, ONLY: DBSELPOL 
      IMPLICIT NONE  
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'MYRCAL.I' 
      INCLUDE 'TEMPC4.I' 
      INCLUDE 'IVTYPE.I' 
C 
C     Declare parameter list 
C 
      INTEGER, INTENT(IN) :: IDX 
      INTEGER, INTENT(IN) :: MY 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: IP 
      INTEGER, INTENT(IN) :: IH 
C 
C     Declare local variables 
C 
      INTEGER IVTMP5 
C 
C  Initialize variables to be used in upcoming loops. 
C 
      CALL BIGIN1(IH) 
C 
C  Skip calculations if travel fraction is 0.0, zeroing emission factor, 
C  or if the pollutant is not selected. 
C 
      IF(TF(IDX,IV).EQ.0.0) GOTO 999 
C 
      IF(DBSELPOL(IP).EQ.1) GOTO 999 
C 
C  Reset CO offset switch and calculate temperature correction factor (TCF). 
C  There is no TCF in MOBILE4 for diesel IV (IVTNP5=0).  When the ambient 
C  temperature is the FTP standard 75 degrees (MDLOHI=0), TCF represents the 
C  RVP CF for RVP MYG cases 2 and 3 for LDGV/T & HDGV, otherwise it remains the 
C  initialized default neutral value (1.0) for all IV (1-8). 
C 
      LOWCO=1 
      IVTMP5=VM6M5B(IV) 
      IF(IVTMP5.GE.1) CALL BIGTCF(IP,MY,IV,IVTMP5) 
C 
  999 RETURN 
      END 
