      SUBROUTINE ATHCHK(ICY,INERR) 
C 
C  This routine verifies the HAPS data and puts the values into 
C  the ATHBER array. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX5.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
C 
      INTEGER IATH 
      INTEGER ICY 
      INTEGER IMAP 
      INTEGER INERR 
      INTEGER IRB 
      INTEGER IV 
      INTEGER IXV 
      INTEGER JDX 
      INTEGER MY 
      INTEGER NJDX 
      INTEGER NMAP 
C 
C  Initialize the holding array to the missing code 
C 
      ATHBER = -999. 
C 
C  Initialize the pollutant mapping array 
C 
      NMAP=0 
      ATHPMAP = 0 
C 
C  Loop for each input ATH 
C 
      DO IATH=1,NATH 
C 
C  Map the pollutant number 
C 
      IF(NMAP.EQ.0) THEN 
      GOTO 5 
      ELSE 
      DO IMAP=1,NMAP 
      IF(ATHPOL(IATH).EQ.ATHPMAP(IMAP)) GOTO 10 
      END DO 
      END IF 
C 
   5  NMAP=NMAP+1 
      IMAP=NMAP 
      ATHPMAP(IMAP)=ATHPOL(IATH) 
C 
C  Set the vehicle class 
C 
   10 IV = ATHVTYPE(IATH) 
C 
C  Set the exhaust/evap index 
C 
      IF(ATHEKIND(IATH).EQ.'EXH') THEN 
      IXV = 1 
      ELSE 
      IXV = 2 
      END IF 
C 
C  Set the ratio/ber index 
C 
      IF(ATHRATIO(IATH).EQ.'RATIOVOC') IRB=1 
      IF(ATHRATIO(IATH).EQ.'RATIOTOG') IRB=2 
      IF(ATHRATIO(IATH).EQ.'BER')      IRB=3 
      IF(ATHRATIO(IATH).EQ.'RATIOPM')  IRB=4 
C 
C  Process the model years. If the model years are outside the 
C  25 year band, skip. 
C 
      DO MY=ATHBEGMY(IATH),ATHENDMY(IATH) 
C 
      JDX=ICY-MY+1 
C 
      IF(JDX.LT.1 .OR. JDX.GT.25) CYCLE 
C 
C  Load the array, checking for overlaps 
C 
      IF(ATHBER(1,JDX,IXV,IRB,IV,IMAP).EQ.-999. .AND. 
     *   ATHBER(2,JDX,IXV,IRB,IV,IMAP).EQ.-999.) THEN 
      ATHBER(1,JDX,IXV,IRB,IV,IMAP)=ATHRATZML(IATH) 
      ATHBER(2,JDX,IXV,IRB,IV,IMAP)=ATHDET(IATH) 
      ELSE 
      WRITE(IOUERR,2) MY,TRIM(ATHPNAME(IATH)) 
    2 FORMAT(' Error: Conflicting HAPS data given for MY ',I4,/ 
     *       ' and pollutant ',A) 
      GOTO 999 
      END IF 
C 
      END DO ! MY 
C 
      END DO ! IATH 
C 
C  With the array now loaded, check to be sure all of the 
C  model years for a given combination are filled. 
C 
      DO IMAP=1,NMAP 
      DO IV=1,MAXVEH 
      DO IRB=1,4 
      DO IXV=1,2 
      NJDX=0 
      DO JDX=1,25 
      IF(ATHBER(1,JDX,IXV,IRB,IV,IMAP).NE.-999.) NJDX=NJDX+1 
      END DO 
      IF(NJDX.NE.0 .AND. NJDX.NE.25) THEN 
      WRITE(IOUERR,3) IMAP,IV,IRB,IXV 
   3  FORMAT(' Error: All of the model year data for HAPS are not'/ 
     &      ' filled in. IMAP=',I2,' IV=',I2,' IRB=',I1,' IXV=',I1) 
      GO TO 999 
      END IF 
      END DO 
      END DO 
      END DO 
      END DO 
C 
      GOTO 99 
C 
  999 INERR=INERR+1 
C 
   99 RETURN 
      END 
       
       
