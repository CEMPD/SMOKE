      SUBROUTINE ATHCHK(ICY,INERR) 
C 
C  This routine verifies the HAPS data and puts the values into 
C  the ATHBER array. 
C 
C  Changes: 
C  AIR Task 28: Added ATHEAP for AT dynamic array support. 
C 
C  Called by: ATCALX 
C 
C  Calls: None 
C 
C 
      USE ATHEAP 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX5.I'  ! NMAP,ATHPOL,ATHPMAP,ATHSNAME,ATHPNAME,ATHEKIND, 
                         ! ATHRATIO,ATHBER,ATHRATZML,ATHDET,ATHIV 
      INCLUDE 'IOUCOM.I' ! IOUERR 
      INCLUDE 'IVTYPE.I' ! MAXVEH 
C 
C  Declare parameter list 
C 
      INTEGER ICY   ! Calendar year 
      INTEGER INERR ! Error 
C 
C  Declare local variables 
C 
      INTEGER IATH  ! HAPS pollutant index 
      INTEGER IMAP  ! HAPS data mapped into the ATHBER storage array 
      INTEGER IRB   ! Ratio/BER index 1=VOC, 2=TOG, 3=BEF, 4=PM 
      INTEGER IV    ! Vehicle class 
      INTEGER JDX   ! age 
      INTEGER JXV   ! HAPS eXhaust/eVap index 
      INTEGER MY    ! Model year 
      INTEGER NJDX  ! Number of valid ages 
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
   5   NMAP=NMAP+1 
       IMAP=NMAP 
       ATHPMAP(IMAP)=ATHPOL(IATH) 
       ATHSNAME(IMAP)=ATHPNAME(IATH)(1:12) 
C 
C  Set the vehicle class 
C 
   10  IV = ATHVTYPE(IATH) 
C 
C  Set the exhaust/evap index 
C 
       IF(ATHEKIND(IATH).EQ.'EXH') THEN 
        JXV = 1 
       ELSE 
        JXV = 2 
       END IF 
C 
C  Set the ratio/ber index. IRB=1-3 are AT-based, 4 is PM-based 
C 
       IF(ATHRATIO(IATH).EQ.'RATIOVOC') IRB=1 
       IF(ATHRATIO(IATH).EQ.'RATIOTOG') IRB=2 
       IF(ATHRATIO(IATH).EQ.'BEF')      IRB=3 
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
        IF(ATHBER(1,JDX,JXV,IRB,IV,IMAP).EQ.-999. .AND. 
     *     ATHBER(2,JDX,JXV,IRB,IV,IMAP).EQ.-999.) THEN 
         ATHBER(1,JDX,JXV,IRB,IV,IMAP)=ATHRATZML(IATH) 
         ATHBER(2,JDX,JXV,IRB,IV,IMAP)=ATHDET(IATH) 
        ELSE 
         WRITE(IOUERR,2) MY,TRIM(ATHPNAME(IATH)) 
    2    FORMAT(' Error: Conflicting HAPS data given for MY ',I4,/ 
     *         ' and pollutant ',A) 
         GOTO 999 
        END IF 
C 
       END DO ! MY 
C 
      END DO ! IATH 
C 
C  With the array now loaded, check to be sure all of the 
C  model years for a given combination are filled. If so, set 
C  a flag for later use, otherwise send an error if only 
C  partially filled. 
C 
      DO IMAP=1,NMAP 
       DO IV=1,MAXVEH 
        DO IRB=1,4 
         DO JXV=1,2 
          NJDX=0 
          ATHIV(JXV,IRB,IMAP,IV)=.FALSE. 
          DO JDX=1,25 
           IF(ATHBER(1,JDX,JXV,IRB,IV,IMAP).NE.-999.) NJDX=NJDX+1 
          END DO 
         IF(NJDX.EQ.0) CYCLE 
         IF(NJDX.EQ.25) THEN 
          ATHIV(JXV,IRB,IMAP,IV)=.TRUE. 
         ELSE 
          WRITE(IOUERR,3) IMAP,IV,IRB,JXV 
   3      FORMAT(' Error: All of the model year data for HAPS are not'/ 
     &           ' filled in. IMAP=',I2,' IV=',I2,' IRB=',I1,' JXV=',I1) 
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
