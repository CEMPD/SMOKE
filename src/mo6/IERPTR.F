      INTEGER FUNCTION IERPTR(ISR,MY,IP,IV,KEYER) 
C 
C  IERPTR selects the array pair and the model year group pointers  
C  into the intercept and slope arrays for the emission rate  
C  equations. 
C 
C  Called by EF_MCDSLPRE81, BEFHDV 
C 
C  Calls QUITER. 
C 
C  Changes: (Last change first) 
C  22 Jun 00 @EPA-BAG    Used INTENT qualifier in parameter declarations 
C  08 Sep 99 @DynTel-BAG 2-660 Rearranged variable declarations/comments 
C   6 JAN 99 @DynTel-HXQ 2-668 Pre-1981 LDGV & LDGT basic exhaust 
C             engine start and running 
C  17 Nov 98 @DynTel-ZK 2-000 explicit Typing 
C  14 Sep 98 @DynTel-RJD 2-660  
C  27 Aug 98 @DynTel-RJD 0-632 Removed fuel delivery system parameter 
C            IFDS from the argument list. 
C  12 Aug 98 @DynTel-RJD 2-663 Added fuel delivery system parameter  
C            IFDS to the argument list. 
C             
C 
C  Input on call: 
C 
C    argument list: ISR,MY,IP,IV 
C    common blocks: 
C    /BASEQ2/ MYGERB,MAXERB 
C    /BASEQ6/ MYGERU,NUMERU 
C    /FLAGS2/ NEWFLG 
C    /VVNAM/  VVLDGAS (IVTYPE) 
C    /REGION/ IREJN 
C 
C  Output on return: 
C    argument list: KEYER 
C    function: IERPTR 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  NERU     I    number of user supplied ER equations for the IPxIVxIREJN 
C                case passed in on the call to IERPTR 
C 
C  Notes: 
C 
C  STOREQ looks for the match in first the user supplied (if any) and then the 
C  hardcoded (initialized in BLOCK DATA) equation myg arrays. 
C 
C  The user is responsible for insuring entered ER records do not have 
C  overlapping mygs.  IERPTR returns the index to the first myg match. 
C 
C  12 May 98 @DynTel-RJD 2-652 Added the ISR variable to the input parameter 
C            list. Expanded NUMERU and MYGERU arrays to accommodate separate 
C            start and running emissions. 
C  14 Jul 98 @DynTel-RJD 2-660 Partitioned basic emission rate coefficients 
C            into 3 catagories. Hardcoded coefficients for 81-93 LDGV  
C            vehicles (KEYER=0), hardcoded coefficients for all other 
C            vehicles (KEYER=1) and user input coefficients (KEYER=2). The 
C            KEYER=0 case is new. 
C 
      IMPLICIT NONE 
      INCLUDE 'BASEQ2.I' 
      INCLUDE 'BASEQ6.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'REGION.I' 
C 
C     Declare external functions 
C 
C 
C     Declare parameter list 
C      
      INTEGER, INTENT(IN) :: ISR 
      INTEGER, INTENT(IN) :: MY 
      INTEGER, INTENT(IN) :: IP 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(OUT) :: KEYER 
C 
C     Declare local variables 
C 
      INTEGER NERU,IG 
      KEYER = 0 
C      
      NERU=NUMERU(IP,IV,IREJN,ISR) 
C 
C     use mobile5 default BER's 
C 
      IF(NEWFLG.EQ.1.OR.NERU.EQ.0) GOTO 20 
C 
C  Might be a match in the user supplied myg's. 
C 
      KEYER=2 
      DO 10 IG=1,NERU 
         IF(MY.GE.MYGERU(IG,1,IP,IV,IREJN,ISR).AND. 
     *      MY.LE.MYGERU(IG,2,IP,IV,IREJN,ISR)) GOTO 90 
   10 CONTINUE 
C 
C  Look in the Mobile5 Block Data hardcoded myg's. 
C 
   20 IF (VVLDGAS(IV).EQ.1.AND.MY.GE.1981.AND.MY.LE.1993) THEN 
        KEYER=0 
        IG=MY-1980 
        GOTO 90 
      ELSE 
        KEYER=1 
        DO 30 IG=1,MAXERB 
          IF(MY.LE.MYGERB(IG,IP,IV,IREJN)) GOTO 90 
   30   CONTINUE 
      ENDIF 
C 
C  Set default. 
C 
      IG=1 
      CALL QUITER(0.,11,61,MY) 
C 
   90 IERPTR=IG 
C 
      RETURN 
      END 
