      REAL FUNCTION SOAKCF(MYSOK,CAT,IH,IWEEK,ST_INDX) 
C 
C  SOAKCF (previously SOAKWGT) calculates the a soak time averaged soak function. 
C 
C  Called by HRLOOP 
C 
C  Calls  
C 
C  Changes: (Last change first) 
C    
C  04 Dec 00 @EPA-bag Task X6, removed MY, JDX, and IDX from MYCODE.I; removed 
C            MYCODE.I and added model year, MYSOK, to the parameter list 
C  02 Jun 00 AIR Task 03: Changed code to only calculate current hour's 
C                         soak correction factor. 
C  26 apr 00 @EPA -bg    removed IV from parameter list; IV is only used in 
C                        diagnostic write statements which are commented out 
C  28 may 99 @DynTel -bg 2-650 Added code to write out results for a check with 
C            EPA documents. Output code is commented out. 
C  06 may 99 @DynTel -bg 2-650 Made both weekday and weekend SOAKDIST both 
C            fractions rather than weekday SOAKDIST percentages and weekend 
C            SOAKDIST fractions. So, there is no need for change of SOAKDIST 
C            numbers into fractions. Also, increased the size of the hourly 
C            index on SOAKDIST to 24. 
C   6 May 99 @DynTel-MLA 2-650 Added IHS because the old range of the IH 
C             parameter (1-14) has been expanded. The new range is 1-24. 
C  07 JAN 99 @DynTel -ZK 2-000 Explicit Typing                                   
C  20 Oct 98 @DynTel-BAG 2-000 Modified to decrease the MOBILE6 run time. 
C            No longer calls CALSOAK or SOAKADJ. CALSOAK now named SOAKFCN 
C            is called in the pollutant loop of HRLOOP. California soak 
C            function weighting by the fractions of vehicles with particular 
C            catalyst types is now done in SOAKCF by calling P3WPOX in 
C            HRLOOP (in the model year loop, IDX) and then passing the 
C            fraction, CAT, as a parameter. Model year, MY, is passed as a 
C            common block variable in MYCODE so that P3WPOX uses the same 
C            model year. 
C  29 Sep 98 @DynTel-BAG 2-000 Variable, ST, used in array was not declared  
C            as a integer.  
C 
C  Input on call or from calls: 
C 
C     parameter list: ICY, IWEEK 
C     common blocks: 
C     /SKDIST/ SOAKDIST 
C     /SOAK/   SOAKFUN 
C 
C  Output on return: 
C 
C    common block:  
C 
C  Local variable / array dictionary: 
C 
C   Name    Type              Description 
C  ------   ----  ------------------------------------------------ 
C  ST_INDX   I     soak times index (1,69) 
C  IV        I     vehicle type index 
C  IH        I     time of day index 
C  SOAKCF    R     Final soak time correction factors, i.e., the 
C                  soak time averaged soak function. The weighting in this 
C                  average is the percentage of vehicles, SOAKDIST, soaking 
C                  for a particular ST_INDX at a given hour, IH, of the day 
C                  (weekday or weekend, IWEEK) 
C 
C  Notes. 
C   
C  02 Apr 98 @DynTel-RJD 6-650 Created this subroutine for STARTS PROGRAM.      
C   4 Aug 98 @DynTel-MLA 6-650 Merged this subroutine into Mobile. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'SKDIST.I' 
      INCLUDE 'SOAK.I' 
C 
C     Declare parameter list 
C 
      INTEGER, INTENT(IN) :: MYSOK 
      INTEGER, INTENT(IN) :: IH 
      INTEGER, INTENT(IN) :: IWEEK 
      INTEGER, INTENT(IN) :: ST_INDX 
      REAL,    INTENT(IN) :: CAT 
C 
C     Declare local variables 
C 
      REAL    SOAKADJ 
C 
C 
C 
        IF (MYSOK.LE.1974) THEN 
          SOAKADJ = SOAKFUN(ST_INDX,1) 
        ELSE 
          SOAKADJ = (1.-CAT)*SOAKFUN(ST_INDX,1)+CAT* 
     *                       SOAKFUN(ST_INDX,2) 
        END IF 
        SOAKCF = SOAKDIST(IH,ST_INDX,IWEEK) * SOAKADJ 
C 
      RETURN 
      END 
