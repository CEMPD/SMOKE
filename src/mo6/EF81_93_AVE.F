      REAL FUNCTION EF81_93_AVE(MY,IP,IV,VMTAGE,IFDS) 
      USE USERDAT 
C 
C  EF81_93_AVE returns the average RUNNING emission factor. 
C  These emissions are determined from a fit of the average 
C  of high and normal emissions as a function of vehicle 
C  mileage/age. They are used when finding the fraction of 
C  high emissions vehicles for I&M credits and for weighting 
C  the high and normal emissions numbers(see HRLOOP.FOR). 
C 
C  Called by EF_HI_FRCN 
C 
C  Calls EF81_93,IERPTR 
C 
C  Changes: (Last change first) 
C 
C  22 Jun 00 EPA-BAG  removed JDX from EF81_93_AVE parameter list 
C  09 Jun 00 EPA-BAG  removed JDX from EF81_93 parameter list 
C  08 Sep 99 @DynTel-bag 2-660 Removed this code from BEREF.FOR and renamed 
C            to EF81_93_AVE.FOR. 
C 
C  Input on call or from calls: 
C 
C    argument list: MY,IP,IV,VMTAGE,IFDS 
C    common blocks: 
C    /LDGBER/ LDGRC,ITLDG 
C 
C  Output on return: 
C 
C    function: EF81_93_AVE 
C    common blocks: 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ----------------------------------------------------- 
C  KMILES   R    accumulated mileage in thousands of miles. 
C  ZPOINT   R    zero mileage point = intercept of ef curve 
C  SLOPE1   R    slope of ef curve (mileage < CRITP1) 
C  CRITP1   R    first critical point on ef curve 
C  SLOPE2   R    slope of ef curve (mileage > CRITP1) 
C  CRITP2   R    second critical point on ef curve 
C  SLOPE3   R    slope of ef curve (mileage > CRITP2) 
C 
C  Notes: 
C      
      IMPLICIT NONE 
      INCLUDE 'LDGBER.I' 
C 
C 
C     Declare external functions 
C 
      REAL    EF81_93 
      INTEGER IERPTR 
C 
C     Declare parameter list 
C 
      INTEGER, INTENT(IN) :: MY 
      INTEGER, INTENT(IN) :: IP 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: IFDS 
      REAL   , INTENT(IN) :: VMTAGE 
C 
C     Declare local variables 
C 
      REAL    KMILES 
      REAL    ZPOINT 
      REAL    CRITP1 
      REAL    CRITP2 
      REAL    SLOPE1 
      REAL    SLOPE2 
      REAL    SLOPE3 
      INTEGER LMY 
      INTEGER KEYER 
      INTEGER IVTAM 
      INTEGER ITG 
      INTEGER IGER 
      INTEGER, PARAMETER :: NORMAL=1,HIGH=2 
      INTEGER, PARAMETER :: RUN = 1 
C 
      TYPE (AVG_REC),POINTER :: USER_AVG 
      LOGICAL EXISTS 
C 
      LMY=MY 
      KMILES=10.0*VMTAGE 
C 
C  Emission factors can be calculated from data in BD04, BD40 or from  
C  user input data. The integer variable KEYER distinguishes between  
C  these cases. IGER is an integer index for the model year. 
C     
      IGER=IERPTR(RUN,LMY,IP,IV,KEYER) 
C 
C  Block data input BD40. The ITLDGV array stores a technolgy group  
C  pointer for each value of IGER and IFDS. 
C 
      IF(KEYER.EQ.0) THEN 
        IVTAM=1 
        IF (IV.NE.1) IVTAM=2 
C 
C Test if a user supplied record exists for ifds,ivtam,ip. User is a 
c returned pointer if such record exists. 
C 
        CALL CHK_AVG(IFDS,IVTAM,IP,MY,USER_AVG,EXISTS) 
C 
C If such a record exists use it instead of stored values 
C 
        IF (EXISTS) THEN 
          ZPOINT=USER_AVG%ZML 
          SLOPE1=USER_AVG%ADJ 
          SLOPE2=USER_AVG%DR1 
          SLOPE3=USER_AVG%DR2 
          CRITP1=USER_AVG%CPT1 
          CRITP2=USER_AVG%CPT2 
C 
        ELSE 
C 
          ITG=ITLDG(IGER,IFDS,IVTAM) 
          ZPOINT=LDGRC(1,ITG,IP,IVTAM) 
          SLOPE1=LDGRC(2,ITG,IP,IVTAM) 
          SLOPE2=LDGRC(3,ITG,IP,IVTAM) 
          SLOPE3=LDGRC(4,ITG,IP,IVTAM) 
          CRITP1=LDGRC(5,ITG,IP,IVTAM) 
          CRITP2=LDGRC(6,ITG,IP,IVTAM) 
C 
        ENDIF 
C 
        IF (KMILES.LE.CRITP1) THEN 
          EF81_93_AVE=ZPOINT+SLOPE1*KMILES 
        ELSEIF (KMILES.GE.CRITP1.AND.KMILES.LE.CRITP2) THEN 
          EF81_93_AVE=ZPOINT+SLOPE1*CRITP1+SLOPE2*(KMILES-CRITP1) 
        ELSE 
          EF81_93_AVE=ZPOINT+SLOPE1*CRITP1+SLOPE2*(CRITP2-CRITP1) 
     *                           +SLOPE3*(KMILES-CRITP2) 
        END IF 
C 
C  For very high accumulated mileages the average emissions may become 
C  higher than the high basic emission factor. In this case the average 
C  emissions will be set equal to the high emissions. 
C 
        IF(EF81_93_AVE .GE. 
     *     EF81_93(MY,IP,IV,VMTAGE,IFDS,RUN,HIGH)) THEN 
          EF81_93_AVE=EF81_93(MY,IP,IV,VMTAGE,IFDS,RUN,HIGH) 
        ENDIF 
C 
C  The average emissions may become lower than the normal basic emission 
C  factor. In this case the average emissions will be set equal to the 
C  normal emissions. 
C 
        IF(EF81_93_AVE .LE. 
     *     EF81_93(MY,IP,IV,VMTAGE,IFDS,RUN,NORMAL)) THEN 
          EF81_93_AVE=EF81_93(MY,IP,IV,VMTAGE,IFDS,RUN,NORMAL) 
        ENDIF 
C 
      END IF 
C 
      RETURN 
      END 
