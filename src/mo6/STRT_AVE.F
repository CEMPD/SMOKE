      REAL FUNCTION STRT_AVE(MY,IP,IV,VMTAGE,IFDS) 
C 
C  STRT_AVE returns the average START emission factor for 
C  1981-93 Model Year Light duty gas vehicles.  These emission 
C  factors are used in the calculation of start methane emissions. 
C 
C  Called by BIGMET 
C 
C  Calls EF81_93,IERPTR 
C 
C  Changes: (Last change first) 
C 
C  22 Jun 00 EPA-BAG  removed JDX from EF81_93_AVE parameter list 
C  09 Jun 00 EPA-BAG  removed JDX from EF81_93 parameter list 
C  24 Apr 00 @EPA-ELG  New Subroutine 
C 
C  Input on call or from calls: 
C 
C    argument list: MY,IP,IV,VMTAGE,IFDS 
C    common blocks: 
C    /LDGBER/ LDGRC,ITLDG 
C 
C  Output on return: 
C 
C    function: STRT_AVE 
C    common blocks: 
C 
C  Local variable / array dictionary: 
C 
C   Name      Type              Description 
C  ------     ----  ----------------------------------------------------- 
C  NUM_HIGH    R    Fraction of High Emitters in both start and running for 1981+ 
C  RUN_AVE     R    Average Running emissions 
C  RUN_HIGH    R    Average High emitter emission level for running 
C  RUN_NORM    R    Average Normal emitter emission level for running 
C  STRT_AVE    R    Average Start emissions 
C  STRT_HIGH   R    Average High emitter emission level for starts 
C  STRT_NORM   R    Average Normal emitter emission level for starts 
C 
C  Notes: 
C      
      IMPLICIT NONE 
C 
C 
C     Declare external functions 
C 
      REAL    EF81_93 
      REAL    EF81_93_AVE 
C 
C     Declare parameter list 
C 
      INTEGER, INTENT(IN) :: IFDS 
      INTEGER, INTENT(IN) :: IP 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: MY 
C 
      REAL NUM_HIGH 
      REAL RUN_AVE 
      REAL RUN_HIGH 
      REAL RUN_NORM 
      REAL STRT_HIGH 
      REAL STRT_NORM 
      REAL VMTAGE 
C 
C     Declare local variables 
C 
      INTEGER, PARAMETER :: NORMAL=1, HIGH=2 
      INTEGER, PARAMETER :: RUN = 1, STRT=2 
C 
C 
      RUN_AVE   = EF81_93_AVE(MY,IP,IV,VMTAGE,IFDS) 
      RUN_HIGH  = EF81_93(MY,IP,IV,VMTAGE,IFDS,RUN,HIGH) 
      RUN_NORM  = EF81_93(MY,IP,IV,VMTAGE,IFDS,RUN,NORMAL) 
      NUM_HIGH  = (RUN_AVE-RUN_NORM) / (RUN_HIGH-RUN_NORM) 
C 
       NUM_HIGH  = MAX(0.0,NUM_HIGH) 
       NUM_HIGH  = MIN(1.0,NUM_HIGH) 
C 
      STRT_HIGH = EF81_93(MY,IP,IV,VMTAGE,IFDS,STRT,HIGH) 
      STRT_NORM = EF81_93(MY,IP,IV,VMTAGE,IFDS,STRT,NORMAL) 
C 
      STRT_AVE  = STRT_NORM*(1.0-NUM_HIGH) + STRT_HIGH*NUM_HIGH 
C 
      RETURN 
      END 
