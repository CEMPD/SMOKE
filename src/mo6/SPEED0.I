C 
C  Changes: (Last change first) 
C 
C  13 Jul 01 @EPA-djb Wish235, Add Average Speed command. 
C 
C  variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------ 
C  AVGFLG   I    Indicates when an average speed is to be used. 
C 
C                0 : Use speed distributions (default). 
C                1 : Use user supplied average speed value. 
C 
C  DEFAVG   R    Contains average speeds for default national average 
C                speed VMT arrays. 
C 
C  USRFAC   I    Index indicating which facility option chosen by 
C                the user determined from Average Speed command 
C                parameters. 
C 
C                1 : All VMT occurs on freeways, not including ramps. 
C                2 : All VMT occurs on freeways, including ramps. 
C                3 : All VMT occurs on arterial/collector roadways. 
C                4 : Average speed includes a mix of roadway types 
C                    specified by the user in USRVMT. 
C                5 : All VMT is assigned to local roadways. 
C                6 : All VMT is assigned to freeway ramps. 
C 
C  USRSPD   R    The average speed indicated in the Average Speed command. 
C 
C  USRVMT   R    The VMT distribution among the four roadway types 
C                specified by the user in the Average Speed command. 
C 
C                1 : Freeway VMT fraction (not including ramps). 
C                2 : Arterial/Collector VMT fraction. 
C                3 : Local roadway VMT fraction. 
C                4 : Freeway ramp VMT fraction. 
C 
C 
C These data are initialized in BDSPEED1.FOR 
C 
      REAL, DIMENSION(4)     :: DEFAVG 
      REAL                   :: USRSPD 
      REAL, DIMENSION(4)     :: USRVMT 
      INTEGER                :: USRFAC 
      INTEGER                :: AVGFLG 
C 
      COMMON /SPEED0/ USRSPD,USRVMT,USRFAC,DEFAVG,AVGFLG 
