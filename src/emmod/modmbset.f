
        MODULE MODMBSET
             
        INTEGER, ALLOCATABLE, PUBLIC :: MCREFSORT ( :,: ) ! sorted MCREF data  
        INTEGER, ALLOCATABLE, PUBLIC :: MCREFIDX ( :,: )  ! index into MCREF by ref. county
        
        INTEGER, ALLOCATABLE, PUBLIC :: MVREFSORT ( :,: ) ! sorted MVREF data
        
        INTEGER, ALLOCATABLE, PUBLIC :: SCENLIST ( :,: )  ! scenario number and local-as-arterial
                                                          ! flag for each source
        
        INTEGER, PUBLIC :: NREFC           ! no. of reference counties
        INTEGER, PUBLIC :: NINVC           ! no. of unique counties in inventory
        
        INTEGER, PUBLIC :: NREFFLAGS = 3   ! no. of settings flags in MVREF file
        
        CHARACTER*300, ALLOCATABLE :: M6LIST( : )  ! contents of M6LIST file
        
        INTEGER, ALLOCATABLE, PUBLIC :: COUNTYSRC ( : )   ! county FIP code for each source - read
                                                          ! from SPDSUM file (handles spatial averaging)
        
        REAL, ALLOCATABLE, PUBLIC :: EMISSIONS( : ) ! array to hold M6 results
        INTEGER, PUBLIC :: EMISPOS = 0     ! current position in EMISSIONS array
        
C.........  Various constants for file options, road types, etc.
        INTEGER, PARAMETER, PUBLIC :: RURALINTERSTATE = 1   ! rural interstate
        INTEGER, PARAMETER, PUBLIC :: RURALPRINCART   = 2   ! rural principle arterial
        INTEGER, PARAMETER, PUBLIC :: RURALMINORART   = 6   ! rural minor arterial
        INTEGER, PARAMETER, PUBLIC :: RURALMAJORCOLL  = 7   ! rural major collector
        INTEGER, PARAMETER, PUBLIC :: RURALMINORCOLL  = 8   ! rural minor collector
        INTEGER, PARAMETER, PUBLIC :: RURALLOCAL      = 9   ! rural local
        INTEGER, PARAMETER, PUBLIC :: URBANINTERSTATE = 11  ! urban interstate
        INTEGER, PARAMETER, PUBLIC :: URBANFREEWAY    = 12  ! urban freeway
        INTEGER, PARAMETER, PUBLIC :: URBANPRINCART   = 14  ! urban principle arterial
        INTEGER, PARAMETER, PUBLIC :: URBANMINORART   = 16  ! urban minor arterial
        INTEGER, PARAMETER, PUBLIC :: URBANCOLL       = 17  ! urban collector
        INTEGER, PARAMETER, PUBLIC :: URBANLOCAL      = 19  ! urban local
       
        INTEGER, PARAMETER, PUBLIC :: FREEWAY  = 1    ! MOBILE6 freeway sources
        INTEGER, PARAMETER, PUBLIC :: ARTERIAL = 2    ! MOBILE6 arterial sources
        INTEGER, PARAMETER, PUBLIC :: LOCAL    = 3    ! MOBILE6 local sources
        INTEGER, PARAMETER, PUBLIC :: RAMP     = 4    ! MOBILE6 ramp sources
        INTEGER, PARAMETER, PUBLIC :: NONE     = 5    ! MOBILE6 non-facility
        
        INTEGER, PARAMETER, PUBLIC :: DAILY    = 1    ! daily temperature profiles
        INTEGER, PARAMETER, PUBLIC :: WEEKLY   = 2    ! weekly temperature averaging
        INTEGER, PARAMETER, PUBLIC :: MONTHLY  = 3    ! monthly temperature averaging
        INTEGER, PARAMETER, PUBLIC :: EPISLEN  = 4    ! episode length temp. averaging
        
        END MODULE MODMBSET
        