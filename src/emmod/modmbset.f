
        MODULE MODMBSET
     
        INTEGER, ALLOCATABLE, PUBLIC :: MCREFSORT ( :,: ) ! sorted MCREF data  
        INTEGER, ALLOCATABLE, PUBLIC :: MCREFIDX ( :,: )  ! index into MCREF by ref. county
        
        INTEGER, ALLOCATABLE, PUBLIC :: MVREFSORT ( :,: ) ! sorted MVREF data
        
        INTEGER, ALLOCATABLE, PUBLIC :: SCENLIST ( : )    ! scenario number for each source
        
        INTEGER, PUBLIC :: NREFC           ! no. of referenc counties
        INTEGER, PUBLIC :: NINVC           ! no. of unique counties in inventory
        
        INTEGER, PUBLIC :: NREFFLAGS = 3   ! no. of settings flags in MVREF file
        
        CHARACTER*300, ALLOCATABLE :: M6LIST( : )  ! contents of M6LIST file
        
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
        
        INTEGER, PARAMETER, PUBLIC :: DAILY    = 1    ! daily temperature profiles
        INTEGER, PARAMETER, PUBLIC :: WEEKLY   = 2    ! weekly temperature averaging
        INTEGER, PARAMETER, PUBLIC :: MONTHLY  = 3    ! monthly temperature averaging
        INTEGER, PARAMETER, PUBLIC :: METLEN   = 4    ! met file length temp. averaging
        
        END MODULE MODMBSET
        