      MODULE ATHEAP 
C 
C  This module sets up the large DBAT, ATHBER, and ADDITIONAL HAPS arrays 
C  to allocatable status. 
C 
C  Changes: (Last change first) 
C 
C 08 Mar 02: AIR Task 28: Moved ADDITIONAL HAPS related arrays to ATWORK. Updated Comments 
C 28 Nov 01: AIR Task 28: Moved ATHBER to ATWORK. 
C 04 Sep 01: AIR Task 28: New MODULE. 
C 
C  ATHBER    R     Air Toxic HAPS basic ratios/ZMLs and BERs 
C  DBAT      R     The Air Toxic exhaust and evap emission factors, by hour, road type, vehicle class 
C  EFAT      R     The Air Toxic exhaust and evap emission factors, by vehicle class 
C  VAT       R     The Air Tosic exhaust and evap "All Vehicle" emission factors 
C 
C  ATHBER(IRZD,JDX,JXV,IRB,IV,IMAP) where: 
C     IRZD = 1=Ratio/ZML (Ratio is in mg/mi per g/mi, ZML is in mg/mi), 2=Deterioration Rate, mg/mi per 10K miles 
C     JDX  = Age, 1 to 25 years 
C     JXV  = 1=eXhaust, 2=eVap  
C     IRB  = Ratio/BER index 1=VOC, 2=TOG, 3=BEF, 4=PM 
C     IV   = Vehicle class, 1-28 
C     IMAP = Pollutant ID Number, mapped into 1-50 placeholders 
C 
C  DBAT(IDX,IH,IROAD,IXV,IAT,IV) where: 
C     IDX   = Age, 1 to 25 years 
C     IH    = Hour, 1-24 
C     IROAD = Road Type, 1-5 
C     IAT   = Air Toxic type: 1=Benz, 2=MTBE, 3=Buta, 4=Form, 5=Acet, 6=Acrolein 
C     IV    = Vehicle Class, 1-28 
C 
C  EFAT(IXV,IAT+NMAP,IV) where: 
C     IXV  = 1=eXhaust, 2=eVap  
C     IAT+NMAP = number of Air Toxics plus ADD HAPS pollutants 
C     IV   = Vehicle class, 1-28 
c 
C  VAT(IXV,IAT+NMAP) 
C     IXV  = 1=eXhaust, 2=eVap  
C     IAT+NMAP = number of Air Toxics plus ADD HAPS pollutants 
C 
C  ATHIV(JXV,IRB,IMAP,IV) 
C     JXV  = 1=eXhaust, 2=eVap  
C     IRB  = Ratio/BER index 1=VOC, 2=TOG, 3=BEF, 4=PM 
C     IMAP = Pollutant ID Number, mapped into 1-50 placeholders 
C     IV   = Vehicle class, 1-28 
C 
C  The following single-dimensioned arrays are associated with ADDITIONAL HAPS: 
C 
C  ATHBEGMY  I     Air Toxic HAPS beginning model year 
C  ATHDET    R     Air Toxic HAPS deterioration rates, mg/10k miles 
C  ATHEKIND  C     Air Toxic HAPS emission kind 
C  ATHENDMY  I     Air Toxic HAPS ending model year 
C  ATHPMAP   I     Air Toxic HAPS pollutant name mapping index 
C  ATHPNAME  C     Air Toxic HAPS pollutant name 
C  ATHPOL    I     Air Toxic HAPS pollutant index 
C  ATHRATIO  C     Air Toxic HAPS ratio to MOBILE6 VOC/TOG emission factors 
C  ATHRATZML R     Air Toxic HAPS ratio or zero mile level 
C  ATHSNAME  C     Air Toxic HAPS shortened pollutant name, makes room for Exhaust/Evap type extensions 
C  ATHVTYPE  I     Air Tosic HAPS vehicle type 
C  ATHIV     L     Air Toxic HAPS logical data available switch 
C 
C 
      REAL,         DIMENSION(:,:)        , ALLOCATABLE :: VAT 
      REAL,         DIMENSION(:,:,:)      , ALLOCATABLE :: EFAT 
      REAL,         DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: ATHBER 
      REAL,         DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: DBAT 
C 
      CHARACTER*4,  DIMENSION(:),           ALLOCATABLE :: ATHEKIND 
      CHARACTER*8,  DIMENSION(:),           ALLOCATABLE :: ATHRATIO 
      CHARACTER*12, DIMENSION(:),           ALLOCATABLE :: ATHSNAME 
      CHARACTER*32, DIMENSION(:),           ALLOCATABLE :: ATHPNAME 
      INTEGER,      DIMENSION(:),           ALLOCATABLE :: ATHBEGMY 
      INTEGER,      DIMENSION(:),           ALLOCATABLE :: ATHENDMY 
      INTEGER,      DIMENSION(:),           ALLOCATABLE :: ATHPMAP 
      INTEGER,      DIMENSION(:),           ALLOCATABLE :: ATHPOL 
      INTEGER,      DIMENSION(:),           ALLOCATABLE :: ATHVTYPE   
      LOGICAL,      DIMENSION(:,:,:,:),     ALLOCATABLE :: ATHIV 
      REAL,         DIMENSION(:),           ALLOCATABLE :: ATHDET 
      REAL,         DIMENSION(:),           ALLOCATABLE :: ATHRATZML 
C 
      END MODULE 
C 
C---------------------------------------------------------------------- 
C 
      SUBROUTINE ATALLOC 
C 
C  This routine allocates the large DBAT, ATHBER, and ADDITIONAL HAPS arrays 
C 
C  Changes: (Last change first) 
C 
C 08 Mar 02: AIR Task 28: Moved ADDITIONAL HAPS related arrays to ATWORK. 
C 28 Nov 01: AIR Task 28: Moved ATHBER to ATWORK. 
C 04 Sep 01: AIR Task 28: New Subroutine. 
C 
      USE ATHEAP 
C 
      ALLOCATE (ATHBER(2,25,2,4,28,50)) 
      ALLOCATE (DBAT  (25,24,5,6,6,28)) 
      ALLOCATE (EFAT  (6,56,28)) 
      ALLOCATE (VAT   (6,56)) 
      ALLOCATE (ATHIV (2,4,50,28)) 
C 
      ALLOCATE (ATHBEGMY (1400)) 
      ALLOCATE (ATHDET   (1400)) 
      ALLOCATE (ATHEKIND (1400)) 
      ALLOCATE (ATHENDMY (1400)) 
      ALLOCATE (ATHPMAP  (1400)) 
      ALLOCATE (ATHPNAME (1400)) 
      ALLOCATE (ATHPOL   (1400)) 
      ALLOCATE (ATHRATIO (1400)) 
      ALLOCATE (ATHRATZML(1400)) 
      ALLOCATE (ATHSNAME (1400)) 
      ALLOCATE (ATHVTYPE (1400)) 
C 
      RETURN 
      END 
C 
C---------------------------------------------------------------------- 
C 
      SUBROUTINE ATDEALLOC 
C 
C  This routine deallocates the large DBAT, ATHBER, and ADDITIONAL HAPS arrays 
C 
C  Changes: (Last change first) 
C 
C 08 Mar 02: AIR Task 28: Moved ADDITIONAL HAPS related arrays to ATWORK. 
C 28 Nov 01: AIR Task 28: Moved ATHBER to ATWORK. 
C 04 Sep 01: AIR Task 28: New Subroutine. 
C 
      USE ATHEAP 
C 
      DEALLOCATE (ATHBER) 
      DEALLOCATE (DBAT) 
      DEALLOCATE (EFAT) 
      DEALLOCATE (VAT) 
C 
      DEALLOCATE (ATHBEGMY) 
      DEALLOCATE (ATHDET) 
      DEALLOCATE (ATHEKIND) 
      DEALLOCATE (ATHENDMY) 
      DEALLOCATE (ATHIV) 
      DEALLOCATE (ATHPMAP) 
      DEALLOCATE (ATHPNAME) 
      DEALLOCATE (ATHPOL) 
      DEALLOCATE (ATHRATIO) 
      DEALLOCATE (ATHRATZML) 
      DEALLOCATE (ATHSNAME) 
      DEALLOCATE (ATHVTYPE) 
C 
      RETURN 
      END 
