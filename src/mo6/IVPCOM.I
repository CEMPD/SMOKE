C 
C  Changes: (Last change first) 
C 
C 08 May 00 AIR Task 03: Added IVPTRR array for mapping BIGRVP calls. 
C  2 Feb 00 @DynTel-ddj, 1-032  Added IVPTRC array for mapping 
C           TIER2_50K array. 
C 15 Nov 96 DynTel@MLA Expanded the arrays to accomodate the new Heavy Duty 
C           Gasoline, Heavy Duty Diesel, and LDGT vehicle categories. 
C 
C  7 Apr 97 - DynTel@MLA  Expanded arrays IVPTRA, IVPTRB, IVPTRL, and IVPTRT  
C             to accomodate the new BUS vehicle types. 
C 
C  The data is initialized in BD16 
C 
      INTEGER IVPTRL 
      INTEGER IVPTRT 
      INTEGER IVPTRA 
      INTEGER IVPTRB 
      INTEGER IVPTRC 
      INTEGER IVPTRS 
      INTEGER IVPTRR 
C 
      COMMON /IVPCOM/ IVPTRL(28),IVPTRT(28),IVPTRA(28),IVPTRB(28), 
     &                IVPTRC(28),IVPTRS(28),IVPTRR(28) 
