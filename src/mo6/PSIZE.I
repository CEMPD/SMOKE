C  PSIZE.I 
C 
C Changes: 
C 
C 12 Feb 01 AIR Task 11: New include for MOBILE6. 
C 12 Feb 01 AIR Task 11: Removed dust and particle size variables. 
C 04 Jan 01 AIR Task 08: Included MOIST in common PSIZE. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Common block dictionary: 
C 
C   Name    Type             Description 
C  ------   ----   ------------------------------------------------------- 
C  /PSIZE/: 
C 
C  PSBRK     R    fraction of PM for break wear 
C  PSDIE     R    fraction of PM for diesel fuel 
C  PSL       R    fraction of PM for leaded gasoline 
C  PSNL      R    fraction of PM for unleaded gasoline, no Catalyst 
C  PSNLCT    R    fraction of PM for unleaded gasoline, Catalyst 
C  PSTIRE    R    fraction of PM for tire wear 
C 
C 
C  The data are initialzed in BD02P5. 
C 
C 
      REAL PSBRK 
      REAL PSDIE 
      REAL PSL 
      REAL PSNL 
      REAL PSNLCT 
      REAL PSTIRE 
C 
      COMMON /PSIZE/ PSL, PSNL, PSNLCT, PSDIE, PSBRK, PSTIRE 
C 
