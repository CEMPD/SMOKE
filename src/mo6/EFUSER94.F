      MODULE EFUSER94 
 
      ! Values stored in the emissions record are the first and last 
      ! model years for the emission factor, the standards category, 
      ! the ZML, and the deterioration rate. 
 
C     Changes: 
C 
C     15 Sep 99 @DynTel-ddj 2-000  Added STANDARD index to USEREFS 
C               array. 
       
      TYPE :: EMISSREC 
        INTEGER                  ::  MYFIRST 
        INTEGER                  ::  MYLAST 
        INTEGER                  ::  STDCAT 
        REAL                     ::  ZML 
        REAL                     ::  DR1 
        !!! REAL                 ::  DR2 
        TYPE(EMISSREC), POINTER  ::  NEXT 
      END TYPE EMISSREC 
 
      ! The pointer in the PTRNODE Record points to a linked list 
      ! of EMISSREC records. 
 
      TYPE :: PTRNODE 
        TYPE(EMISSREC), POINTER ::  ROOT 
      END TYPE PTRNODE 
 
      ! The dimensions of the USEREF array are vehicle type, pollutant, 
      ! high/normal flag, start/running flag and standard. 
 
      TYPE (PTRNODE), DIMENSION(:,:,:,:,:), ALLOCATABLE  ::  USEREFS 
  
      END MODULE EFUSER94 
