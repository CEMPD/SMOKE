C 
C  Changes: (Last change first) 
C 
C 25 Jan 2000 @EPA-djb added arrays, AGED(25,30),DAF(30). 
C 03 Mar 98 - DynTel@RJD 2-204  Changed the name of the array storing the Model 
C           Year Registration fractions from JANMYR to MEVMYR to reflect the 
C           fact that January is no longer the only month of evaluation that 
C           the Mobile model can handle. 
C 
C 14 Apr 97 - DynTel@MLA  Expanded arrays MEVMYR, TF, TFMYM, and XMYM to  
C           accomodate the new BUS vehicle types. 
C 
C  3 Dec 96 DynTel@MLA Expanded arrays MEVMYR, TF, TFMYM, and XMYM to  
C           accomodate the new Heavy Duty Gasoline, Heavy Duty Diesel,  
C           and LDGT vehicle categories. 
C 
C 
C The data is initialized in BD27. 
C 
      REAL XMYM,MEVMYR,TF,TFMYM, AGED, DAF 
      COMMON /MYRCAL/ XMYM(25,28),MEVMYR(25,28),TF(25,28),TFMYM(25,28), 
     *                AGED(25,30),DAF(30) 
