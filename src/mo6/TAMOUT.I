C 
C  Changes: (Last change first) 
C 
C 26 Apr 00 @EPA-elg   Added arrays TEX and TAMBAS for Tampering effects. 
C 26 Nov 97 DynTel@MLA Expanded arrays TAMBAG, TCC, TDU, and THS to  
C           accomodate the new Heavy Duty Gasoline and LDGT vehicle  
C           categories. 
C 
C 10 Apr 97 - DynTel@MLA  Expanded arrays TAMBAG, TCC, TDU, and THS to  
C           accomodate the new GAS BUS vehicle type. 
C 
C 
C The data is initialized in BD25. 
C 
      REAL TAMBAG, THS, TDU, TCC, TAMBAS, TEX 
      COMMON /TAMOUT/ TAMBAG(3,3,25,14),THS(2,25,14), 
     *          TDU(25,14),TCC(25,14), 
     *          TAMBAS(3,3,25,14), TEX(28,3,2,25) 
