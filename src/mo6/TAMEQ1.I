C 
C  Changes: (Last change first) 
C 
C 26 Nov 97 DynTel@MLA Expanded arrays TAMZML, TAMDR, and TAMA50 to  
C           accomodate the new Heavy Duty Gasoline and LDGT vehicle  
C           categories. 
C 
C 10 Apr 97 - DynTel@MLA  Expanded arrays TAMA50, TAMDR, and TAMZML to  
C           accomodate the new GAS BUS vehicle type. 
C 
C 
C The data is initialized in BD25. 
C 
      REAL TAMZML, TAMDR, TAMA50 
      COMMON /TAMEQ1/ TAMZML(9,14,3,2),TAMDR(9,14,3,2),TAMA50(9,14,3,2)                                 
