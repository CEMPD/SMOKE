C 
C  18 Jan 01 @EPA-bag 2007 rule for heavy duty vehicles; increased 
C         size of MYGERB 
C  22 Oct. 98 - DynTel@HXQ  Expanded arrays MYGERB for new model years      
C 
C   5 Nov 96 DynTel@MLA  Expanded array MYGERB to accomodate the new  
C                        Heavy Duty Gasoline, Heavy Duty Diesel, and 
C                        LDGT vehicle types. 
C 
C  2 Apr 97 - DynTel@MLA  Expanded array MYGERB to accomodate the new BUS  
C             vehicle categories. 
C 
C  The data is initialized in BD03. 
C 
      INTEGER MYGERB,MAXERB 
      COMMON /BASEQ2/ MYGERB(35,3,28,2),MAXERB 
