C 
C  Changes: (Last change first) 
C 
C 19 May 99 DynTel@HXQ  Delete arrays BMYTPD 
C  3 Dec 96 DynTel@MLA Expanded arrays BMYMPD, BMYTPD, and BMTPDF to  
C           accomodate the new Heavy Duty Gasoline, Heavy Duty Diesel,  
C           and LDGT vehicle categories. 
C 
C 14 Apr 97 - DynTel@MLA  Expanded arrays BMTPDF, BMYMPD, and BMYTPD to  
C           accomodate the new BUS vehicle types. 
C 
C The data is initialized in BD27. 
C  
      REAL BMYMPD,BMTPDF 
      COMMON /CEVBMY/ BMYMPD(25,28),BMTPDF(25,28) 
