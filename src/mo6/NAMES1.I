C  /NAMES1/: 
C 
C  Changes: (Last change first) 
C 
C  Common Block Dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C  NAMFLG   C*8   Control flag names for prompts and error messages 
C  NAMMMR   C*8   Prompts for replacement AMAR & JULMYR, by JDX & IV 
C  NAMRVP   C*8   Name of RVP falling outside allowable range 
C  NAMTEM   C*4   'mini' & 'maxi' prefixes for fuel tank temp error messages 
C 
C  Notes: 
C 
C  28 Aug 97 DynTel@MLA 2-622  NAMFLG was expanded from 26 to 27 to  
C            accomodate SSFLAG (the spreadsheet flag), which was added  
C            for Mobile6. 
C 
C  18 May 98 DynTel@MLA 2-622  Added RFGFLG to NAMFLG. 
C  20 May 98 DynTel@MLA 2-622  Added WINFLG to NAMFLG. 
C  24 Jun 98 DynTel@BAG 2-000  NAMFLG array reduced to 9 elements. 
C 
C The data is initialized in BD24 
C 
      CHARACTER*8 NAMFLG,NAMMMR,NAMRVP       
      CHARACTER*4 NAMTEM       
      COMMON /NAMES1/ NAMFLG(10),NAMMMR(2,3,8),NAMRVP(4),NAMTEM(2) 
