C 
C  The File save common block for speed correction 
C 
C 
C 13 Feb 00 @EPA-elg  Modified FVMT and GLFVMT for additional IV array dimension. 
C  1 March 1999 @DynTel-HQ 2-698 speed correction related arrays 
C 
C  Local array subscripts: 
C 
C IH                     HOUR 
C FA                     Facility type 
C SP                     Speed distribution 
C COE                    coefficient 
C LE                     Level 
C 
C SVMT(4)                 STVMT(FA) 
C SHVMT(24)               SHVMT(IH) 
C SFVMT(4,24,28)          SFVMT(FA,IH,IV) 
C SSVMT(2,14,24)          SSVMT(FA,SP,IH) 
C SLOCAL(3,4)             SLOCAL(COE,IP) 
C SRAMP(3,4)              SRAMP(COE,IP) 
C SOCECO(2,4)             SSOCECO(LE,IP) 
C SOCY(3,4)               SOCY(COE,IP) 
C SOCX(3,4)               SOCX(COE,IP) 
C SLEVEL(3,4)             SLevel(LE,IP) 
C SEMSCF(4,4,3,2)         SEMSCF(IP, SP, LE, COE) 
C 
C 
      REAL STVMT, SHVMT, SFVMT, SSVMT, SLOCAL, SRAMP,  
     * SOCECO, SOCY, SOCX, SLEVEL, SEMSCF 
      COMMON /SAVE09/ STVMT(4), SHVMT(24), SFVMT(4, 24, 28), 
     * SSVMT(2, 14, 24), SLOCAL(3,4), SRAMP(3,4), SOCECO(2,4), 
     * SOCY(3,4), SOCX(3,4), SLEVEL(3,4), SEMSCF(2,3,4,4) 
