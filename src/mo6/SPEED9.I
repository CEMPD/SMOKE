C  The File include the variable for speed correction 
C 
C  Local array subscripts: 
C 
C 13 Feb 00 @EPA-elg  Modified FVMT and GLFVMT for additional IV array dimension. 
C 28 Nov 00 @EPA-bag  bug301 added GLHVMT, GLFVMT, GLSVMT for run AND scenario 
C           level user input 
C 30 May 00 @EPA-elg  Eliminated BDSPEED1.i and moved arrays SFCAC AND SFCF 
C                     to Speed9.i 
C 27 Jan 00 @DynTel-MLQA 1-011 Changed IV dimension of BEFIDL to 28 
C           from 27 for LDDT34. 
C 12 Jan 00 @DynTel-ddj 1-013  Added ramp and local average speed 
C           variables AVSPR and AVSPL. 
C  1 March 1999 @DynTel-HQ 2-698 speed correction related arrays 
C 
C These data are initialized in BDSPEED2.FOR 
C 
C IH                     HOUR 
C FA                     Facility type 
C SP                     Speed distribution 
C COE                    coefficients 
C LE                     Level 
C 
C VMT(4)                 TVMT(FA) 
C HVMT(24)               HVMT(IH) 
C GHVMT(24)              GHVMT(IH) 
C FVMT(4,24,28)          FVMT(FA,IH,IV) 
C SVMT(2,14,24)          SVMT(FA,SP,IH) 
C LOCAL(3,4)             LOCAL(COE,IP) 
C RAMP(3,4)              RAMP(COE,IP) 
C OCECO(2,4)             OCECO(COE,IP) 
C OCY(3,4)               OCY(COE,IP) 
C OCX(3,4)               OCX(COE,IP) 
C OCEF(4,4)              OCEF(COE,IP) 
C OCLEVEL(3,4)           OCLEVEL(COE,IP) 
C LEVEL(3,4)             Level(LE,IP) 
C SCFF7(3,4)             SCFF7(LE,IP) 
C SCFAC7(3,4)            SCFAC7(LE,IP) 
C SCFF(3,13,4)           SCFF(LE,SP,IP) 
C SCFAC(3,13,4)          SCFAC(LE,SP,IP) 
C EMSCF(4,4,3,2)         EMSCF(IP, SP, LE, COE) 
C M5SCF(14,3)            M5SCF(SP,IP) 
C BEFIDL(3,25,28,24)     BEFIDL(IP,IDX,IV,IH) 
C 
C SFCAC(3,13,4) 
C SFCF(3,13,4) 
C 
C 
      REAL TVMT, HVMT, FVMT, SPEED, SVMT, LOCAL, RAMP, AVESPD, 
     * OCECO, OCY, OCX, OCEF, OCLEVEL, SCFAC, LEVEL, SCFF, EMSCF, 
     * SCFF7,SCFAC7, M5SCF,SCFL,SCFR,BEFIDL, AVSPR, AVSPL, 
     * SFCAC, SFCF 
C 
      REAL GLHVMT, GLFVMT, GLSVMT 
C 
      COMMON  /SPEED9/ TVMT(4), HVMT(24),FVMT(4, 24, 28), 
     * SPEED(15), SVMT(2, 14, 24), AVESPD(14,2), 
     * LOCAL(3,4), RAMP(3,4), 
     * OCECO(2,4), OCY(3,4), OCX(3,4), OCEF(4,4), 
     * OCLEVEL(3,4),SCFF7(3,4),SCFAC7(3,4),SCFAC(3,14,4), 
     * LEVEL(3,4), SCFF(3,14,4), M5SCF(14,3),EMSCF(2,3,4,4), 
     * SCFL, SCFR, BEFIDL(3,25,28,24), AVSPR, AVSPL, 
     * SFCAC(3,13,4), SFCF(3,13,4), GLHVMT(24), GLFVMT(4,24,28), 
     * GLSVMT(2,14,24) 
