      SUBROUTINE PMPARSIZ 
C 
C     This subroutine computes the fraction of particles less than 
C     or equal to the particle size cutoff input by the user, PSZCUT. 
C 
C Changes: 
C 
C 12 Feb 01 AIR Task 11: Removed dust variables. 
C 04 Jan 01 AIR Task 08: Corrected documentation. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C     Called by:  PMCALX 
C 
C     Calls: none 
C 
C     Input on call: 
C     /PSIZE/ PSZCUT 
C     /PSBLK/ PSDIST(5,2,8) 
C 
C     Output on return: 
C     /PSIZE/ PSL, PSNL, PSNLCT, PSDIE, PSBRK, PSTIRE 
C     /PSDUST/ PSDUNP,PSDPVD 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  MXDATA   I     The maximum number of size bins for the array PSDIST 
C  MXSRCE   I     The maximum number of categories of sources of PM in 
C                    the array PSDIST. 
C  PSCALC   R     The minimum PS cut off for fugitive dust is now 2.5 
C 
C 
C  Notes: 
C  There are eight categories of sources: 
C 
C      1   PSL         fraction for leaded gasoline 
C      2   PSNLCT      fraction for unleaded gasoline, Catalyst 
C      3   PSNL        fraction for unleaded gasoline, no Catalyst 
C      4   PSDIE       fraction for diesel fuel 
C      5   PSBRK       fraction for break wear 
C      6   PSTIRE      fraction for tire wear 
C 
C     The fractions are linearly interpolated from the data in psdist, 
C 
C          PSDIST=PSDIST(ipts,I, category) 
C     where 
C          ipts = the number of particle sizes for which there is 
C                 data about the fraction of total particles less 
C                 than or equal to those sizes (MXDATA=5) 
C          I = 1  corresponds to the particle size 
C          I = 2  for sources 1-6: corresponds to the fraction of particles 
C                 less than or equal to that size 
C                 for sources 7-8: base emission factor (k-factor) at the 
C                 given particle size 
C          category = the source category, i.e. the same eight listed 
C                 above in the same order. 
C 
C 
C  Jul-21-1994 @ CSC-bsg Changes output of fugitive dust to reflect 
C  the new lower limit of Particle Size Cutoff if less than 2.5 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'LEAD.I' 
      INCLUDE 'PART1.I' 
      INCLUDE 'PSIZE.I' 
C 
C  Declare local variables/arrays. 
C 
      INTEGER I 
      INTEGER IDATA 
      INTEGER ISRC 
      INTEGER MXSIZE 
      INTEGER NPTR 
      REAL    PS 
      REAL    PSCALC 
      REAL    SLOPE 
C 
      INTEGER,PARAMETER :: MXDATA = 5 
      INTEGER,PARAMETER :: MXSRCE = 6 
C 
C 
      DO 50 ISRC = 1, MXSRCE 
       PSCALC=PSZCUT 
        PS = 0.0 
C 
       IF (PSZCUT .LE. PSDIST(1,1,ISRC)) THEN 
          SLOPE = (PSDIST(2,2,ISRC) - PSDIST(1,2,ISRC)) / 
     *              (PSDIST(2,1,ISRC) - PSDIST(1,1,ISRC)) 
          PS = PSDIST(1,2,ISRC) - SLOPE * (PSDIST(1,1,ISRC) - PSCALC) 
       ELSE 
          DO 10 IDATA = 1, MXDATA 
             SLOPE = 0.0 
             NPTR = 1 
             MXSIZE = MXDATA 
             DO 5 I = 2, MXDATA 
              IF (PSDIST(I,1,ISRC).EQ.0.0 .AND. 
     *                       PSDIST(I-1,1,ISRC).GT.0.) MXSIZE = I-1 
5              CONTINUE 
             DO 20 I = 2, MXSIZE 
              IF (PSDIST(I,1,ISRC) .GE. PSCALC) THEN 
                 NPTR = I 
                 GO TO 30 
              ENDIF 
   20          CONTINUE 
   30          SLOPE = (PSDIST(NPTR,2,ISRC) - PSDIST(NPTR-1,2,ISRC)) / 
     *                 (PSDIST(NPTR,1,ISRC) - PSDIST(NPTR-1,1,ISRC)) 
             PS = SLOPE*(PSCALC - PSDIST(NPTR-1,1,ISRC)) + 
     *                              PSDIST(NPTR-1,2,ISRC) 
   10       CONTINUE 
       ENDIF 
C 
       IF (ISRC .EQ. 1) PSL = PS 
       IF (ISRC .EQ. 2) PSNLCT = PS 
       IF (ISRC .EQ. 3) PSNL = PS 
       IF (ISRC .EQ. 4) PSDIE = PS 
       IF (ISRC .EQ. 5) PSBRK = PS 
       IF (ISRC .EQ. 6) PSTIRE = PS 
C 
   50 CONTINUE 
C 
      RETURN 
      END 
