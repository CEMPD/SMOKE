C  Changes : (Last change first) 
C 
C  10 Jul 01 @EPA-djb Wish130, Expanded STARTSPERDAY array to include age. 
C  21 Nov 00 @EPA-BAG BUG#300 Add VA021 for SOAKDIST 
C  23 May 00 @EPA-djb Add VA019 & VA020 for STARTSPERDAY and STARTDIST. 
C  29 Sep 98 @DynTel-HQ 2-678  Add CONVFLG for conversion factor  
C 
C  The data is initialized in BDSAVE.FOR 
C  These are saved conversion factors 
C 
      REAL VA018, VA019, VA020, VA021 
      COMMON /SAVE08/VA018(56,19),VA019(25,28,2),VA020(24,2), 
     &               VA021(24,70,2) 
