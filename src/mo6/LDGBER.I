C Data is initialized in BD40 
C 
C  Changes: (Last change first) 
C 
C  2 Feb 00 @DynTel-ddj 1-032 Increased IV dimension on CRT050K to 8. 
C 17 Feb 99 @Dyntel-RJD 2-694 Moved ASEBZML's to LDGIMC.I. 
C 17 Feb 99 @Dyntel-RJD 2-694 Added CRT050K and CRT100K arrays. 
C 03 Sep 98  @DynTel-RJD 2-660 Removed redundant arrays. 
C 02 Jul 98 @DynTel-RJD 2-660 Contains coefficients for the calculation 
C            of separated start and running emissions. High emitters with 
C            faulty (but not tampered) emission control systems are also 
C            separated out. 
C 
      REAL NSEBZML,NSEBDR,HSEBZML,NREBZML,NREBDR,HREBZML, 
     *     LDGRC,CRT050K,CRT100K 
      INTEGER ITLDG 
      COMMON /LDGBER/ NSEBZML(13,3,5,3),NSEBDR(13,3,5,3), 
     *     HSEBZML(13,3,5,3),NREBZML(13,3,5,3), 
     *     NREBDR(13,3,5,3),HREBZML(13,3,5,3),LDGRC(6,7,3,2), 
     *     ITLDG(13,3,2),CRT050K(11,3,8),CRT100K(11,3,5) 
