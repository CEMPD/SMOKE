C 
C The data is initialized in BD27. 
C 
C  Changes: (Last change first) 
C 
C  16 Apr 01 @EPA-elg Added POST96IM Flag to control diagnostic I/M command. 
C  20 Nov 00 @EPA-elg Added CAAFLG for No clean air act. 
C  18 Jul 00 @EPA-bag bug#54 removed TEMFLG 
C  19 Jun 98 @DynTel-bag 2-000  Removed OUTFMT. Only two types of 
C           output are allowed in MOBILE6 
C  4 Jun 98 @DynTel-mla 2-659  Removed LOCFLG since the LAP record is 
C           obsolete in Mobile6. 
C 
      INTEGER ATPFLG,TPDFLG,RLFLAG,CAAFLG,CAAMY,POST96IM 
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,CAAFLG,CAAMY,POST96IM 
