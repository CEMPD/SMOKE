C 
C  Changes: (Last change first) 
C 
C  7 jan 00 @DynTel-ddj  1-020  New include FILE for SFTP arrays. 
C 
C The initialization of the data could be found in BD_SFTP. 
C 
      REAL SFTP_FRAC, SFTP_BENE_AC, SFTP_BENE 
C 
      COMMON /SFTP / SFTP_FRAC(2,4), SFTP_BENE_AC(2,2,5,2), 
     &               SFTP_BENE(2,3,5,2) 
