C Version "@(#)$Id$ $Source$ $Date$ 

C  This is the location of the data for non-FTP effects
C
C  Common Block Array Subscripts:
C
C  REDUCT(3)       -  REDUCT ( IP )
C  NFTPV(4)        -  NFTPV  ( IP, 4=METHANE )
C  PHASSC(0:6)     -  PHASSC ( IYEAR )
C
C  Common Block Dictionary:
C
C   Name   Type              Description
C  ------  ----  --------------------------------------------------
C  NFTPYR   I    First model year of phase-in of new certification
C  NONFTP   I    Flag to indicate user-supplied driving behavior
C                1=no non-FTP effects
C                2=default non-FTP effects
C                3=user-supplied driving behavior for non-FTP effects
C  REDUCT   R    Reduction of non-FTP emission effects during phase-in
C  NFTPV    R    non-FTP emission effects
C  PHASSC   R    Phase-in schedule
C
C
C  Notes:
C
C  Aug-18-1994 @ CSC-bsg Request 462 NFTPDT was added for MOBILE5b to get the 
C  value of non-FTP offset.  Uses the variable NFTPYR from the user when the 
C  NONFTP Flag is 3, NONFTP is set in GETVER/HOR and depends on PRTFLG as read 
C  from the user.
C  Sep-20-1994 @ CSC-bsg Request 467 non-FTP offset number are stored in
C  COMMON.  In-use driving now depends on REDUCT, and PHASSC as well as
C  its original input parameters: vehicle speed etc...
C  Sep-30-1994 @ CSC-bsg Changed the IYEAR definitions and the common variable
C  PHASSC to go from 0 to 6 instead of 1 to 7.
C  Oct-6-1994 @ CSC-bsg Request 467 Comments were added.
C
C
      INTEGER NFTPYR,NONFTP
      REAL REDUCT,NFTPV,PHASSC
      COMMON /NFTPDT/ NFTPYR,NONFTP,REDUCT(3),NFTPV(4),PHASSC(0:6)
