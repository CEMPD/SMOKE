C Version "@(#)$Id$ $Source$ $Date$ 

C   This is the location for LEV Program Parameter Record.
C
C   Common Block Dictionary:
C
C    Name    Type                  Description
C   ------   ----  -----------------------------------------------------
C   LEVFLG    I    Flag indicating LEV user specification.
C                  1=none
C                  2=yes
C   IMFLG     I    Flag indicating whether appropriate I/M program is
C                  to be assumed in modeling.
C                  1=do not follow guidance memo I/M credits
C                  2=follow guidance memo I/M credits
C   LEVYRS    I    LEV Years as read in when Region=4.
C   APPLEV    I    Select application of LEV stds for LDGT2.
C                  'b', 1 => do not apply LEV standards to LDGT2
C                       2 => Apply LEV standards to LDGT2.
C
C   Notes:
C   
C   October-3-1994 @CSC-pme Request 2-472.  Variable APPLEV was added.
C
      INTEGER LEVFLG,IMFLG,LEVYRS,APPLEV
      COMMON /LEVBLK/ LEVFLG,IMFLG,LEVYRS,APPLEV                                                               
