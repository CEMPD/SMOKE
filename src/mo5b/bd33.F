C Version "@(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD33
C
C  BLOCK DATA Subprogram 33: /REGION/,/SCENE1/,/SCENE2/
C
C  BLOCK DATA 33 defines common blocks related to the first scenario record.
C
C  Common block array subscripts:
C
C  FEET(2)  -  FEET ( IREJN )
C  SPD(2)   -  SPD  ( IV )
C
C  Common block dictionary:
C
C   Name   Type              Description
C  ------  ----  --------------------------------------------------
C  /REGION/:
C
C  IREJN    I    user supplied region code: 1 = low altitude, 2 = high
C                1 = low altitude region
C                2 = high altitude region
C  ALT      R    altitude in feet
C  INITPR   I    switch: 1 = 1st scenario's output pass, print run parameters
C                        0 = 2nd or later, do not 'CALL OUTHD<OUTFMT>'
C
C  /SCENE1/:
C
C  SPD      R    user supplied average speeds of 8 vehicle types
C  PCCN     R    user supplied % in cold start mode by non-catalyst
C                equipped vehicle types
C  PCHC     R    user supplied % in hot start mode by catalyst
C                equipped vehicle types
C  PCCC     R    user supplied % in cold start mode by catalyst
C                equipped vehicle types
C
C  /SCENE2/:
C
C  FCC      R    fraction of vmt in cold start mode, vehicle catalyst equipped
C  FCN      R    fraction of vmt in cold start mode, non-catalyst equipped
C  FHC      R    fraction of vmt in hot start mode, catalyst equipped
C  FHN      R    fraction of vmt in hot start mode, non-catalyst equipped
C
C  Notes:
C
C  None.
C
C
      INCLUDE 'LEVBLK.I'
      INCLUDE 'REGION.I'
      INCLUDE 'SCENE1.I'
      INCLUDE 'SCENE2.I'
C
C  Only FEET and INITPR need to be initialized.
C
      DATA FEET/500.,5500./,INITPR/1/
C
C  The rest are read in or set during execution, so dummy zeroes are used here.
C
      DATA IREJN,ALT,LEVFLG,LEVYRS/3*0,1994/
      DATA SPD,PCCN,PCHC,PCCC/11*0.0/
      DATA FCC,FCN,FHC,FHN/4*0.0/
C
      END
