C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION PCLEFT(MY,ICY,IP,IV)
C
C  PCLEFT determines the basic emission factor multiplicative
C  adjustment ( = percentage left of BEF = PCLEFT ) for the
C  effects of an inspection / maintanance (I/M) program, after
C  checking whether or not I/M applies to the factor being
C  computed by the calling function.
C
C  Called by BEF, BEFID1, BEFID2, BEFID3, BEFID4, BEFID5,
C  BEFID6, BEFIDL, and GETLEV.
C
C  Calls PCLFT1, and PCLFT2
C
C  Input on call:
C
C    parameter list: MY,ICY,IP,IV
C
C  Output on return:
C
C    function: PCLEFT
C
C
C
C  Notes:
C
C  Non-TECH4+ vehicles use 19th value for 20-24 benefit year.
C  TECH4+ vehicles have arrays expanded to MAXYRS-1.
C  PCLEFT was modified to include the IM240 check, and NOx
C  credits.
C  PCLEFT was modified to fix the index into I/M credits.
C  PCLEFT was modified to fix the biennial phase-in for
C  MOBILE5 corrected.  PCLEFT was modified again and turned into
C  a driver routine to decide to return the value of PCLFT1 or
C  PCLFT2 depending or ICY.  The old code is maintained in
C  PCLFT2.
C
C
      PCLEFT=PCLFT1(MY,ICY,IP,IV)
C
      RETURN
      END
