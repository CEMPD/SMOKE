C Version "@(#)$Id$ $Source$ $Date$ 

cmh LMOS change from Mobile5a
c   Change length of LEVNAM to allow full directory pathname
c
c     CHARACTER   LEVNAM*40
      CHARACTER   LEVNAM*256
      COMMON /LEVINF/ LEVNAM
