C Version @(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD24
C
C  BLOCK DATA Subprogram 24: /NAMES1/,/STRING/,/USDATA/
C
C  BLOCK DATA 24 initializes character strings used in prompts & output formats
C  and zeroes pointers into character string arrays.
C
C  Common block array subscripts:
C
C  NAMFLG(18)     -  NAMFLG ( IFG )
C  NAMMMR(2,3,8)  -  NAMMMR ( ICH, IG, IV )
C  NAMRVP(4)      -  NAMRVP ( IRVP )
C  NAMTEM(2)      -  NAMTEM ( ITEM )
C  NOYES(2)       -  NOYES  ( <ILDT,INO,IMFLAG,ATPFLG> )
C  USNAME(4,4)    -  USNAME ( ICH, IUSD )
C  IUSD(4)        -  IUSD   ( ICT )
C
C  Common block dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------
C  /NAMES1/:
C
C  NAMFLG   C*8  control flag names for prompts and error messages
C  NAMMMR   C*8  prompts for replacement AMAR & JULMYR, by JDX & IV
C  NAMRVP   C*8  name of RVP falling outside allowable range
C  NAMTEM   C*4  'mini' & 'maxi' prefixes for fuel tank temp error messages
C
C  /STRING/:
C
C  NOYES    C*4  character string for no/yes
C  COMMA    C*4  character string for comma
C  PERIOD   C*4  character string for period
C  COLON    C*1  character string for colon
C
C  /USDATA/:
C
C  USNAME   R    names of one-time data supplied by user
C  NUSD     I    number of names of user supplied one-time data for this run
C  IUSD     I    pointers to the names of user supplied one-time data for run
C
C  Notes:
C
C  BD 24 was modified in Version 01 to correct the definition order of USNAME
C  and IUSD.
C
C
      CHARACTER*8 NAMFLG,NAMMMR
      CHARACTER*8 USNAME
      CHARACTER*8 NAMRVP
      CHARACTER*4 NAMTEM
      CHARACTER*4 NOYES,COMMA,PERIOD
      CHARACTER*1 COLON
C
      COMMON /NAMES1/ NAMFLG(20),NAMMMR(2,3,8),NAMRVP(4),NAMTEM(2)
      COMMON /STRING/ NOYES(2),COMMA,PERIOD,COLON
      COMMON /USDATA/ USNAME(4,4),NUSD,IUSD(4)
C
C  /NAMES1/: These names are used in input routines and in QUITER.
C
      DATA NAMFLG /' TAMFLG ', ' SPDFLG ', ' VMFLAG ', ' MYMRFG ',
     1     ' NEWFLG ', ' IMFLAG ', ' ALHFLG ', ' ATPFLG ', ' RLFLAG ',
     2     ' LOCFLG ', ' TEMFLG ', ' OUTFMT ', ' PRTFLG ', ' IDLFLG ',
     3     ' NMHFLG ', ' HCFLAG ', ' COLDFG ', ' OXYFLG ', ' DSFLAG ',
     4     ' RFGLAG '/
C
      DATA NAMMMR/
     *'1-10 for',' LDGV:  ','11-20 fo','r LDGV: ','21-25 fo','r LDGV: ',
     *'1-10 for',' LDGT1: ','11-20 fo','r LDGT1:','21-25 fo','r LDGT1:',
     *'1-10 for',' LDGT2: ','11-20 fo','r LDGT2:','21-25 fo','r LDGT2:',
     *'1-10 for',' HDGV:  ','11-20 fo','r HDGV: ','21-25 fo','r HDGV: ',
     *'1-10 for',' LDDV:  ','11-20 fo','r LDDV: ','21-25 fo','r LDDV: ',
     *'1-10 for',' LDDT:  ','11-20 fo','r LDDT: ','21-25 fo','r LDDT: ',
     *'1-10 for',' HDDV:  ','11-20 fo','r HDDV: ','21-25 fo','r HDDV: ',
     *'1-10 for',' MC:    ','11-20 fo','r MC:   ','21-25 fo','r MC:   '/
C
      DATA NAMRVP /'Period 1','Period 2','Refuelng','Weatherd'/
C
      DATA NAMTEM/'mini','maxi'/
C
C  /STRING1/:
C
      DATA NOYES/'No  ','Yes '/
      DATA COMMA/',   '/,PERIOD/'.   '/,COLON/':'/
C
C  /USDATA/:
C
      DATA USNAME/' tamperi','ng and m','isfuelin','g rates ',
     *            '  basic ','exhaust ','emission','s rates ',
     *            '  mileag','e accrua','l distri','butions ',
     *            ' veh reg','istratio','n distri','butions '/
C
      DATA NUSD,IUSD/5*0/
C
      END
