
        MODULE MODUNCERT

!***********************************************************************
!  Module body starts at line
!
!  DESCRIPTION:
!     This module contains the public data used for the uncertainty reader
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 8/2001 by A. Holland
!     Modified 5/02 by G. Cano
!
!***************************************************************************
!
! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
!                System
! File: @(#)$Id$
!
! COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
! All Rights Reserved
!
! See file COPYRIGHT for conditions of use.
!
! Environmental Programs Group
! MCNC--North Carolina Supercomputing Center
! P.O. Box 12889
! Research Triangle Park, NC  27709-2889
!
! env_progs@mcnc.org
!
! Pathname: $Source$
! Last updated: $Date$ 
!
!****************************************************************************

        INCLUDE 'EMPRVT3.EXT'

!.........  Packet-specific parameters
        INTEGER, PARAMETER, PUBLIC :: NUCPCKT = 3   ! no. of uncert pkts.
        INTEGER, PARAMETER, PUBLIC :: UCPKTLEN  = 19! max length of pkt name
        INTEGER, PARAMETER, PUBLIC :: FA_IDX  = 1   ! pos. factor assign. pkt
        INTEGER, PARAMETER, PUBLIC :: EMP_IDX  = 2  ! pos. empirical pkt in mstr
        INTEGER, PARAMETER, PUBLIC :: PAR_IDX  = 3  ! pos. parametric pkt

        CHARACTER(LEN=UCPKTLEN), PARAMETER, PUBLIC :: 
     &                          UCPCKTS( NUCPCKT ) = 
     &                                  ( / '/FACTOR ASSIGNMENT/',
     &                                      '/EMPIRICAL/        ',
     &                                      '/PARAMETRIC/       '  / )

!.........  Define types needed for module
        TYPE :: FATYPE

            SEQUENCE

            INTEGER       :: INDX                   ! position of method pkt.

            CHARACTER(LEN=FIPLEN3) :: CFIP          ! FIPS code
            CHARACTER(LEN=SCCLEN3) :: TSCC          ! source category code
            CHARACTER(LEN=IOVLEN3) :: CPOL          ! pollutant name
            CHARACTER*1            :: METH          ! method
	    CHARACTER(LEN=DSCLEN3) :: NDIST	    ! name of distribution
	    CHARACTER*2		   :: APRCH	    ! approach
	    CHARACTER(LEN=PLTLEN3) :: PLT	    ! plant ID
	    CHARACTER(LEN=CHRLEN3) :: CHAR1	    ! plant characteristic 1
            CHARACTER(LEN=CHRLEN3) :: CHAR2         ! plant characteristic 2
            CHARACTER(LEN=CHRLEN3) :: CHAR3         ! plant characteristic 3
            CHARACTER(LEN=CHRLEN3) :: CHAR4         ! plant characteristic 4
            CHARACTER(LEN=CHRLEN3) :: CHAR5         ! plant characteristic 5
	    
	END TYPE
	
	TYPE :: EMPTYPE  
	
	    INTEGER        :: EPKTENT		    ! no. of emp. pkt. entries   
	    
	    CHARACTER*2	  	   :: ETYPE	    ! empirical dist. type
	    CHARACTER(LEN=DSCLEN3) :: EMPNAM	    ! empirical dist. name

        END TYPE
	
	TYPE :: PARTYPE
	
	    INTEGER        :: NUMP		    ! number of parameters
	
	    CHARACTER*1		   :: PTYPE	    ! parametric dist. type
	    CHARACTER(LEN=DSCLEN3) :: PARNAM	    ! parametric dist. name
	    
        END TYPE	

!.........  UNCERT file characteristics
        INTEGER, PUBLIC :: MXEMPDAT = 0   ! max no. of emission factor vals.
        INTEGER, PUBLIC :: MXPARDAT = 0   ! max no. of parameters
        INTEGER, PUBLIC :: NLINE_UC = 0   ! no. of lines in the file
        INTEGER, PUBLIC :: NFPCKT  = 0    ! no. of factor assignment pkts.
        INTEGER, PUBLIC :: NEPCKT  = 0    ! no. of empirical packets
        INTEGER, PUBLIC :: NPPCKT  = 0    ! no. of parametric packets 
	INTEGER, PUBLIC :: FPKTENT = 0	  ! no. of factor assignment entries

        LOGICAL, PUBLIC :: UC_ERROR = .FALSE.  ! true: error found reading file
        LOGICAL, PUBLIC, ALLOCATABLE :: USEPOLL( : )  !  true: pol in pkt

!.........  UNCERT characteristics arrays

        TYPE( FATYPE ) , PUBLIC, ALLOCATABLE :: FAPCKT( : )  ! fac. assig. pkt.
        TYPE( EMPTYPE ), PUBLIC, ALLOCATABLE :: EMPPCKT( : ) ! empirical pkt.
        TYPE( PARTYPE ), PUBLIC, ALLOCATABLE :: PARPCKT( : ) ! parametric pkt.

	INTEGER, PUBLIC, ALLOCATABLE :: FAINDX( : )       ! fac. assig. ent. no.

	REAL, PUBLIC, ALLOCATABLE :: EFVAL ( :, : )       ! emission fac. value
	REAL, PUBLIC, ALLOCATABLE :: PROB ( :, : )        ! probability
	REAL, PUBLIC, ALLOCATABLE :: PARAMET ( :, : )     ! parameters

!.........  Temporary packet-specific settings

	TYPE( FATYPE ) , PUBLIC :: FA_
	TYPE( EMPTYPE ), PUBLIC :: EMP_
	TYPE( PARTYPE ), PUBLIC :: PAR_

        INTEGER, PUBLIC :: PKT_IDX  = 0       ! index to ALLPCKTS for current
        INTEGER, PUBLIC :: PKTEND   = 0       ! ending line number of packet
        INTEGER, PUBLIC :: PKTSTART = 0       ! starting line number of packet
	INTEGER, PUBLIC :: EMPENTN  = 0	      ! empirical pkt. entry number
	INTEGER, PUBLIC :: FAENTN   = 0       ! factor assignment entry number

	REAL, PUBLIC :: EMISFAC = 0	      ! tmp. emis. fac. value
	REAL, PUBLIC :: EMPROB = 0	      ! tmp. probablity	
	REAL, PUBLIC :: PARA ( 100 )          ! tmp. parameters

        LOGICAL, PUBLIC :: INPACKET = .FALSE. ! true: currently in a packet
	LOGICAL, PUBLIC :: INFAPKT  = .FALSE. ! true: in factor assignment pkt
	LOGICAL, PUBLIC :: INEMPPKT = .FALSE. ! true: in empirical pkt
	LOGICAL, PUBLIC :: INPARPKT = .FALSE. ! true: in parametric pkt
	LOGICAL, PUBLIC :: EMPSTART = .FALSE. ! true: at start of empirical pkt
	LOGICAL, PUBLIC :: FASTART  = .FALSE. ! true: at start of factor pkt
	
        INTEGER, PUBLIC :: PKTCOUNT ( NUCPCKT ) ! no. of packets of given type
        LOGICAL, PUBLIC :: PKTSTATUS( NUCPCKT ) ! currently active packet

        CHARACTER(LEN=UCPKTLEN), PUBLIC :: PCKTNAM
        
!.........  Output file variables

	CHARACTER(LEN=IOVLEN3), PUBLIC, ALLOCATABLE :: UONAMES( : )  ! all possible variable names
        CHARACTER(LEN=IOVLEN3), PUBLIC, ALLOCATABLE :: UNAMES( : )   ! uncert var names

        INTEGER                         UPINVAR        ! no. of inventory vars from pointer file
        INTEGER                         NUOVAR         ! no. of total possible variables
        INTEGER                      :: UCOUNT = 0     ! sources with uncertainty

        INTEGER, PUBLIC, ALLOCATABLE :: APRCH( :, : )  ! approach
        INTEGER, PUBLIC, ALLOCATABLE :: EPTYP( :, : )  ! parametric/emp. type
        INTEGER, PUBLIC, ALLOCATABLE :: INSRC( : )     ! uncertainty sources
        INTEGER, PUBLIC, ALLOCATABLE :: METHOD( :, : ) ! empirical or parametric
        INTEGER, PUBLIC, ALLOCATABLE :: NUMEP( :, : )  ! no. of parameters/
        INTEGER, PUBLIC, ALLOCATABLE :: SRCNUM( : )    ! source number array
        INTEGER, PUBLIC, ALLOCATABLE :: UNCIDX( :, : ) ! empirical/parametric reference no.
                                                       ! empirical entires         
        REAL,    PUBLIC, ALLOCATABLE :: PARMS( :, : )  ! parameters 
        REAL,    PUBLIC, ALLOCATABLE :: EMFVAL( :, : ) ! emission factor values
        REAL,    PUBLIC, ALLOCATABLE :: PROBVAL( :, : )! probability values

!.........  Other variables
	LOGICAL, PUBLIC              :: UCFLAG = .FALSE. ! true: use uncertainty
	LOGICAL, PUBLIC              :: GUCFLAG = .FALSE.! true: use GRDMAT uncertainty
	LOGICAL, PUBLIC              :: TUCFLAG = .FALSE.! true: use TEMPORAL uncertainty
        
        LOGICAL, PUBLIC, ALLOCATABLE :: USTAT( : )     ! true: source is affected

        CHARACTER*4                  :: UCAT = ' '     ! category during uncertainty runs

        CHARACTER*9, PARAMETER :: PDF( 5 ) = ( / 'NORMAL   ',
     &                                           'LOGNORMAL',
     &                                           'GAMMA    ',
     &                                           'WEIBULL  ',
     &                                           'BETA     ' / )
        CHARACTER*16, PUBLIC :: UTNAME( 7 )            ! area or point temporal input

        LOGICAL, ALLOCATABLE, PUBLIC :: UACTVTY( : )   ! T:uncertainty activities
        LOGICAL, ALLOCATABLE, PUBLIC :: UEANAM( : )    ! T:uncertainty pollutants
        LOGICAL, ALLOCATABLE, PUBLIC :: UEAREAD( : )   ! T:uncertainty polls. & acts.
        LOGICAL, ALLOCATABLE, PUBLIC :: UEINAM( : )    ! T:uncertainty pollutants
        LOGICAL, ALLOCATABLE, PUBLIC :: UEFINAM( : )   ! T:uncertainty emission factors

        INTEGER, PARAMETER,   PUBLIC :: STATDIST = 5   ! Normal, Lognormal, Gamma, Weibull, Beta, resp.
        INTEGER, PARAMETER,   PUBLIC :: RMXLEN = 4     ! max number of realization is 9999

        INTEGER, PUBLIC, ALLOCATABLE :: NU_EXIST( :,: ) ! cert data to read in uncert mode
        INTEGER, PUBLIC, ALLOCATABLE :: UI_EXIST( :,: ) ! uncertainty data exists, read

        INTEGER                      :: DRAWRLZN = 0   ! realizations to draw (samples)
        INTEGER                         NEMP           ! number of empirical statistics options
        INTEGER                         NPAR           ! number of parametric statistics options
        INTEGER                         NPRB           ! number of probability value entries
        INTEGER                      :: UNIACT = 0     ! number of uncertainty activities
        INTEGER                      :: UNIPOL = 0     ! number of uncertainty pollutants
        INTEGER                         UNIPPA         ! number of uncertainty polls. & acts.
        INTEGER                      :: UNPACT = 0     ! number of variables per uncertainy activity
        INTEGER                      :: UNPPOL = 0     ! number of uncertainty pollutants & activities
        INTEGER                      :: UNSRC = 0      ! sources with uncertainty        
        INTEGER                      :: UNVAR = 0      ! number of uncertainty variables
        INTEGER                         VEMP           ! number of emp. stat. variables
        INTEGER                         VPAR           ! number of PAR. stat. variables
        INTEGER                         VPRB           ! number of probability value variables

        REAL,    PUBLIC, ALLOCATABLE :: SAMPGEN( :, : )! numbers generated for sampling


        END MODULE MODUNCERT
