
        SUBROUTINE WRIDAOUT( CATEGORY, FILNAM, FILUNIT, NSRC, NSEG, 
     &                       POLBUF, STATUS )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      Read and write source characteristic information to IDA output file.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C**************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS:
        INTEGER         FIND1
        INTEGER         GETIFDSC
        INTEGER         TRIMLEN

        EXTERNAL        FIND1, GETIFDSC, TRIMLEN

C...........   SUBROUTINE ARGUMENTS
        CHARACTER*(*)   CATEGORY           ! Source category
        CHARACTER*(*)   FILNAM             ! File name
        INTEGER         FILUNIT            ! Unit number for character part
        INTEGER         NSRC               ! Number of sources
        INTEGER         NSEG               ! Number of segments of pollutant buf
        CHARACTER*(*)   POLBUF( NSRC,NSEG )! Pollutant buffer
        INTEGER         STATUS             ! Exit status

C...........   Local parameters, indpendent
        INTEGER, PARAMETER :: PLANTLEN = 15
        INTEGER, PARAMETER :: POINTLEN = 15
        INTEGER, PARAMETER :: STACKLEN = 12
        INTEGER, PARAMETER :: BOILRLEN = 5
        INTEGER, PARAMETER :: SEGMTLEN = 2
        INTEGER, PARAMETER :: PDESCLEN = 40
        INTEGER, PARAMETER :: SCCLEN   = 10

        LOGICAL IDACOLS( 7 )

C...........   Common (all source categories) source characteristics
        INTEGER, ALLOCATABLE:: IFIP  ( : )  !  source FIPS state/county codes
        INTEGER, ALLOCATABLE:: INTDUM( : )  !  dummy integer array

C...........   Area source characteristics
C NOTE: Fill in later for area sources

C...........   Mobile source characteristics
C NOTE: Fill in later for mobile sources

C...........   Point source characteristics
        INTEGER, ALLOCATABLE:: ISCC  ( : )  !  source SCC
        INTEGER, ALLOCATABLE:: ISIC  ( : )  !  source SIC
        INTEGER, ALLOCATABLE:: IORIS ( : )  !  source ORIS ID code
        INTEGER, ALLOCATABLE:: INVYR ( : )  !  inventory year
        REAL   , ALLOCATABLE:: XLOCA ( : )  !  UTM X-location (m)
        REAL   , ALLOCATABLE:: YLOCA ( : )  !  UTM Y-location (m)
        REAL   , ALLOCATABLE:: STKHT ( : )  !  stack height   (m)
        REAL   , ALLOCATABLE:: STKDM ( : )  !  stack diameter (m)
        REAL   , ALLOCATABLE:: STKTK ( : )  !  exhaust temperature (deg K)
        REAL   , ALLOCATABLE:: STKVE ( : )  !  exhaust velocity    (m/s)

        CHARACTER(LEN=BLRLEN3), ALLOCATABLE:: CBLRID( : ) !  boiler ID
        CHARACTER(LEN=DSCLEN3), ALLOCATABLE:: CPDESC( : ) !  plant description

        CHARACTER(LEN=SRCLEN3), ALLOCATABLE:: CSOURC( : ) !  concatonated source

C...........   IDA output variables (names same as IDA format description)

        INTEGER         STID, CYID, BEGYR, ENDYR, HOURS, START, DAYS
        INTEGER         WEEKS, ISIC

        REAL            STKHGT, STKDIAM, STKTEMP, STKFLOW, STKVEL
        REAL            BOILCAP, WINTHRU, SPRTHRU, SUMTHRU, FALTHRU
        REAL            THRUPUT, MAXRATE, HEATCON, SULFCON, ASHCON
        REAL            NETDC, LATC, LONC

        CHARACTER(LEN=1)        CAPUNITS, OFFSHORE  
        CHARACTER(LEN=PLANTLEN) PLANTID  
        CHARACTER(LEN=POINTLEN) POINTID
        CHARACTER(LEN=STACKLEN) STACKID
        CHARACTER(LEN=BOILRLEN) BLRID
        CHARACTER(LEN=SEGMTLEN) SEGMENT
        CHARACTER(LEN=PDESCLEN) PLANT
        CHARACTER(LEN=SCCLEN  ) SCC

C...........   Other local variables

        INTEGER         L, L1, L2, S        ! counters and indices

        INTEGER         COID     ! tmp country code
        INTEGER         FIP      ! tmp FIPS state and county code
        INTEGER         IOS      ! i/o status
        INTEGER         IS, IE   ! start & end position of pols in VNAMES3D
        INTEGER         LMAX     ! maximum length of POLINE
        INTEGER         LCOID    ! previous country ID in loop
        INTEGER         LYEAR    ! previous year in loop
        INTEGER         NCHAR    ! number of strings returned from PARSCSRC
        INTEGER         NNONPV   ! number of non-pollutant-specific variables
        INTEGER         NP       ! actual number of entries in POLINE
        INTEGER         NVPERP   ! number of variables per pollutant
        INTEGER         YEAR     ! tmp 4-digit year

        REAL            CTOF     ! celius to farenheit
        REAL            M2FT     ! meters to feet

        CHARACTER*65, ALLOCATABLE:: POLINE( : )   !  pollutant list strings
        CHARACTER*100 CHARS( 7 )     !  source fields for output
        CHARACTER*300 MESG           !  message buffer

        LOGICAL IDACOLS( 7 )
        DATA    IDACOLS / 6*.TRUE., .FALSE. /

        CHARACTER*16 :: PROGNAME = 'WRIDAPOL' ! program name

C***********************************************************************
C   begin body of subroutine WRIDAPOL

        STATUS = 0

C.........  Get file header from I/O API NetCDF inventory file
        IF( .NOT. DESC3( FILNAM ) ) THEN

            STATUS = 1
            MESG = 'Could not read description for "' //
     &             FILNAM( 1:TRIMLEN( FILNAM ) ) // '"'
            CALL M3MSG2( MESG )
            RETURN
            
        ENDIF

C.........  Write pollutant names to single buffer for header
        NNONPV = GETIFDSC( FDESC3, '/NON POLLUTANT/' )
        NVPERP = GETIFDSC( FDESC3, '/PER POLLUTANT/' )

        LMAX = LEN( POLINE( 1 ) )  ! length of POLINE string
        IS   = NONPV + 1          
        IE   = MXVARS3
        L    = 0
        NP   = 0
        DO I = IS, IE, NVPERP      ! Find out how many rows

            L1 = TRIMLEN( VNAME3D( I ) ) + 1
            IF( L + L1 .GT. LMAX ) THEN
               NP = NP + 1
               L = 0
            ENDIF
            L = L + L1

        ENDDO

        ALLOCATE( POLINE( NP ), STAT=IOS ) ! Allocate memory 
        CALL CHECKMEM( IOS, 'POLINE', PROGNAME )  

        NP = 1                             ! Fill in POLINE
        POLINE( 1 ) = VNAME3D( IS )
        DO I = IS + NVPERP, IE, NVPERP

            L1 = TRIMLEN( VNAME3D( I ) ) + 1
            L2 = TRIMLEN( POLINE( NP ) )
            IF( L2 + L1 .GT. LMAX ) THEN
               NP = NP + 1
            ENDIF

            POLINE( NP ) = POLINE ( NP ) ( 1:L2 ) // ' ' // 
     &                     VNAME3D( I )  ( 1:L1 )
        ENDDO   

C.........  Allocate memory for shared arrays
        ALLOCATE( IFIP( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IFIP', PROGNAME )  
        ALLOCATE( INTDUM( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INTDUM', PROGNAME )  

C.........  Processing for area sources
        IF( CATEGORY .EQ. 'AREA' ) THEN

C..............  Allocate memory 

C..............  Read area-source characteristics 
C            CALL RAREACHR( FILNAM, NSRC, IFIP, ??? )

C..............  Write area-source characteristics to output file

C.........  Allocate memory for mobile-source arrays and read from file
        ELSEIF( CATEGORY .EQ. 'MOBILE' ) THEN

C..............  Allocate memory 

C..............  Read mobile-source characteristics 
C            CALL RMOBLCHR( FILNAM, NSRC, IFIP, ??? )

C..............  Write mobile-source characteristics to output file 

C.........  Allocate memory for point-source arrays and read from file
        ELSEIF( CATEGORY .EQ. 'POINT' ) THEN

C..............  Allocate memory 
            ALLOCATE( ISCC( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ISCC', PROGNAME )  
            ALLOCATE( ISIC( NPSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ISIC', PROGNAME )
            ALLOCATE( IORIS( NPSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'IORIS', PROGNAME )
            ALLOCATE( INVYR( NPSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INVYR', PROGNAME )
            ALLOCATE( XLOCA( NPSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'XLOCA', PROGNAME )
            ALLOCATE( YLOCA( NPSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'YLOCA', PROGNAME )
            ALLOCATE( STKHT( NPSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'STKHT', PROGNAME )
            ALLOCATE( STKDM( NPSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'STKDM', PROGNAME )
            ALLOCATE( STKTK( NPSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'STKTK', PROGNAME )
            ALLOCATE( STKVE( NPSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'STKVE', PROGNAME )
            ALLOCATE( CBLRID( NPSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CBLRID', PROGNAME )  
            ALLOCATE( CPDESC( NPSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CPDESC', PROGNAME )  
            ALLOCATE( CSOURC( NPSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CSOURC', PROGNAME )  

C..............  Read point-source characteristics 
            CALL RPNTSCHR( FILNAM, FILUNIT, NSRC, IFIP, ISCC, ISIC, 
     &                     IORIS, INTDUM, INTDUM, INVYR, XLOCA, 
     &                     YLOCA, STKHT, STKDM, STKTK, STKVE,
     &                     CBLRID, CPDESC, CSOURC )

C.............  Set values of variables that SMOKE does not (yet) support
            BEGYR    = -9
            ENDYR    = -9
            BOILCAP  = -9.
            CAPUNITS = ' '
            WINTHRU  = 0.
            SPRTHRU  = 0.
            SUMTHRU  = 0.
            FALTHRU  = 0.
            HOURS    = 0
            START    = 0
            DAYS     = 0
            WEEKS    = 0
            THRUPUT  = -9.
            MAXRATE  = -9.
            HEATCON  = -9.
            SULFCON  = -9.
            ASHCON   = -9.
            NETDC    = -9.
            OFFSHORE = ' '

C.............  Compute conversion constants
            M2FT  = 1./FT2M
            CTOF  = 1./FTOC

C.............  Write point-source characteristics to output file
            LCOID = -9
            LYEAR = -9
            DO S = 1, NSRC
                CALL PARSCSRC( CSOURC( S ), IDACOLS, CHARS, NCHAR )

C.................  Truncate character string variables
                PLANTID  = CHARS( 2 ) 
                POINTID  = CHARS( 3 )
                STACKID  = CHARS( 4 )
                SEGMENT  = CHARS( 5 )
                SCC      = CHARS( 6 )

                BLRID    = CBLRID( S )
                PLANT    = CPDESC( S )

C.................  Convert units of stack parameters
                STKHGT  = STKHT( S ) * M2FT
                STKDIAM = STKDM( S ) * M2FT
                STKTEMP = ( STKTK( S ) - CTOK ) / CTOF + 32.
                STKVEL  = STKVE( S ) * M2FT
                STKFLOW = STKVEL * 0.25 * PI * STKDIAM * STKDIAM

C.................  Store others in temporary variables
                COID = IFIP( S ) / 10000
                FIP  = IFIP( S ) - COID * 10000
                STID = FIP / 1000 
                CYID = FIP - STID * 1000

                WRITE( CORIS, '(I6)' ) IORIS( S )

                SIC  = ISIC ( S )                
                YEAR = INVYR( S )                
                LATC = YLOCA( S )                
                LONC = XLOCA( S )                

C.................  Write out header
                IF( COID .NE. LCOID .OR. YEAR .NE. LYEAR ) THEN
                    LCOID = COID
                    LYEAR = YEAR

                    WRITE( CYEAR, '(I4)' ) YEAR

                    K = FIND1( COID, MXCNTRY3, CNTRYCD3 )

                    IF( K .GT. 0 ) THEN

                        WRITE( OUTUNIT, 93000 ) 
     &                     '#TYPE     Point Source Inventory',
     &                     '#COUNTRY  ' // CNTRYNM3( K ),
     &                     '#YEAR     ' // CYEAR,
     &                     '#DESC     Output from SMOKE',
     &                     '#POLID    ' POLINE( 1 ),
     &                   ( POLINE( I ), I = 2, NP )

                    ELSE
                        STATUS = 1
                        WRITE( MESG,94010 ) 'Invalid country code', K,
     &                         'found at source', S
                        CALL M3MESG( MESG )
                        GO TO 101          ! To end of loop (skip write)

                    ENDIF
                ENDIF

C.................  Write out main entries
                WRITE( OUTUNIT, 93340 ) STID, CYID, PLANT, POINT,
     &                 STACK, CORIS, BOILER, SEGMENT, PLNTDESC, SCC,
     &                 BEGYR, ENDYR, STKHGT, STKDIAM, STKTEMP, STKFLOW, 
     &                 STKVEL, BOILCAP, CAPUNITS, WINTHRU, SPRTHRU, 
     &                 FALTHRU, HOURS, START, DAYS, WEEKS, THRUPUT, 
     &                 MAXRATE, HEATCON, SULFCON, ASHCON, NETDC, SIC,
     &                 LATC, LONC, OFFSHORE, 
     &               ( POLBUF( S,I ), I = 1, NSEG )

101         ENDDO

        ENDIF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( 500A )

C93320   FORMAT(  area

C93330   FORMAT(  mobile

93340   FORMAT( I2, I3, A15, A15, A12, A6, A5, A2, A40, A10, I4, I4,    ! point
     &          F4.0, F6.2, F4.0, F10.2, F9.2, F8.2, A1, 4F2.0, 2I2, I1,
     &          I2, F11.1, F12.3, F8.2, F5.2, F5.2, F9.3, I4, 2F9.4, A1,
     &          <NSEG>( A ) )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END

