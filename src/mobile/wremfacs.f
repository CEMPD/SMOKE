
        SUBROUTINE WREMFACS( FNAME, NUMSRC, SDATE )

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC, ONLY: IRCLAS, IVTYPE

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC
        
C...........   This module contains emission factor tables and related
        USE MODEMFAC, ONLY: NEFS, EFSNAM, EFSDSC
        
        USE MODMBSET, ONLY: SCENLIST, EMISSIONS, M6NONE
        
        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'M6CNST3.EXT'   !  Mobile6 constants
        
C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER         CVTRDTYPE
        INTEGER         CVTVEHTYPE
        INTEGER         EFPOSITION
        
        EXTERNAL   CVTRDTYPE, CVTVEHTYPE, EFPOSITION

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: FNAME    ! logical name of emission factors file
        INTEGER,      INTENT (IN) :: NUMSRC   ! total number of sources
        INTEGER,      INTENT (IN) :: SDATE    ! episode start date

C...........   LOCAL VARIABLES and their descriptions:

C...........   Local allocatable arrays
        REAL,    ALLOCATABLE :: SRCEFS( :,: )  ! per-source emission factors
        INTEGER, ALLOCATABLE :: EFIDX( : )     ! index of source numbers

C.........   Other local variables
        INTEGER    I, L, L2, LJ     ! counters and indices        
        INTEGER    IHR, ISRC        ! indices for assigning efs to sources
        INTEGER    IEFTYPE          ! ditto
        INTEGER    IPOL             ! ditto
        INTEGER    IOS              ! i/o status
        INTEGER    SCENNUM          ! scenario number of source
        INTEGER    VTYPE            ! vehicle type of source
        INTEGER    FTYPE            ! facility type of source
        INTEGER    CURRFTYPE        ! facility type for ef
        INTEGER    JDATE            ! output date
        INTEGER    JTIME            ! output time
        INTEGER    EFPOS            ! current position in emission factor arrays
        INTEGER    EMISPOS          ! position in master emission factor array
        INTEGER    POLEF            ! ef/pollutant combo position in SRCEFS
        INTEGER    VARPOL           ! pollutant specified by current variable
        INTEGER    VAREMIS          ! emission process specified by current variable
        
        LOGICAL       :: LASAFLAG = .FALSE. ! true: treat local roads as arterial
        LOGICAL       :: SRCMATCH = .FALSE. ! true: matched source to emission factors
        LOGICAL, SAVE :: INITIAL = .TRUE.   ! true: first time through subroutine

        CHARACTER*300          MESG      !  message buffer 
        
        CHARACTER*16  :: PROGNAME = 'WREMFACS' ! program name
        
C***********************************************************************
C   begin body of program WREMFACS

        LJ = LEN_TRIM( ETJOIN )

C.........  Allocate arrays for per-source efs
        IF( INITIAL ) THEN     
            ALLOCATE( SRCEFS( NUMSRC, NEFS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EMISFACS', PROGNAME )
            ALLOCATE( EFIDX( NUMSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EFIDX', PROGNAME )
            
            INITIAL = .FALSE.
        END IF

        DO IHR = 1, 24
            EFPOS = 1
            
            SRCEFS    = 0.
            EFIDX     = 0
        
            DO ISRC = 1, NSRC
                SCENNUM = SCENLIST( ISRC,1 )
            
                IF( SCENNUM == 0 ) CYCLE 
                
C.................  Check local-as-arterial setting for this source                
                IF( SCENLIST( ISRC,2 ) == 2 ) THEN
                    LASAFLAG = .TRUE.
                ELSE
                    LASAFLAG = .FALSE.
                END IF
                    
                FTYPE = CVTRDTYPE( IRCLAS( ISRC ), LASAFLAG )
                VTYPE = CVTVEHTYPE( IVTYPE( ISRC ) )
                
                DO IEFTYPE = 1,MXM6EPR
 
C.....................  Set facility type to NONE for all ef types
C                       except exhaust running and evp running                
                    IF( IEFTYPE == 1 .OR. IEFTYPE == 6 ) THEN
                        CURRFTYPE = FTYPE
                    ELSE
                    	CURRFTYPE = M6NONE
                    END IF 
                    
                    DO IPOL = 1,MXM6POLS
                    
C.........................  Determine position of EF/pollutant combo                    
                        POLEF = M6POL2EF( IEFTYPE, IPOL )

C.........................  Make sure this is a valid combination
                        IF( POLEF == -1 ) CYCLE

C.........................  Find position in master emissions array
                        EMISPOS = EFPOSITION( SCENNUM, IPOL, VTYPE, 
     &                                        IEFTYPE, CURRFTYPE, 
     &                                        IHR, .FALSE. )
                      
                        SRCEFS( EFPOS, POLEF ) = EMISSIONS( EMISPOS )
                        
                    END DO    ! end pollutant loop              
                END DO    ! end emission factor loop
                
                EFIDX( EFPOS ) = ISRC	
                EFPOS = EFPOS + 1
            END DO     ! end source loop
 
C.............  Write emission factors to file
            JDATE = SDATE
            JTIME = ( IHR - 1 ) * 10000

            IF( .NOT. WRITE3( FNAME, 'SOURCES', JDATE, JTIME,
     &                        EFIDX ) ) THEN
     	        MESG = 'Could not write source numbers ' //
     &                 ' to "' // FNAME( 1:LEN_TRIM( FNAME ) ) // '".'
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            END IF

            DO I = 1, NEFS
                L  = INDEX( EFSNAM( I ), ETJOIN )
                L2 = LEN_TRIM( EFSNAM( I ) )

C.................  Match pollutant in variable name to number               
                SELECT CASE( EFSNAM( I )( L+LJ:L2 ) )
                CASE( 'CO' )
                    VARPOL = 2
                CASE( 'NOX' )
                    VARPOL = 3
                CASE( 'SO4' )
                    VARPOL = 4
                CASE( 'O25' )
                    VARPOL = 5
                CASE( 'E25' )
                    VARPOL = 6
                CASE( 'G25' )
                    VARPOL = 7
                CASE( 'SO2' )
                    VARPOL = 8
                CASE( 'NH3' )
                    VARPOL = 9
                CASE( 'B25' )
                    VARPOL = 10
                CASE( 'T25' )
                    VARPOL = 11
                CASE( 'BEN' )
                    VARPOL = 12
                CASE( 'MTB' )
                    VARPOL = 13
                CASE( 'BUT' )
                    VARPOL = 14
                CASE( 'FOR' )
                    VARPOL = 15
                CASE( 'ACE' )
                    VARPOL = 16
                CASE( 'ACR' )
                    VARPOL = 17
                CASE( 'O10' )
                    VARPOL = 18
                CASE( 'E10' )
                    VARPOL = 19
                CASE( 'G10' )
                    VARPOL = 20
                CASE( 'B10' )
                    VARPOL = 21
                CASE( 'T10' )
                    VARPOL = 22
                CASE DEFAULT   ! any variety of hydrocarbon
                    VARPOL = 1
                END SELECT

C.................  Match emission process to number                
                SELECT CASE( EFSNAM( I )( 1:L-1 ) )
                CASE( 'EXR' )
                    VAREMIS = 1
                CASE( 'EXS' )
                    VAREMIS = 2
                CASE( 'HOT' )
                    VAREMIS = 3
                CASE( 'DNL' )
                    VAREMIS = 4
                CASE( 'RST' )
                    VAREMIS = 5
                CASE( 'EVR' )
                    VAREMIS = 6
                CASE( 'CRC' )
                    VAREMIS = 7
                CASE( 'RFL' )
                    VAREMIS = 8
                CASE( 'BRK' )
                    VAREMIS = 9
                CASE( 'TIR' )
                    VAREMIS = 10
                CASE DEFAULT
                    MESG = 'Unrecognized emission process ' //
     &                     EFSNAM( I )( 1:L-1 )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END SELECT
                
C.................  Determine position of EF/pollutant combo                    
                POLEF = M6POL2EF( VAREMIS, VARPOL )     
                
                IF( .NOT. WRITE3( FNAME, EFSNAM( I )( 1:L2 ), JDATE, 
     &                            JTIME, SRCEFS( :,POLEF ) ) ) THEN
     	            MESG = 'Could not write ' // 
     &	                   EFSDSC( I )( 1:LEN_TRIM( EFSDSC( I ) ) ) // 
     &	                   'to "' // FNAME( 1:LEN_TRIM( FNAME ) ) // 
     &                     '".'
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF                
            END DO     ! end ef/pollutant combo loop
            
        END DO     ! end hour loop
        
        END SUBROUTINE WREMFACS