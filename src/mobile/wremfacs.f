
        SUBROUTINE WREMFACS( FNAME, NUMSRC, SDATE )

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC, ONLY: IRCLAS, IVTYPE

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC
        
C...........   This module contains emission factor tables and related
        USE MODEMFAC, ONLY: NEFS, EFSNAM, EFSDSC
        
        USE MODMBSET, ONLY: SCENLIST, EMISSIONS, EMISPOS, NONE
        
        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

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
        
        DO IHR = 1,24
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
                
                DO IEFTYPE = 1,8
 
C.....................  Set facility type to NONE for all ef types
C                       except exhaust running and evp running                
                    IF( IEFTYPE == 1 .OR. IEFTYPE == 6 ) THEN
                        CURRFTYPE = FTYPE
                    ELSE
                    	CURRFTYPE = NONE
                    END IF 
                    
                    DO IPOL = 1,3

C.........................  Skip non-exhaust EFs for non-HC pollutants
                        IF( IPOL /= 1 .AND. IEFTYPE > 2 ) CYCLE

C.........................  Find position in emissions array
                        EMISPOS = EFPOSITION( SCENNUM, IPOL, VTYPE, 
     &                                        IEFTYPE, CURRFTYPE, IHR )

C.........................  Determine position of EF/pollutant combo                    
                        IF( IEFTYPE == 1 .OR. IEFTYPE == 2 ) THEN
                            POLEF = ( IEFTYPE - 1 )*3 + IPOL
                        ELSE
                            POLEF = IEFTYPE + 4
                        END IF
                        
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
                CASE DEFAULT
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
                CASE DEFAULT
                    MESG = 'Unrecognized emission process ' //
     &                     EFSNAM( I )( 1:L-1 )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END SELECT
                
C.................  Determine position of EF/pollutant combo                    
                IF( VAREMIS == 1 .OR. VAREMIS == 2 ) THEN
                    POLEF = ( VAREMIS - 1 )*3 + VARPOL
                ELSE
                    POLEF = VAREMIS + 4
                END IF                
                
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