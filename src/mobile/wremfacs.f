
        SUBROUTINE WREMFACS( FNAME, NUMSRC, SDATE, VOLNAM )

C***********************************************************************
C  subroutine body starts at line 109
C
C  DESCRIPTION:
C       Finds correct emission factor for each source and writes 
C       emission factors to output files
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     10/01: Created by C. Seppanen
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***********************************************************************
C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC, ONLY: IRCLAS, IVTYPE

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC
        
C...........   This module contains emission factor tables and related
        USE MODEMFAC, ONLY: NEFS, MXETYPE, EMTNAM, NSUBPOL, SUBPOLS, 
     &                      OUTPUTHC, SCENLIST, EMISSIONS, NTOTHAPS, 
     &                      HAPNAMES, HAPEFS
        
        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'M6CNST3.EXT'   !  Mobile6 constants
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions
        
C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER         CVTRDTYPE
        INTEGER         CVTVEHTYPE
        INTEGER         INDEX1
        
        EXTERNAL   CVTRDTYPE, CVTVEHTYPE, INDEX1

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: FNAME    ! logical name of emission factors file
        INTEGER,      INTENT (IN) :: NUMSRC   ! total number of sources
        INTEGER,      INTENT (IN) :: SDATE    ! episode start date
        CHARACTER(*), INTENT (IN) :: VOLNAM   ! volatile pollutant name

C...........   LOCAL VARIABLES and their descriptions:

C...........   Local allocatable arrays
        REAL,    ALLOCATABLE, SAVE :: SRCEFS( : )   ! per-source emission factors
        INTEGER, ALLOCATABLE, SAVE :: EFIDX( : )    ! index of source numbers

C.........   Other local variables
        INTEGER    I, J, K          ! counters and indices        
        INTEGER    IHR              ! hour loop index
        INTEGER    ISRC             ! source loop index
        INTEGER    IEF              ! emission process loop index
        INTEGER    IPOL             ! pollutant loop index
        INTEGER    IOS              ! I/O status
        INTEGER    SCENNUM          ! scenario number of source
        INTEGER    VTYPE            ! vehicle type of source
        INTEGER    FTYPE            ! facility type of source
        INTEGER    JDATE            ! output date
        INTEGER    JTIME            ! output time
        INTEGER    EFPOS            ! current position in emission factor arrays
        INTEGER    EMISPOS          ! position in master emission factor array
        INTEGER    POLEF            ! ef/pollutant combo position in SRCEFS
        INTEGER    VARPOL           ! pollutant specified by current variable
        
        REAL       EFVAL            ! temporary emission factor value
        REAL       RAMPEF           ! ramp emission factor, used when creating freeway composite
        
        LOGICAL       :: RLASAFLAG = .FALSE.! true: treat rural local roads as arterial
        LOGICAL       :: ULASAFLAG = .FALSE.! true: treat urban local roads as arterial
        LOGICAL       :: USEHAP  = .FALSE.  ! true: use user-defined toxic value
        LOGICAL, SAVE :: INITIAL = .TRUE.   ! true: first time through subroutine

        CHARACTER(3)       EFNAME    !  emission process name
        CHARACTER(11)      POLLNAME  !  pollutant name
        CHARACTER(16)      CURRVAR   !  current variable being written
        CHARACTER(300)     MESG      !  message buffer 
        
        CHARACTER(16) :: PROGNAME = 'WREMFACS' ! program name
        
C***********************************************************************
C   begin body of subroutine WREMFACS

C.........  Allocate arrays for per-source efs
        IF( INITIAL ) THEN     
            ALLOCATE( SRCEFS( NUMSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EMISFACS', PROGNAME )
            ALLOCATE( EFIDX( NUMSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EFIDX', PROGNAME )
            
            INITIAL = .FALSE.
        END IF

        DO IHR = 1, 24  ! loop over hours
          DO IEF = 1, MXM6EPR  ! loop over emission processes
          
C.............  Set facility type to NONE for all ef types except
C               exhaust running, evp running, brake, and tire wear
            IF( M6FAC2EF( IEF, M6NONE ) == 1 ) THEN
                FTYPE = M6NONE
            ELSE
            	FTYPE = 1
            END IF 

C.............  Loop over all pollutants including user-defined                
            DO IPOL = 1,MXM6POLS + NTOTHAPS
            
C.................  Make sure this is a valid ef/pollutant combination
              IF( IPOL <= MXM6POLS ) THEN
                  IF( M6POL2EF( IEF, IPOL ) == -1 ) CYCLE
              ELSE

C...................  For user-defined toxics, skip CRC, BRK, and TIR
                  IF( IEF == 7 .OR. IEF == 9 .OR. IEF == 10 ) CYCLE
              END IF

C...............  Reset array position and ef arrays                
              EFPOS  = 0
              SRCEFS = 0.  ! array
              EFIDX  = 0   ! array
              
              DO ISRC = 1, NSRC  ! loop over sources

C.................  Get scenario number for current source        
                SCENNUM = SCENLIST( ISRC,1 )            
                IF( SCENNUM == 0 ) CYCLE 

C.................  Store array position for current source                
                EFPOS = EFPOS + 1
                EFIDX( EFPOS ) = ISRC
                
C.................  Check local-as-arterial setting for this source
C                       1 - Model both rural and urban local roads as local roads
C                       2 - Model both rural and urban local roads as arterial roads
C                       3 - Model rural local roads as arterial, urban local roads as local
C                       4 - Model rural local roads as local, urban local roads as arterial
                SELECT CASE( SCENLIST( ISRC,2 ) )
                CASE( 1 )
                    RLASAFLAG = .FALSE.
                    ULASAFLAG = .FALSE.
                CASE( 2 )
                    RLASAFLAG = .TRUE.
                    ULASAFLAG = .TRUE.
                CASE( 3 )
                    RLASAFLAG = .TRUE.
                    ULASAFLAG = .FALSE.
                CASE( 4 )
                    RLASAFLAG = .FALSE.
                    ULASAFLAG = .TRUE.
                END SELECT

C.................  Set road and vehicle type
                IF( FTYPE /= M6NONE ) THEN
                    FTYPE = CVTRDTYPE( IRCLAS( ISRC ), RLASAFLAG,
     &                                 ULASAFLAG )
                END IF
                
                VTYPE = CVTVEHTYPE( IVTYPE( ISRC ) )

C.................  Check if vehicle type is valid for this process
                IF( IPOL <= MXM6POLS ) THEN
                    IF( M6VEH2EF( IEF, VTYPE ) == -1 ) CYCLE
                END IF

C.................  Store appropriate emission factor in source-based array
                IF( IPOL <= MXM6POLS ) THEN

                    SRCEFS( EFPOS ) = 
     &                 EMISSIONS( IEF )%PTR( SCENNUM,
     &                                       M6POL2EF( IEF,IPOL ), 
     &                                       M6VEH2EF( IEF,VTYPE ),
     &                                       M6FAC2EF( IEF,FTYPE ),IHR )

C.....................  Adjust for ramp emissions if freeway
                    IF( FTYPE == M6FREEWAY ) THEN
                        RAMPEF = 
     &                     EMISSIONS( IEF )%PTR( SCENNUM,
     &                                           M6POL2EF( IEF,IPOL ), 
     &                                           M6VEH2EF( IEF,VTYPE ),
     &                                           M6FAC2EF( IEF,M6RAMP ),
     &                                           IHR )
                        SRCEFS( EFPOS ) = 
     &                      ADJUST_FREEWAY_EF( SRCEFS( EFPOS ), RAMPEF )
                    END IF
                ELSE
                    SRCEFS( EFPOS ) = 
     &                 HAPEFS( SCENNUM, VTYPE, IPOL - MXM6POLS, 
     &                         IEF, IHR, FTYPE )
     
C.....................  Adjust for ramp emissions if freeway
                    IF( FTYPE == M6FREEWAY ) THEN
                        RAMPEF = 
     &                      HAPEFS( SCENNUM, VTYPE, IPOL - MXM6POLS, 
     &                              IEF, IHR, M6RAMP )
                        SRCEFS( EFPOS ) =
     &                      ADJUST_FREEWAY_EF( SRCEFS( EFPOS ), RAMPEF )
                    END IF
                END IF
                
C.................  Create NONHAP hydrocarbon value if needed
                IF( IPOL == 1 ) THEN
                
C.....................  Loop through subtraction pollutants                
                    DO I = 1, NSUBPOL
                        USEHAP = .FALSE.
                        VARPOL = 0
                        
                        DO J = 1,MXM6POLS
                            IF( SUBPOLS( I ) == M6POLS( J ) ) THEN
                                VARPOL = J
                                EXIT
                            END IF
                        END DO
                        
                        IF( VARPOL == 0 ) THEN
C.............................  Loop through user-defined toxics
                            DO J = 1, NTOTHAPS
                                IF( SUBPOLS( I ) == HAPNAMES( J ) ) THEN
                                    VARPOL = J
                                    USEHAP = .TRUE.
                                    EXIT
                                END IF
                            END DO
                        END IF

C.........................  Still couldn't match pollutant name, quit with error                            
                        IF( VARPOL == 0 ) THEN
                            MESG = 'Unrecognized pollutant ' //
     &                                 SUBPOLS( I )
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                        END IF

C.........................  Get toxic emission factor                                
                        IF( USEHAP ) THEN
                            EFVAL = HAPEFS( SCENNUM, VTYPE, VARPOL, 
     &                                      IEF, IHR, FTYPE )

C.............................  Adjust for ramp emissions
                            IF( FTYPE == M6FREEWAY ) THEN
                                RAMPEF = HAPEFS( SCENNUM, VTYPE, VARPOL,
     &                                           IEF, IHR, M6RAMP )
                                EFVAL = 
     &                              ADJUST_FREEWAY_EF( EFVAL, RAMPEF )
                            END IF
                        
                        ELSE
                        
C.............................  Check that this is a valid pol/ef combo
                            IF( M6POL2EF( IEF, VARPOL ) == -1 ) CYCLE

                            EFVAL = EMISSIONS( IEF )%PTR( SCENNUM,
     &                                       M6POL2EF( IEF,VARPOL ), 
     &                                       M6VEH2EF( IEF,VTYPE ),
     &                                       M6FAC2EF( IEF,FTYPE ),IHR )
     
C.............................  Adjust for ramp emissions
                            IF( FTYPE == M6FREEWAY ) THEN
                                RAMPEF = EMISSIONS( IEF )%PTR( SCENNUM,
     &                                       M6POL2EF( IEF,VARPOL ), 
     &                                       M6VEH2EF( IEF,VTYPE ),
     &                                       M6FAC2EF( IEF,M6RAMP ),
     &                                       IHR )
                                EFVAL = 
     &                              ADJUST_FREEWAY_EF( EFVAL, RAMPEF )
                            END IF
                            
                        END IF

C.........................  Subtract toxic ef from hydrocarbon ef                        
                        SRCEFS( EFPOS ) = SRCEFS( EFPOS ) - EFVAL

C.........................  Make sure the value is not negative
                        IF( SRCEFS( EFPOS ) < 0. ) THEN
                            SRCEFS( EFPOS ) = 0.
                        END IF
                        
                    END DO  ! end subtract pollutant loop
                END IF
              END DO  ! end source loop
                     
C...............  Set current pollutant name based on index
              IF( IPOL <= MXM6POLS ) THEN
                  POLLNAME = M6POLS( IPOL )
                  IF( POLLNAME == 'HC' ) THEN
                      IF( OUTPUTHC /= ' ' ) THEN
                          POLLNAME = OUTPUTHC
                      ELSE
                          POLLNAME = TRIM( VOLNAM )
                      END IF
                  END IF
              ELSE
C...................  Check for user-defined pollutants
                  IF( IPOL <= MXM6POLS + NTOTHAPS ) THEN
                      POLLNAME = HAPNAMES( IPOL - MXM6POLS )
                  ELSE
                      MESG = 'INTERNAL ERROR: Unexpected pollutant'
                      CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                  END IF
              END IF

C...............  Set emission process name based on index
              EFNAME = M6PROCS( IEF )

C...............  Build emission process/pollutant combination name
              CURRVAR = EFNAME // ETJOIN // TRIM( POLLNAME )
              
C...............  Make sure combo is in the MEPROC file
              K = INDEX1( TRIM( CURRVAR ), MXETYPE, EMTNAM( :,1 ) )
              IF( K <= 0 ) CYCLE
              
C...............  Determine write time
              JDATE = SDATE
              JTIME = ( IHR - 1 ) * 10000
              
C...............  Write variable to file    
              IF( .NOT. WRITESET( FNAME, TRIM( CURRVAR ), ALLFILES, 
     &                            JDATE, JTIME, 
     &                            SRCEFS( 1:EFPOS ) ) ) THEN
     	          MESG = 'Could not write ' // TRIM( CURRVAR ) // 
     &	                 'to "' // FNAME( 1:LEN_TRIM( FNAME ) ) // '".'
                  CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
              END IF
                              
            END DO  ! end pollutant loop
          END DO  ! end emission process loop
 
C.............  Write source numbers to file
          IF( .NOT. WRITESET( FNAME, 'SOURCES', ALLFILES, JDATE, 
     &                        JTIME, EFIDX( 1:EFPOS ) ) ) THEN
     	      MESG = 'Could not write source numbers ' //
     &               ' to "' // FNAME( 1:LEN_TRIM( FNAME ) ) // '".'
              CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
          END IF

        END DO     ! end hour loop

C.........  Reset array for storing user-defined HAPs
        IF( ALLOCATED( HAPEFS ) ) HAPEFS = 0.  ! array
        
C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS
        
C.............  This internal subprogram calculates a composite emission
C               factor, based on the freeway and ramp emission factors.
            REAL FUNCTION ADJUST_FREEWAY_EF( FREEWAY_FAC, RAMP_FAC )
            
C.............  PRE: FREEWAY_FAC and RAMP_FAC are assigned
C.............  POST: FNCVAL == (1-Ramp VMT)*FREEWAY_FAC + (Ramp VMT)*RAMP_FAC

C.............  Function arguments
            REAL, INTENT (IN) :: FREEWAY_FAC    ! freeway emission factor
            REAL, INTENT (IN) :: RAMP_FAC       ! ramp emission factor
            
            ADJUST_FREEWAY_EF = ( 1-RAMPVMT ) * FREEWAY_FAC + 
     &                          ( RAMPVMT )   * RAMP_FAC
            
            END FUNCTION ADJUST_FREEWAY_EF
        
        END SUBROUTINE WREMFACS
        
