
        SUBROUTINE SETNONHAP

C**************************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine is integrates criteria and toxics pollutant 
C      emissions by creating NONHAP values.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 11/02 by C. Seppanen
C
C**************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: POLVAL, IPOSCOD, NPCNT, CSOURC

C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: INVSTAT, MXIDAT, INVDNAM, INVDVTS,
     &                      ITMSPC, ITEXPL

C...........   This module contains the cross-reference tables
        USE MODXREF, ONLY: LNONHAP
        
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC, NPPOL, NCHARS, NEM, NOZ

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        
C...........   EXTERNAL FUNCTIONS and their descriptions
        CHARACTER*2     CRLF
        INTEGER         INDEX1
        INTEGER         ENVINT
        
        EXTERNAL        CRLF, INDEX1, ENVINT
        
C.........  Pollutant names
        CHARACTER(LEN=IOVLEN3),PARAMETER :: VOCNAM = 'VOC'
        CHARACTER(LEN=IOVLEN3),PARAMETER :: TOGNAM = 'TOG'
        CHARACTER(LEN=IOVLEN3),PARAMETER :: NHVNAM = 'NONHAPVOC'
        CHARACTER(LEN=IOVLEN3),PARAMETER :: NHTNAM = 'NONHAPTOG'
        CHARACTER(LEN=4      ),PARAMETER :: NOIEND = '_NOI'

C.........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: TMPIDX( : )      ! sorting index
        INTEGER, ALLOCATABLE :: TMPPOSCOD( : )   ! tmp pollutant code array
        
        REAL,    ALLOCATABLE :: TMPPOLVAL( :,: ) ! tmp emissions array
        
C.........   Other local variables
        INTEGER  I,J,K,S,L2   ! counters and indices
        
        INTEGER  IDXSIZE      ! size of sorting index
        INTEGER  IOS          ! I/O status
        INTEGER  CPOL         ! current pollutant number
        INTEGER  PPOL         ! previous pollutant number
        INTEGER  CURRPOS      ! current position in POLVAL and IPOSCOD arrays
        INTEGER  VNMPOS       ! position of VOC in pollutant names array
        INTEGER  TNMPOS       ! position of TOG in pollutant names array
        INTEGER  NONVPOS      ! position of NONHAPVOC in poll name array
        INTEGER  NONTPOS      ! position of NONHAPTOG in poll name array
        INTEGER  VOCPOS       ! location of VOC entry in srcs array
        INTEGER  TOGPOS       ! location of TOG entry in srcs array
        INTEGER  TOXPOS       ! position of NOI toxic in poll name array
        INTEGER  MXWARN       ! maximum number of warnings
        INTEGER :: NWARN = 0   ! current number of warnings
        INTEGER :: NCRNOTOX = 0! number of sources with criteria but no toxics
        INTEGER :: NTOXNOCR = 0! number of sources with toxics but no criteria
        INTEGER  STIDX        ! starting index into pollutant array for current source
        INTEGER  ENDIDX       ! ending index into pollutant array
        
        REAL     VOCEANN      ! summed annual VOC emissions
        REAL     VOCEOZN      ! summed ozone season VOC emissions
        REAL     TOGEANN      ! summed annual TOG emissions
        REAL     TOGEOZN      ! summed ozone season TOG emissions
        
        LOGICAL  FNDPOL       ! true: found toxic pollutant to be processed
        LOGICAL  FNDVOC       ! true: found VOC entry in inventory
        LOGICAL  FNDTOG       ! true: found TOG entry in inventory
        LOGICAL  EFLAG        ! true: error occured
        LOGICAL  LASTFLAG     ! true: entry is last for current source
        LOGICAL  NEEDSORT     ! true: need to resort pollutants for current source
        
        CHARACTER(LEN=IOVLEN3) POLNAM   ! temporary pollutant name
        CHARACTER(LEN=100    ) BUFFER   ! message buffer
        CHARACTER(LEN=256    ) MESG     ! message buffer 
        
        CHARACTER*16 :: PROGNAME = 'SETNONHAP' ! program name

C***********************************************************************
C   begin body of subroutine SETNONHAP

C.........  Get maximum number of warnings
        MXWARN = ENVINT( WARNSET, ' ', 100, IOS )
        NWARN = 0

C.........  Find position of pollutant names in pollutant array
        VNMPOS  = INDEX1( VOCNAM, MXIDAT, INVDNAM )
        TNMPOS  = INDEX1( TOGNAM, MXIDAT, INVDNAM )
        NONVPOS = INDEX1( NHVNAM, MXIDAT, INVDNAM )
        NONTPOS = INDEX1( NHTNAM, MXIDAT, INVDNAM )
        
C.........  Check that either VOC or TOG is in the inventory
        IF( ( VNMPOS == 0 .OR. INVSTAT( VNMPOS ) == 0 ) .AND.
     &      ( TNMPOS == 0 .OR. INVSTAT( TNMPOS ) == 0 )       ) RETURN
        
C.........  Check if any toxics pollutants are to be subtracted
        FNDPOL = .FALSE.
        DO I = 1, MXIDAT
            IF( INVDVTS( I ) /= 'N' ) THEN
                IF( INVSTAT( I ) /= 0 ) THEN
                    FNDPOL = .TRUE.
                    EXIT
                END IF
            END IF
        END DO

        IF( .NOT. FNDPOL ) RETURN

        CURRPOS = 0

C.........  Loop through sources
        DO I = 1, NSRC

C.............  Initialize values for this source
            VOCPOS = 0
            TOGPOS = 0
            FNDVOC = .FALSE.
            FNDTOG = .FALSE.
            
            VOCEANN = 0.
            VOCEOZN = 0.
            TOGEANN = 0.
            TOGEOZN = 0.
            
            EFLAG = .FALSE.
            LASTFLAG = .FALSE.

C.............  Process source if it is not integrated
            IF( ALLOCATED( LNONHAP ) .AND. .NOT. LNONHAP( I ) ) THEN

C.................  Loop through pollutants for this source
                DO J = 1, NPCNT( I )

C.....................  Increment current position in arrays            
                    CURRPOS = CURRPOS + 1
                    
C.....................  Store pollutant for this position
                    CPOL = IPOSCOD( CURRPOS )
                    
C.....................  If pollutant is not part of VOC or TOG, cycle
                    IF( INVDVTS( CPOL ) == 'N' ) CYCLE

C.....................  If pollutant is not a model species, set it to zero
                    IF( .NOT. ITMSPC( CPOL ) ) THEN
                        POLVAL( CURRPOS,NEM ) = 0.0
                        POLVAL( CURRPOS,NOZ ) = 0.0

C..................... Otherwise, if pollutant is not an explicit species, rename to NOI
                    ELSE IF( .NOT. ITEXPL( CPOL ) ) THEN
                
C.........................  Create NOI name
                        POLNAM = INVDNAM( CPOL )
                        IF( LEN_TRIM( POLNAM ) > 11 ) THEN
                            POLNAM = POLNAM(1:11)
                        END IF
                        POLNAM = TRIM( POLNAM ) // NOIEND

C.........................  Find NOI name in pollutant names array                        
                        TOXPOS = INDEX1( POLNAM, MXIDAT, INVDNAM )

C.........................  If found, set pollutant for current source to NOI pollutant
                        IF( TOXPOS > 0 ) THEN
                            IPOSCOD( CURRPOS ) = TOXPOS
                            INVSTAT( TOXPOS ) = 2
                        ELSE
                            MESG = 'ERROR: Non-integrated toxic ' //
     &                             'pollutant ' // TRIM( POLNAM ) //
     &                             ' was not found in the ' //
     &                             'INVTABLE file.'
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                        END IF
                     
                    END IF  ! check if pollutant is model or explicit species
                    
                END DO  ! loop through pollutants
                
C.................  Skip rest of loop since we're done with this source
                CYCLE

            END IF

C.............  Loop through pollutants for this source            
            DO J = 1, NPCNT( I )

C.................  Increment current position in arrays            
                CURRPOS = CURRPOS + 1

C.................  Store pollutant for this position
                CPOL = IPOSCOD( CURRPOS )

C.................  If current pollutant is VOC or TOG entry, save position
                IF( CPOL == VNMPOS ) THEN
                    VOCPOS = CURRPOS
                END IF
                
                IF( CPOL == TNMPOS ) THEN
                    TOGPOS = CURRPOS
                END IF
            
C.................  Check if this is last pollutant for this source
                IF( J == NPCNT( I ) ) THEN
                    LASTFLAG = .TRUE.

C.................  Otherwise, if not part of VOC or TOG, skip rest of loop
                ELSE IF( INVDVTS( CPOL ) == 'N' ) THEN
                    CYCLE
                    
                END IF

C.................  Sum toxic emissions for this source
                IF( INVDVTS( CPOL ) == 'V' ) THEN
                    VOCEANN = VOCEANN + POLVAL( CURRPOS,NEM )
                    VOCEOZN = VOCEOZN + POLVAL( CURRPOS,NOZ )
                    FNDVOC = .TRUE.
                END IF
                
                IF( INVDVTS( CPOL ) == 'T' .OR. 
     &              INVDVTS( CPOL ) == 'V' ) THEN
                    TOGEANN = TOGEANN + POLVAL( CURRPOS,NEM )
                    TOGEOZN = TOGEOZN + POLVAL( CURRPOS,NOZ )
                    FNDTOG = .TRUE.
                END IF

C.................  If this is not the last entry for source, cycle
C                   Otherwise, check conditions and subtract toxic
C                   emissions from criteria values
                IF( .NOT. LASTFLAG ) CYCLE

                LASTFLAG = .FALSE.

C.................  Format information for this source
                CALL FMTCSRC( CSOURC( I ), NCHARS, BUFFER, L2 )
                
C.................  Give warning if source has toxics but no criteria
                IF( VOCPOS == 0 .AND. FNDVOC ) THEN
                    IF( NWARN <= MXWARN ) THEN
                        MESG = 'WARNING: Found toxic emissions ' //
     &                         'but no VOC emissions for source:' // 
     &                         CRLF() // BLANK5 // BUFFER( 1:L2 )
                        CALL M3MESG( MESG )
                        NWARN = NWARN + 1
                    END IF
                
                    NTOXNOCR = NTOXNOCR + 1
                    EFLAG = .TRUE.
                END IF
                
                IF( TOGPOS == 0 .AND. FNDTOG ) THEN
                    IF( NWARN <= MXWARN ) THEN
                        MESG = 'WARNING: Found toxic emissions ' //
     &                         'but no TOG emissions for source:' // 
     &                         CRLF() // BLANK5 // BUFFER( 1:L2 )
                        CALL M3MESG( MESG )
                        NWARN = NWARN + 1
                    END IF
                    
                    NTOXNOCR = NTOXNOCR + 1
                    EFLAG = .TRUE.
                END IF
                
C.................  Give warning if source has criteria but no toxics
                IF( VOCPOS /= 0 .AND. .NOT. FNDVOC ) THEN
                    IF( NWARN <= MXWARN ) THEN
                        MESG = 'WARNING: Found VOC emissions ' //
     &                         'but no toxic emissions for source:' // 
     &                         CRLF() // BLANK5 // BUFFER( 1:L2 )
                        CALL M3MESG( MESG )
                        NWARN = NWARN + 1
                    END IF
                    
                    NCRNOTOX = NCRNOTOX + 1
                    EFLAG = .TRUE.
                END IF
                
                IF( TOGPOS /= 0 .AND. .NOT. FNDTOG ) THEN
                    IF( NWARN <= MXWARN ) THEN
                        MESG = 'WARNING: Found TOG emissions ' //
     &                         'but no toxic emissions for source:' //
     &                         CRLF() // BLANK5 // BUFFER( 1:L2 )
                        CALL M3MESG( MESG )
                        NWARN = NWARN + 1
                    END IF
                    
                    NCRNOTOX = NCRNOTOX + 1
                    EFLAG = .TRUE.
                END IF
                
C.................  Skip rest of loop if an error has occured
                IF( EFLAG ) CYCLE
                
C.................  Subtract toxic emissions from criteria emissions  
                IF( VOCPOS /= 0 ) THEN                  
                    POLVAL( VOCPOS,NEM ) = 
     &                  POLVAL( VOCPOS,NEM ) - VOCEANN
                    POLVAL( VOCPOS,NOZ ) =
     &                  POLVAL( VOCPOS,NOZ ) - VOCEOZN

C.....................  Check that annual NONHAP value is not negative
                    IF( POLVAL( VOCPOS,NEM ) < 0. ) THEN
                        IF( NWARN <= MXWARN ) THEN
                            MESG = 'WARNING: Total annual toxic ' //
     &                             'emissions greater than annual ' //
     &                             'VOC emissions for source:' // 
     &                             CRLF() // BLANK5 // BUFFER( 1:L2 )
                            CALL M3MESG( MESG )
                            NWARN = NWARN + 1
                        END IF
                        
                        POLVAL( VOCPOS,NEM ) = .0
                    END IF
     
C.....................  Check that ozone season NONHAP value is not negative
                    IF( POLVAL( VOCPOS,NOZ ) < 0. ) THEN
                        IF( NWARN <= MXWARN ) THEN
                            MESG = 'WARNING: Total ozone season ' //
     &                             'toxic emissions greater than ' //
     &                             'ozone season VOC emissions for ' //
     &                             'source:' // CRLF() // BLANK5 // 
     &                             BUFFER( 1:L2 )
                            CALL M3MESG( MESG )
                            NWARN = NWARN + 1
                        END IF
                        
                        POLVAL( VOCPOS,NOZ ) = .0
                    END IF
                    
C.....................  Rename VOC to NONHAPVOC
                    IPOSCOD( VOCPOS ) = NONVPOS
                    INVSTAT( NONVPOS ) = 2
     
                END IF
                
                IF( TOGPOS /= 0 ) THEN
                    POLVAL( TOGPOS,NEM ) = 
     &                  POLVAL( TOGPOS,NEM ) - TOGEANN
                    POLVAL( TOGPOS,NOZ ) =
     &                  POLVAL( TOGPOS,NOZ ) - TOGEOZN
     
C.....................  Check that annual NONHAP value is not negative
                    IF( POLVAL( TOGPOS,NEM ) < 0. ) THEN
                        IF( NWARN <= MXWARN ) THEN
                            MESG = 'WARNING: Total annual toxic ' //
     &                             'emissions greater than annual ' //
     &                             'TOG emissions for source:' // 
     &                             CRLF() // BLANK5 // BUFFER( 1:L2 )
                            CALL M3MESG( MESG )
                            NWARN = NWARN + 1
                        END IF
                        
                        POLVAL( TOGPOS,NEM ) = .0
                    END IF
                
C.....................  Check that ozone season NONHAP value is not negative
                    IF( POLVAL( TOGPOS,NOZ ) < 0. ) THEN
                        IF( NWARN <= MXWARN ) THEN
                            MESG = 'WARNING: Total ozone season ' //
     &                             'toxic emissions greater than ' //
     &                             'ozone season TOG emissions for ' //
     &                             'source:' // CRLF() // BLANK5 // 
     &                             BUFFER( 1:L2 )
                            CALL M3MESG( MESG )
                            NWARN = NWARN + 1
                        END IF
                        
                        POLVAL( TOGPOS,NOZ ) = .0
                    END IF
                    
C....................  Rename TOG to NONHAPTOG
                   IPOSCOD( TOGPOS ) = NONTPOS
                   INVSTAT( NONTPOS ) = 2
     
                END IF

            END DO  ! loop through pollutants
            
        END DO  ! loop through sources

        IF( NCRNOTOX > 0 .OR. NTOXNOCR > 0 ) THEN
            WRITE( MESG,94010 ) 'During processing, the following ' //
     &         'number of sources were encountered: ' // 
     &         CRLF() // BLANK10 // 
     &         'Sources with VOC or TOG emissions but no toxic ' //
     &         'emissions: ', NCRNOTOX, 
     &         CRLF() // BLANK10 //
     &         'Sources with toxic emissions but no VOC or TOG ' //
     &         'emissions: ', NTOXNOCR
            CALL M3MESG( MESG )
        END IF

C.........  Sort POLVAL and IPOSCOD to put new pollutants (NONHAP and NOI) in correct order

C.........  Determine maximum size for sorting index, then allocate memory
        IDXSIZE = MAXVAL( NPCNT )
        
        ALLOCATE( TMPIDX( IDXSIZE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPIDX', PROGNAME )
        ALLOCATE( TMPPOLVAL( IDXSIZE,NPPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPPOLVAL', PROGNAME )
        ALLOCATE( TMPPOSCOD( IDXSIZE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPPOSCOD', PROGNAME )

        STIDX = 1
        
C.........  Loop through sources        
        DO I = 1, NSRC

C.............  Set ending index for this source
            ENDIDX = STIDX + NPCNT( I ) - 1

C.............  Initialize sorting index and check if any sorting needs to be done
            PPOL = 0
            NEEDSORT = .FALSE.
            
            DO J = 1, NPCNT( I )
                CPOL = IPOSCOD( STIDX + J - 1 )

C.................  If current pollutant is lower than previous, need to resort                
                IF( CPOL < PPOL ) THEN
                    NEEDSORT = .TRUE.
                END IF
                
                PPOL = CPOL
                TMPIDX( J ) = J
            END DO

C.............  Make sure this source needs to be resorted
            IF( NEEDSORT ) THEN

C.................  Store current values in temporary arrays
                TMPPOLVAL( 1:NPCNT( I ),: ) = POLVAL ( STIDX:ENDIDX,: )
                TMPPOSCOD( 1:NPCNT( I ) )   = IPOSCOD( STIDX:ENDIDX )
                
C.................  Sort section of IPOSCOD corresponding to this source
                CALL SORTI1( NPCNT( I ), TMPIDX( 1:NPCNT( I ) ), 
     &                       IPOSCOD( STIDX:ENDIDX ) ) 

C.................  Loop through pollutants for current source        
                DO J = 1, NPCNT( I )
            
                    K = TMPIDX( J )
                    
                    POLVAL ( STIDX + J - 1,: ) = TMPPOLVAL( K,: )
                    IPOSCOD( STIDX + J - 1 )   = TMPPOSCOD( K )
            
                END DO

            END IF  ! check if pollutants need sorting

C.............  Increment starting index        
            STIDX = ENDIDX + 1
        
        END DO

C.........  Deallocate local memory
        DEALLOCATE( TMPIDX, TMPPOLVAL, TMPPOSCOD )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE SETNONHAP
       
