
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
        USE MODSOURC 

C.........  This module contains the lists of unique inventory information
        USE MODLISTS

C...........   This module contains the cross-reference tables
        USE MODXREF
        
C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        
C...........   EXTERNAL FUNCTIONS and their descriptions
        INTEGER         INDEX1
        
        EXTERNAL        INDEX1
        
C.........  Pollutant names
        CHARACTER(LEN=IOVLEN3),PARAMETER :: VOCNAM = 'VOC'
        CHARACTER(LEN=IOVLEN3),PARAMETER :: TOGNAM = 'TOG'
        CHARACTER(LEN=IOVLEN3),PARAMETER :: NHVNAM = 'NONHAPVOC'
        CHARACTER(LEN=IOVLEN3),PARAMETER :: NHTNAM = 'NONHAPTOG'
        CHARACTER(LEN=4      ),PARAMETER :: NOIEND = '_NOI'
        
C.........   Other local variables
        INTEGER  I,J,S        ! counters and indices
        INTEGER  POL          ! pollutant number
        INTEGER  CURRPOS      ! current position in POLVAL and IPOSCOD arrays
        INTEGER  VNMPOS       ! position of VOC in pollutant names array
        INTEGER  TNMPOS       ! position of TOG in pollutant names array
        INTEGER  NONVPOS      ! position of NONHAPVOC in poll name array
        INTEGER  NONTPOS      ! position of NONHAPTOG in poll name array
        INTEGER  VOCPOS       ! location of VOC entry in srcs array
        INTEGER  TOGPOS       ! location of TOG entry in srcs array
        INTEGER  TOXPOS       ! position of NOI toxic in poll name array
        
        REAL     VOCEANN      ! summed annual VOC emissions
        REAL     TOGEANN      ! summed annual TOG emissions
        
        LOGICAL  FNDPOL       ! true: found toxic pollutant to be processed
        LOGICAL  FNDVOC       ! true: found VOC entry in inventory
        LOGICAL  FNDTOG       ! true: found TOG entry in inventory
        LOGICAL  EFLAG        ! true: error occured
        LOGICAL  LASTFLAG     ! true: entry is last for current source
        
        CHARACTER(LEN=IOVLEN3) POLNAM  ! temporary pollutant name
        CHARACTER(LEN=256    ) MESG        ! message buffer 
        
        CHARACTER*16 :: PROGNAME = 'SETNONHAP' ! program name

C***********************************************************************
C   begin body of subroutine SETNONHAP

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

C.........  Initialize values for loop
        VOCPOS = 0
        TOGPOS = 0
        FNDVOC = .FALSE.
        FNDTOG = .FALSE.
                    
        VOCEANN = 0.
        TOGEANN = 0.
                    
        EFLAG = .FALSE.
        LASTFLAG = .FALSE.

        CURRPOS = 0

C.........  Loop through sources
        DO I = 1, NSRC

C.............  Loop through pollutants for this source            
            DO J = 1, NPCNT( I )

C.................  Increment current position in arrays            
                CURRPOS = CURRPOS + 1

C.................  Store pollutant for this position
                POL = IPOSCOD( CURRPOS )
            
C.................  Check if this is last pollutant for this source
                IF( J == NPCNT( I ) ) THEN
                    LASTFLAG = .TRUE.
                END IF

C.................  Process source if it is not integrated
                IF( ALLOCATED( LNONHAP ) .AND. .NOT. LNONHAP( I ) ) THEN
            
C.....................  Skip rest of loop if pollutant is not part VOC or TOG
                    IF( INVDVTS( POL ) == 'N' ) CYCLE

C.....................  If pollutant is not a model species, set it to zero
                    IF( .NOT. ITMSPC( POL ) ) THEN
                        POLVAL( CURRPOS,NEM ) = 0.0

                    ELSE

C......................... If pollutant is not an explicit species, rename to NOI
                        IF( .NOT. ITEXPL( POL ) ) THEN
                    
C.............................  Create NOI name
                            POLNAM = INVDNAM( POL )
                            IF( LEN_TRIM( POLNAM ) > 11 ) THEN
                                POLNAM = POLNAM(1:11)
                            END IF
                            POLNAM = TRIM( POLNAM ) // NOIEND

C.............................  Find NOI name in pollutant names array                        
                            TOXPOS = INDEX1( POLNAM, MXIDAT, INVDNAM )

C.............................  If found, set pollutant for current source to NOI pollutant
                            IF( TOXPOS > 0 ) THEN
                                IPOSCOD( CURRPOS ) = TOXPOS
                                INVSTAT( TOXPOS ) = 2
                            ELSE
                                MESG = 'ERROR: Non-integrated toxic ' //
     &                                 'pollutant ' // TRIM( POLNAM ) //
     &                                 ' was not found in the ' //
     &                                 'INVTABLE file.'
                                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                            END IF
                         
                        END IF  ! check if pollutant is an explicit species
                    
                    END IF  ! check if pollutant is a model species
                            
C.................  Otherwise current source is to be integrated           
                ELSE

C.....................  If current pollutant is VOC or TOG entry, save position
                    IF( POL == VNMPOS ) THEN
                        VOCPOS = CURRPOS
                        IF( .NOT. LASTFLAG ) CYCLE
                    END IF
                
                    IF( POL == TNMPOS ) THEN
                        TOGPOS = CURRPOS
                        IF( .NOT. LASTFLAG ) CYCLE
                    END IF

C.....................  Sum toxic emissions for this source
                    IF( INVDVTS( POL ) == 'V' ) THEN
                        VOCEANN = VOCEANN + POLVAL( CURRPOS,NEM )
                        FNDVOC = .TRUE.
                    END IF
                
                    IF( INVDVTS( POL ) == 'T' ) THEN
                        TOGEANN = TOGEANN + POLVAL( CURRPOS,NEM )
                        FNDTOG = .TRUE.
                    END IF

C.....................  If this is last entry for source, check conditions and
C                       subtract toxic emissions from criteria values
                    IF( LASTFLAG ) THEN
                
C.........................  Give warning if source has toxics but no criteria
                        IF( VOCPOS == 0 .AND. FNDVOC ) THEN
                            MESG = 'WARNING: Source has VOC toxic ' //
     &                             'but no criteria emissions'
                            CALL M3MESG( MESG )
                            EFLAG = .TRUE.
                        END IF
                    
                        IF( TOGPOS == 0 .AND. FNDTOG ) THEN
                            MESG = 'WARNING: Source has TOG toxic ' //
     &                             'but no criteria emissions'
                            CALL M3MESG( MESG )
                            EFLAG = .TRUE.
                        END IF
                    
C.........................  Give warning if source has criteria but no toxics
                        IF( VOCPOS /= 0 .AND. .NOT. FNDVOC ) THEN
                            MESG = 'WARNING: Source has VOC criteria' //
     &                             ' but no toxic emissions'
                            CALL M3MESG( MESG )
                            EFLAG = .TRUE.
                        END IF
                    
                        IF( TOGPOS == 0 .AND. FNDTOG ) THEN
                            MESG = 'WARNING: Source has TOG criteria' //
     &                             ' but no toxic emissions'
                            CALL M3MESG( MESG )
                            EFLAG = .TRUE.
                        END IF

C.........................  Make sure an error has not occured
                        IF( .NOT. EFLAG ) THEN
                    
C.............................  Subtract toxic emissions from criteria emissions  
                            IF( VOCPOS /= 0 ) THEN                  
                                POLVAL( VOCPOS,NEM ) = 
     &                              POLVAL( VOCPOS,NEM ) - VOCEANN

C.................................  Check that NONHAP value is not negative
                                IF( POLVAL( VOCPOS,NEM ) < 0. ) THEN
                                    MESG = 'WARNING: Toxic emissions' //
     &                                     ' greater than criteria'
                                    CALL M3MESG( MESG )
                                    POLVAL( VOCPOS,NEM ) = .0
                                END IF
     
C.................................  Rename VOC to NONHAPVOC
                                IPOSCOD( VOCPOS ) = NONVPOS
                                INVSTAT( NONVPOS ) = 2
     
                            END IF
                    
                            IF( TOGPOS /= 0 ) THEN
                                POLVAL( TOGPOS,NEM ) = 
     &                              POLVAL( TOGPOS,NEM ) - TOGEANN
     
C.................................  Check that NONHAP value is not negative
                                IF( POLVAL( TOGPOS,NEM ) < 0. ) THEN
                                    MESG = 'WARNING: Toxic emissions' //
     &                                     ' greater than criteria'
                                    CALL M3MESG( MESG )
                                    POLVAL( TOGPOS,NEM ) = .0
                                END IF
                            
C................................  Rename TOG to NONHAPTOG
                               IPOSCOD( TOGPOS ) = NONTPOS
                               INVSTAT( NONTPOS ) = 2
     
                            END IF
                        
                        END IF  ! check for error
     
C.........................  Reset flags and values
                        VOCPOS = 0
                        TOGPOS = 0
                        FNDVOC = .FALSE.
                        FNDTOG = .FALSE.
                    
                        VOCEANN = 0.
                        TOGEANN = 0.
                    
                        EFLAG = .FALSE.
                        LASTFLAG = .FALSE.
                    END IF  ! check if entry is last for source

                END IF  ! check if source is integrated or not
        
            END DO  ! loop through pollutants
            
        END DO  ! loop through sources

        END SUBROUTINE SETNONHAP
       