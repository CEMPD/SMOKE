
        SUBROUTINE PROCAR2PT( NRAWBP )

C**************************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine process the area-to-point sources, adding X and Y
C      locations and adjusting the annual and ozone-season emissions.
C
C  PRECONDITIONS REQUIRED:
C      
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Split from adjustinv.f 1/03 by C. Seppanen 
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
        USE MODSOURC, ONLY: CSOURC, NPCNT, IPOSCOD, POLVAL, 
     &                      XLOCA, YLOCA

C...........   This module contains the cross-reference tables
        USE MODXREF, ONLY: AR2PTIDX, AR2PTTBL, AR2PTCNT
        
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC, NEM, NOZ, NEF, NRP, NPPOL

C.........  This module contains the arrays for the area-to-point x-form
        USE MODAR2PT, ONLY: AR2PTABL
        
        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions
!        CHARACTER*2     CRLF
!        INTEGER         STR2INT
        
!        EXTERNAL        CRLF, STR2INT

C...........   SUBROUTINE ARGUMENTS
        INTEGER , INTENT (INOUT) :: NRAWBP  ! no. raw records by pollutant

C...........   Local pointers
        INTEGER, POINTER :: OLDNPCNT ( : ) !  number of pollutants per source
        INTEGER, POINTER :: OLDIPOSCOD( : ) !  positn of pol in INVPCOD

        REAL   , POINTER :: OLDPOLVAL( :,: )   ! emission values

        CHARACTER(LEN=ALLLEN3), POINTER :: OLDCSOURC( : ) ! concat src

C...........   Other local variables
        INTEGER         I,J,K,S     ! counters
        INTEGER         IOS         ! I/O error status
        INTEGER         NA2PSRCS    ! no. of area-to-point sources to add
        INTEGER         NA2PRECS    ! no. of sources with pollutants to add
        INTEGER         NEWSRCPOS   ! position in new source arrays
        INTEGER         NEWRECPOS   ! position in new source w/ pollutant arrays
        INTEGER         OLDRECPOS   ! position in old source w/ pollutant arrays
        INTEGER         NREPSRCS    ! total no. of area-to-point sources
        INTEGER         TBLE        ! current area-to-point table number
        INTEGER         ROW         ! current area-to-point row
        INTEGER         OLDNSRC     ! old number of srcs

        REAL            FACTOR      ! factor for area-to-point conversion

        CHARACTER(LEN=256    )  MESG        !  message buffer 

        CHARACTER*16 :: PROGNAME = 'PROCAR2PT' ! program name

C***********************************************************************
C   begin body of subroutine PROCAR2PT
            
C.........  Determine total number of sources to be added; if source only
C           has one location, can use existing arrays and don't need to 
C           create additional memory for it
        NA2PSRCS = 0
        NA2PRECS = 0
        NREPSRCS = 0

        DO S = 1, NSRC
            IF( AR2PTTBL( S ) /= 0 ) THEN
                NA2PSRCS = NA2PSRCS + AR2PTCNT( S ) - 1
                NA2PRECS = NA2PRECS + (AR2PTCNT( S ) - 1)*NPCNT( S )
                NREPSRCS = NREPSRCS + 1
            END IF
        END DO

C.........  Check if any sources need to be processed
        IF( NREPSRCS == 0 ) RETURN

C.........  Store old number of sources
        OLDNSRC = NSRC

        IF( NA2PSRCS > 0 ) THEN

C.............  Update total number of sources and records
            NSRC   = NSRC   + NA2PSRCS
            NRAWBP = NRAWBP + NA2PRECS

C.............  Associate temporary pointers with sorted arrays
            OLDCSOURC  => CSOURC
            OLDNPCNT   => NPCNT
            OLDPOLVAL  => POLVAL
            OLDIPOSCOD => IPOSCOD

C.............  Nullify original sorted arrays
            NULLIFY( CSOURC, POLVAL, NPCNT, IPOSCOD )

C.............  Allocate memory for larger sorted arrays
            ALLOCATE( CSOURC( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CSOURC', PROGNAME )
            ALLOCATE( NPCNT( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NPCNT', PROGNAME )
            ALLOCATE( POLVAL( NRAWBP,NPPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'POLVAL', PROGNAME )
            ALLOCATE( IPOSCOD( NRAWBP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'IPOSCOD', PROGNAME )
            NPCNT = 0         ! array
            POLVAL = BADVAL3  ! array
            
        END IF

C.........  Allocate memory for X and Y locations (in sorted order)
        ALLOCATE( XLOCA( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'XLOCA', PROGNAME )
        ALLOCATE( YLOCA( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'YLOCA', PROGNAME )
        XLOCA = BADVAL3  ! array
        YLOCA = BADVAL3  ! array
        
        NEWSRCPOS = 0
        NEWRECPOS = 0
        OLDRECPOS = 1

C.........  Loop through original sources
        DO S = 1, OLDNSRC

C.............  Check if current source is to be processed
            IF( AR2PTTBL( S ) /= 0 ) THEN

C.................  Set area-to-point table and row for current source                	
                TBLE = AR2PTTBL( S )
                ROW  = AR2PTIDX( S )

C.................  If no sources are actually added, just modify existing values
                IF( NA2PSRCS == 0 ) THEN
                    XLOCA( S ) = AR2PTABL( ROW,TBLE )%LON
                    YLOCA( S ) = AR2PTABL( ROW,TBLE )%LAT
                    
                    FACTOR = AR2PTABL( ROW,TBLE )%ALLOC
                    
C.....................  If not adding sources, factor should be 1.0, so skip 
C                       adjusting emissions
                    IF( FACTOR /= 1. ) THEN                    
                    
                        DO K = OLDRECPOS, OLDRECPOS + NPCNT( S ) - 1
                            
                            IF( POLVAL( K,NEM ) /= BADVAL3 ) THEN
                                POLVAL( K,NEM ) = POLVAL( K,NEM )*FACTOR
                            END IF
                            
                            IF( POLVAL( K,NOZ ) /= BADVAL3 ) THEN
                                POLVAL( K,NOZ ) = POLVAL( K,NOZ )*FACTOR
                            END IF
                            
                        END DO
                        
                    END IF
                    
                    OLDRECPOS = OLDRECPOS + NPCNT( S )
                
                ELSE

C.....................  Loop through all locations for this source
                    DO J = 0, AR2PTCNT( S ) - 1

C.........................  Increment source position and copy source info
                        NEWSRCPOS = NEWSRCPOS + 1
                        CSOURC( NEWSRCPOS ) = OLDCSOURC( S )
                        NPCNT( NEWSRCPOS ) = OLDNPCNT( S )

C.........................  Store X and Y locations
                        XLOCA( NEWSRCPOS ) = AR2PTABL( ROW+J,TBLE )%LON
                        YLOCA( NEWSRCPOS ) = AR2PTABL( ROW+J,TBLE )%LAT

C.........................  Set factor for this location
                        FACTOR = AR2PTABL( ROW+J,TBLE )%ALLOC

C.........................  Loop through pollutants for this source
                        DO K = OLDRECPOS, OLDRECPOS + OLDNPCNT( S ) - 1
                    
C.............................  Increment record position and store pollutant code
                            NEWRECPOS = NEWRECPOS + 1
                            IPOSCOD( NEWRECPOS ) = OLDIPOSCOD( K )

C.............................  Adjust annual and ozone season emissions based on ar2pt factor
                            IF( OLDPOLVAL( K,NEM ) /= BADVAL3 ) THEN
                                POLVAL( NEWRECPOS, NEM ) = 
     &                              OLDPOLVAL( K,NEM ) * FACTOR
                            END IF
                
                            IF( OLDPOLVAL( K,NOZ ) /= BADVAL3 ) THEN
                                POLVAL( NEWRECPOS,NOZ ) = 
     &                              OLDPOLVAL( K,NOZ ) * FACTOR
                            END IF

C.............................  Copy remaining values to new array
                            POLVAL( NEWRECPOS,NEF:NRP ) = 
     &                          OLDPOLVAL( K,NEF:NRP )
                        
                        END DO  ! loop through pollutants
                    
                    END DO  ! loop through area-to-point locations
                    
                    OLDRECPOS = OLDRECPOS + OLDNPCNT( S )

                END IF  ! end check if any sources are added
            
            ELSE

C.................  Not processing current source, but if we're adding sources,
C                   then need to copy information to new arrays
                IF( NA2PSRCS > 0 ) THEN
                
                    NEWSRCPOS = NEWSRCPOS + 1
                    CSOURC( NEWSRCPOS ) = OLDCSOURC( S )
                    NPCNT( NEWSRCPOS ) = OLDNPCNT( S )
                    
                    DO K = OLDRECPOS, OLDRECPOS + OLDNPCNT( S ) - 1
                        
                        NEWRECPOS = NEWRECPOS + 1
                        IPOSCOD( NEWRECPOS ) = OLDIPOSCOD( K )
                        POLVAL( NEWRECPOS,: ) = OLDPOLVAL( K,: )
                        
                    END DO
                
                    OLDRECPOS = OLDRECPOS + OLDNPCNT( S )
                
                END IF
                
            END IF  ! check if source is processed
            
        END DO  ! loop through sources

C.........  Deallocate old arrays
        IF( NA2PSRCS > 0 ) THEN
            DEALLOCATE( OLDCSOURC, OLDNPCNT, OLDPOLVAL, OLDIPOSCOD )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
 
        END SUBROUTINE PROCAR2PT
