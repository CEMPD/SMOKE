
        SUBROUTINE ADJUSTINV( NRAWBP, UDEV, YDEV, CDEV, LDEV )

C**************************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine is just a placeholder to call the area-to-point
C      routines for now.
C
C  PRECONDITIONS REQUIRED:
C      CSOURCA and SRCIDA arrays allocated and defined with Source IDs
C      NSRC is set.
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

C...........   This module contains the cross-reference tables
        USE MODXREF
        
C.........  This module contains the information about the source category
        USE MODINFO

C.........  This module contains the arrays for the area-to-point x-form
        USE MODAR2PT
        
        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions
        CHARACTER*2     CRLF
        INTEGER         STR2INT
        
        EXTERNAL        CRLF, STR2INT

C...........   SUBROUTINE ARGUMENTS
        INTEGER , INTENT (INOUT) :: NRAWBP  ! no. raw records by pollutant
        INTEGER , INTENT (IN)    :: UDEV    ! unit no. for non-HAP exclusions
        INTEGER , INTENT (IN)    :: YDEV    ! unit no. for ar-to-point
        INTEGER , INTENT (IN)    :: CDEV    ! SCC descriptions unit no.
        INTEGER , INTENT (IN)    :: LDEV    ! log file unit no.

C...........   Local pointers
        INTEGER, POINTER :: OLDINDEXA ( : ) !  subscript table for SORTIC
        INTEGER, POINTER :: OLDIFIPA  ( : ) !  raw state/county FIPS code
        INTEGER, POINTER :: OLDTPFLGA ( : ) !  temporal resolution code
        INTEGER, POINTER :: OLDINVYRA ( : ) !  inventory year
        INTEGER, POINTER :: OLDSRCIDA ( : ) !  Source ID
        INTEGER, POINTER :: OLDIPOSCOD( : ) !  positn of pol in INVPCOD

        REAL   , POINTER :: OLDPOLVLA( :,: )   ! emission values

        CHARACTER(LEN=SCCLEN3), POINTER :: OLDCSCCA  ( : ) ! SCC
        CHARACTER(LEN=ALLCAS3), POINTER :: OLDCSOURCA( : ) ! concat src

C...........   Local allocatable arrays
        INTEGER               , ALLOCATABLE :: REPIDX ( : )    ! index for sorting
        REAL                  , ALLOCATABLE :: REPEMIS( :,: )  ! emissions for reporting
        CHARACTER(LEN=ALLCAS3), ALLOCATABLE :: REPSRCS( : )    ! source characteristics

C...........   Other local variables
        INTEGER         I,J,K,S     ! counters
        INTEGER         IOS         ! I/O error status
        INTEGER         LSTATE      ! last state code
        INTEGER         LPOL        ! last pollutant code
        INTEGER         NA2PSRCS    ! no. of area-to-point sources to add
        INTEGER         NREPSRCS    ! total no. of area-to-point sources
        INTEGER         TBLE        ! current area-to-point table number
        INTEGER         ROW         ! current area-to-point row
        INTEGER         OLDNRAWBP   ! old number of srcs x pols
        INTEGER         PE, PS      ! pollutant postn end and start in CSOURCA 
        INTEGER         POS         ! position in unsorted arrays
        INTEGER         REPPOS      ! position in reporting arrays
        INTEGER         SLEN        ! length of source 
        INTEGER         TSTATE      ! temp state code
        INTEGER         TPOL        ! temp pollutant code

        REAL            EANN        ! original annual emissions
        REAL            EOZN        ! original ozone season emissions
        REAL            PREVXLOC    ! previous x location
        REAL            PREVYLOC    ! previous y location

        CHARACTER(LEN=5      )  TPOLPOS     !  Temporary pollutant position
        CHARACTER(LEN=SCCLEN3)  LSCC        !  last SCC code
        CHARACTER(LEN=SCCLEN3)  TSCC        !  temp SCC code
        CHARACTER(LEN=ALLLEN3)  LSRCCHR     !  previous CSOURC
        CHARACTER(LEN=ALLLEN3)  TSRCCHR     !  tmporary CSOURC
        CHARACTER(LEN=3)        LOCID       !  tmporary location number
        CHARACTER(LEN=256    )  MESG        !  message buffer 

        CHARACTER*16 :: PROGNAME = 'ADJUSTINV' ! program name

C***********************************************************************
C   begin body of subroutine ADJUSTINV

C..........  If area-to-point factors file is present...
        IF( YDEV .GT. 0 ) THEN

C.............  Read and preprocess area-to-point factors file
C.............  Result of this call is that the NAR2PT and AR2PTABL 
C               arrays from MODAR2PT and the CHRT09 and ARPT09 arrays
C               from MODLISTS will be populated.
            CALL RDAR2PT( YDEV, CDEV, LDEV )

C.............  Assign area-to-point cross-reference entries to sources
C.............  Result of this call is that the AR2PTTBL, AR2PTIDX, and
C               AR2PTCNT arrays from MODLISTS will be populated
            CALL ASGNAR2PT( NRAWBP )
            
C.............  Determine total number of sources to be added; if source only
C               has one location, can use existing arrays and don't need to 
C               create additional memory for it
            NA2PSRCS = 0
            NREPSRCS = 0

            DO I = 1, NRAWBP
                S = SRCIDA( I )
                IF( AR2PTTBL( S ) /= 0 ) THEN
                    NA2PSRCS = NA2PSRCS + AR2PTCNT( S ) - 1
                    NREPSRCS = NREPSRCS + 1
                END IF
            END DO

C.............  Allocate memory for reporting
            ALLOCATE( REPIDX ( NREPSRCS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'REPIDX', PROGNAME )
            ALLOCATE( REPSRCS( NREPSRCS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'REPSRCS', PROGNAME )
            ALLOCATE( REPEMIS( NREPSRCS,2 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'REPEMIS', PROGNAME )

            OLDNRAWBP = NRAWBP
                
            IF( NA2PSRCS > 0 ) THEN

C.................  Update total number of raw records
                NRAWBP = NRAWBP + NA2PSRCS

C.................  Associate temporary pointers with unsorted arrays
                OLDINDEXA  => INDEXA
                OLDIFIPA   => IFIPA
                OLDTPFLGA  => TPFLGA
                OLDINVYRA  => INVYRA
                OLDCSCCA   => CSCCA
                OLDSRCIDA  => SRCIDA
                OLDIPOSCOD => IPOSCOD
                OLDCSOURCA => CSOURCA
                OLDPOLVLA  => POLVLA

C.................  Nullify original unsorted arrays
                NULLIFY( INDEXA, IFIPA, TPFLGA, INVYRA, CSCCA, 
     &                   SRCIDA, IPOSCOD, CSOURCA, POLVLA )

C.................  Allocate memory for larger unsorted arrays
                CALL SRCMEM( CATEGORY, 'UNSORTED', .TRUE., .FALSE., 
     &                       NRAWBP, NRAWBP, NPPOL )
     
                CALL SRCMEM( CATEGORY, 'UNSORTED', .TRUE., .TRUE.,
     &                       NRAWBP, NRAWBP, NPPOL )
                
                POLVLA = BADVAL3  ! array

C.................  Store old values in new arrays
                INDEXA ( 1:OLDNRAWBP ) = OLDINDEXA
                IFIPA  ( 1:OLDNRAWBP ) = OLDIFIPA
                TPFLGA ( 1:OLDNRAWBP ) = OLDTPFLGA
                INVYRA ( 1:OLDNRAWBP ) = OLDINVYRA
                CSCCA  ( 1:OLDNRAWBP ) = OLDCSCCA
                SRCIDA ( 1:OLDNRAWBP ) = OLDSRCIDA
                IPOSCOD( 1:OLDNRAWBP ) = OLDIPOSCOD
                CSOURCA( 1:OLDNRAWBP ) = OLDCSOURCA
                POLVLA ( 1:OLDNRAWBP,: ) = OLDPOLVLA

C.................  Deallocate old unsorted arrays
                DEALLOCATE( OLDINDEXA, OLDIFIPA,  
     &                      OLDTPFLGA, OLDINVYRA,  OLDCSCCA,  
     &                      OLDSRCIDA, OLDIPOSCOD, OLDCSOURCA, 
     &                      OLDPOLVLA   )
                
            END IF

C.............  Allocate memory for X and Y locations
            ALLOCATE( XLOCAA( NRAWBP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'XLOCAA', PROGNAME )
            ALLOCATE( YLOCAA( NRAWBP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'YLOCAA', PROGNAME )
            XLOCAA = BADVAL3  ! array
            YLOCAA = BADVAL3  ! array
            
C.............  Set position for adding sources to unsorted arrays
            POS = OLDNRAWBP + 1
            REPPOS = 0

C.............  Loop through original sources
            DO I = 1, OLDNRAWBP
                S = SRCIDA( I )
                J = INDEXA( I )

C.................  Check if current source is to be processed
                IF( AR2PTTBL( S ) /= 0 ) THEN
                	
C.....................  Save information for reports
                    REPPOS = REPPOS + 1
                    REPSRCS( REPPOS ) = CSOURCA( J )
                    REPEMIS( REPPOS,1 ) = POLVLA( J,NEM )
                    REPEMIS( REPPOS,2 ) = BADVAL3

C.....................  Set area-to-point table and row for current source                	
                    TBLE = AR2PTTBL( S )
                    ROW  = AR2PTIDX( S )

C.....................  Store X and Y locations
                    XLOCAA( J ) = AR2PTABL( ROW, TBLE )%LON
                    YLOCAA( J ) = AR2PTABL( ROW, TBLE )%LAT

C.....................  Save original annual and ozone season emission values
                    EANN = POLVLA( J, NEM )
                    EOZN = POLVLA( J, NOZ )

C.....................  Adjust annual and ozone season emissions based on allocation factor
                    IF( EANN /= BADVAL3 ) THEN
                        POLVLA( J, NEM ) = 
     &                          EANN * AR2PTABL( ROW,TBLE )%ALLOC
                        REPEMIS( REPPOS,2 ) = POLVLA( J, NEM )
                    END IF
                    
                    IF( EOZN /= BADVAL3 ) THEN
                        POLVLA( J, NOZ ) = 
     &                          EOZN * AR2PTABL( ROW,TBLE )%ALLOC
                    END IF

C.....................  Store location number in CSOURCA
                    WRITE( LOCID, '(I3)' ) 1
                    CSOURCA( J )( POLPOS3-3:POLPOS3-1 ) = 
     &                       ADJUSTR( LOCID )

C.....................  Loop through remaining locations if source has any
                    DO K = 1, AR2PTCNT( S ) - 1
                        ROW = ROW + 1

C.........................  Store source characteristics
                        IFIPA  ( POS ) = IFIPA  ( J )
                        TPFLGA ( POS ) = TPFLGA ( J )
                        INVYRA ( POS ) = INVYRA ( J )
                        CSCCA  ( POS ) = CSCCA  ( J )
                        CSOURCA( POS ) = CSOURCA( J )
                        POLVLA ( POS,NCE ) = POLVLA( J,NCE )
                        POLVLA ( POS,NRE ) = POLVLA( J,NRE )
                        POLVLA ( POS,NRP ) = POLVLA( J,NRP )

C.........................  Store X and Y locations
                        XLOCAA( POS ) = AR2PTABL( ROW,TBLE )%LON
                        YLOCAA( POS ) = AR2PTABL( ROW,TBLE )%LAT

C.........................  Adjust and store emissions
                        IF( EANN /= BADVAL3 ) THEN
                            POLVLA( POS,NEM ) =
     &                              EANN * AR2PTABL( ROW,TBLE )%ALLOC
                            REPEMIS( REPPOS,2 ) = 
     &                          REPEMIS( REPPOS,2 ) + POLVLA( POS,NEM )
                        END IF
                        
                        IF( EOZN /= BADVAL3 ) THEN
                            POLVLA( POS,NOZ ) = 
     &                              EOZN * AR2PTABL( ROW,TBLE )%ALLOC
                        END IF

C.........................  Store location number in CSOURCA
                        WRITE( LOCID, '(I3)' ) K+1
                        CSOURCA( POS )( POLPOS3-3:POLPOS3-1 ) =
     &                           ADJUSTR( LOCID )

C.........................  Increment position in unsorted arrays
                        POS = POS + 1

                    END DO   ! loop through additional locations

                END IF  ! check if source has any locations
            END DO  ! loop through sources

C.............  Process reporting arrays

C.............  Remove county information from source
            REPSRCS( : )( 4:6 ) = ' '

C.............  Sort source information
            DO I = 1, NREPSRCS
                REPIDX( I ) = I
            END DO

            CALL SORTIC( NREPSRCS, REPIDX, REPSRCS )

C.............  Determine total number of sources accounting for multiple
C               pollutants and counties
            LSTATE = 0
            LSCC = EMCMISS3
            LPOL = 0
            DO I = 1, NREPSRCS
                J = REPIDX( I )
            
                TSTATE = STR2INT( REPSRCS( J )( 2:3 ) )
                TSCC   = REPSRCS( J )( SCCPOS3:SCCPOS3+SCCLEN3-1 )
                TPOL   = STR2INT( REPSRCS( J )( POLPOS3:ALLLEN3 ) )
                
                IF( TSTATE /= LSTATE .OR.
     &              TSCC   /= LSCC   .OR.
     &              TPOL   /= LPOL        ) THEN
                    NCONDSRC = NCONDSRC + 1
                END IF
                
                LSTATE = TSTATE
                LSCC   = TSCC
                LPOL   = TPOL
            END DO

C.............  Allocate final size for reporting array
            ALLOCATE( REPAR2PT( NCONDSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'REPAR2PT', PROGNAME )
            REPAR2PT%NFIPS    = 0
            REPAR2PT%ORIGEMIS = 0.
            REPAR2PT%SUMEMIS  = 0.

C.............  Loop through and store final reporting information   
            LSTATE = 0
            LSCC = EMCMISS3
            LPOL = 0
            K    = 0                     
            DO I = 1, NREPSRCS
                J = REPIDX( I )
            
                TSTATE = STR2INT( REPSRCS( J )( 2:3 ) )
                TSCC   = REPSRCS( J )( SCCPOS3:SCCPOS3+SCCLEN3-1 )
                TPOL   = STR2INT( REPSRCS( J )( POLPOS3:ALLLEN3 ) )
                
                IF( TSTATE /= LSTATE .OR.
     &              TSCC   /= LSCC   .OR.
     &              TPOL   /= LPOL        ) THEN
                    
                    K = K + 1
                    REPAR2PT( K )%STATE = TSTATE
                    REPAR2PT( K )%SCC   = TSCC
                    REPAR2PT( K )%POLL  = TPOL
                    REPAR2PT( K )%NFIPS = 1
               
                    IF( REPEMIS( I,1 ) /= BADVAL3 ) THEN
                        REPAR2PT( K )%ORIGEMIS = REPEMIS( I,1 )
                    END IF
               
                    IF( REPEMIS( I,2 ) /= BADVAL3 ) THEN
                        REPAR2PT( K )%SUMEMIS  = REPEMIS( I,2 )
                    END IF
                ELSE
                    REPAR2PT( K )%NFIPS = REPAR2PT( K )%NFIPS + 1
                    
                    IF( REPEMIS( I,1 ) /= BADVAL3 ) THEN
                        REPAR2PT( K )%ORIGEMIS = 
     &                    REPAR2PT( K )%ORIGEMIS + REPEMIS( I,1 )
                    END IF
               
                    IF( REPEMIS( I,2 ) /= BADVAL3 ) THEN
                        REPAR2PT( K )%SUMEMIS  = 
     &                    REPAR2PT( K )%SUMEMIS  + REPEMIS( I,2 )
                    END IF
                END IF
        
                LSTATE = TSTATE
                LSCC   = TSCC
                LPOL   = TPOL
                         
            END DO

C.............  Deallocate temporary arrays
            DEALLOCATE( REPIDX, REPSRCS, REPEMIS )

C.............  Resort inventory and pollutants
            CALL M3MSG2( 'Resorting raw inventory data...' )
            
            DO I = 1, NRAWBP
                INDEXA( I ) = I
            END DO
            
            CALL SORTIC( NRAWBP, INDEXA, CSOURCA )

C.............  Remove location id from CSOURCA array
            CSOURCA( : )( POLPOS3-3:POLPOS3-1 ) = ' '

C.............  Reassign source numbers (copied from procinven.f)
            LSRCCHR = EMCMISS3
            PREVXLOC = BADVAL3
            PREVYLOC = BADVAL3
            S = 0
            SLEN  = SC_ENDP( MXCHRS )
            PS    = SC_BEGP( MXCHRS + 1 )
            PE    = SC_ENDP( MXCHRS + 1 )
            DO I = 1, NRAWBP
                
                J  = INDEXA( I )
           
                TSRCCHR = CSOURCA( J )(  1:SLEN ) ! Source characteristics
                TPOLPOS = CSOURCA( J )( PS:PE   ) ! Pos of pollutant (ASCII)
               
C.................  Update pointer for list of actual pollutants & activities
                K = STR2INT( TPOLPOS )  ! Convert pol/activity position to integer
                IPOSCOD( I ) = K

C.................  Increment source count by comparing this iteration to previous
                IF( TSRCCHR /= LSRCCHR ) THEN
                    S = S + 1
                    LSRCCHR = TSRCCHR
                ELSE
                    IF( XLOCAA( J ) /= PREVXLOC .OR.
     &                  YLOCAA( J ) /= PREVYLOC ) THEN
                        S = S + 1
                    END IF
                END IF

C.................  Save current X and Y locations           
                PREVXLOC = XLOCAA( J )
                PREVYLOC = YLOCAA( J )
                        
C.................  Assign source ID (to use as an index) for all inv X pol/act
                SRCIDA( I ) = S
           
            END DO  ! On sources x pollutants/activities

C.............  Update NSRC
            NSRC = S

        END IF

C..........  If non-HAP exclusions file is present...
        IF( UDEV .GT. 0 ) THEN

C.............  Read and preprocess NONHAPVOC exclusions x-ref
C.............  Only the CHRT* arrays of the MODXREF will be populated,
C               because we only need to identify the sources, not assign
C               anything to them.
            CALL RDXCLUDE( UDEV )

C.............   Assign array for non-HAP exclusions
            CALL ASGNNHAPX( NRAWBP )

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
 
        END SUBROUTINE ADJUSTINV
