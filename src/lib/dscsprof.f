
        SUBROUTINE DSCSPROF( FDEV, NIPOL, EINAM )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine determines the maximum number of species per pollutant,
C      determines the maximum number of profile table entries per pollutant, 
C      and creates a table that stores the species name per pollutant.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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

C...........   Modules for public variables
C...........   This module contains the speciation profile tables
        USE MODSPRO

        IMPLICIT NONE

C...........   Include files

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        
C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         GETFLINE
        INTEGER         INDEX1

        EXTERNAL        CRLF, GETFLINE, INDEX1

C...........   Subroutine arguments (note- outputs MXSPFUL, MXSPEC, and SPCNAMES
C              passed via module MODSPRO)

        INTEGER     , INTENT  (IN) :: FDEV            ! file unit number
        INTEGER     , INTENT  (IN) :: NIPOL           ! number of pollutants
        CHARACTER(*), INTENT  (IN) :: EINAM( NIPOL )  ! pollutant names
   
C...........   Arrays for getting pollutant-specific information from file
        INTEGER       NENTRA ( NIPOL )   ! number of table entries per pollutant
        INTEGER       NSPECA ( NIPOL )   ! number of species per pollutant
        CHARACTER*16  POLNAMA( NIPOL )   ! unsorted pollutant names
        LOGICAL       LMOLAR ( NIPOL )   ! true: moles conversion is not mass

C...........   Arrays for getting species-specific information from file
        INTEGER     , ALLOCATABLE :: INDX1A ( : ) ! sorting index for SPECNMA
        INTEGER     , ALLOCATABLE :: ISPOL  ( : ) ! assoc pol position in EINAM
        CHARACTER*16, ALLOCATABLE :: SPECNMA( : ) ! unsorted species names
                
C...........   Local variables

        INTEGER        I, J, K    ! counters and indices
        INTEGER        ICOUNT     ! tmp counter while populating SPCNAMES
        INTEGER        INPRFTP    ! tmp. profile number
        INTEGER        IOS        ! i/o status
        INTEGER        IPOL       ! pollutant counter
        INTEGER        IREC       ! record counter
        INTEGER        ISP        ! species names counter
        INTEGER        NLINES     ! number of lines in data file
        INTEGER        PPOS       ! tmp position (from INDEX1) of pol in POLNAMA
        INTEGER        SPOS       ! tmp position (from INDEX1) of pol in SPECNMA

        REAL           FAC1, FAC2 ! tmp speciation profile factors

        CHARACTER*5    TMPPRF     ! tmp profile number
        CHARACTER*16   POLNAM     ! pollutant name
        CHARACTER*16   SPECNM     ! tmp species name
        CHARACTER*300  MESG       ! message buffer
        
        CHARACTER*16 :: PROGNAME = 'DSCSPROF' ! program name

C***********************************************************************
C   Begin body of subroutine DSCSPROF

C...........  Make sure routine arguments are valid
        IF( FDEV .LE. 0 .OR. NIPOL .LE. 0 ) THEN
            MESG = 'INTERNAL ERROR: Invalid subroutine arguments'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
       
C...........   Determine length of input file and allocate memory for
C              a temporary species names array, an array that
C              associates a pollutant with each species name, and an
C              index array

        NLINES = GETFLINE( FDEV, 'Speciation profile file' )
        ALLOCATE( SPECNMA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPECNMA', PROGNAME )
        ALLOCATE( ISPOL( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISPOL', PROGNAME )
        ALLOCATE( INDX1A( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDX1A', PROGNAME )
        ALLOCATE( SPCNAMES( MXVARS3,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPCNAMES', PROGNAME )

C...........   Initialize species count per pollutant and flag for indicating
C              true molar conversions (NOTE - for some pollutants like PM10,
C              there is no mole-based factor and outputu should be in units
C              of gm/mole into the mole-base speciation matrix)
        NENTRA  = 0        ! array
        NSPECA  = 0.       ! array
        POLNAMA = ' '      ! array
        LMOLAR  = .TRUE.   ! array

C...........   Read through input file to determine the total number
C              of pollutants in the input file, to determine the
C              number of profiles per pollutant, to store the unique 
C              species names, and to store the units for mass-based and
C              mole-based conversions
        ICOUNT = 1
        IPOL   = 0
        IREC   = 0
        ISP    = 0
        DO I = 1, NLINES
        
            READ( FDEV,93100,END=999,IOSTAT=IOS ) 
     &                             TMPPRF, POLNAM, SPECNM, FAC1, FAC2
     
            IREC = IREC + 1
             
            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010 )
     &              'I/O error', IOS, 'reading speciation profile '//
     &              'file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Search for pollutant in list of valid names, and go to the end
C               of the loop if not found (skip entry)
            J = INDEX1( POLNAM, NIPOL, EINAM )
    
            IF ( J .EQ. 0 ) CYCLE
            
C.............  Search for pollutant unique list of all pollutants
            PPOS = INDEX1( POLNAM, IPOL, POLNAMA )
            
            IF ( PPOS .EQ. 0 ) THEN       ! if current POLNAM is not in
                                           !    POLNAMA, then
                IPOL = IPOL + 1
                POLNAMA( IPOL ) = POLNAM    ! add POLNAM to POLNAMA
                NENTRA ( IPOL ) = 1         ! init for first entry per pol

C.................  If mole-based factor is unity, then no molar transform
                IF( FAC1/FAC2 .EQ. 1. ) LMOLAR( J ) = .FALSE.

                PPOS = IPOL   ! Set for storing species count, below
               
            ELSE     ! if current POLNAM is already in POLNAMA, then

C.................  If a new profile number, then add to count of table entries
C                   for this pollutant
                NENTRA( PPOS ) = NENTRA( PPOS ) + 1
        
            END IF
            
            SPOS = INDEX1( SPECNM, ISP, SPECNMA )
        
            IF ( SPOS .LE. 0 ) THEN    ! if current SPECNM is not in
                                           ! SPECNMA, then
                ISP = ISP + 1
                INDX1A ( ISP )  = ISP
                SPECNMA( ISP )  = SPECNM
                ISPOL  ( ISP )  = J
     
            END IF

C.............  Check if species is already stored for current pollutant, and
C               if not, increment species-per-pollutant counter and 
C               add species to list.
            K = NSPECA( PPOS ) 
            J = INDEX1( SPECNM, K, SPCNAMES( 1,PPOS ) )

            IF( J .LE. 0 ) THEN
            
                K = K + 1

                IF( K .LE. MXVARS3 ) THEN
                    SPCNAMES( K, PPOS ) = SPECNM
                END IF

                NSPECA( PPOS ) = K

            END IF

        END DO

C...........   Determine the max species per full table, and the max species
C              per pollutant
        MXSPFUL = MAXVAL( NENTRA )
        MXSPEC  = MAXVAL( NSPECA )
       
C...........   Allocate memory for species names array and units to use for
C              mole-based transformations. Also initialize.
        DEALLOCATE( SPCNAMES )
        ALLOCATE( SPCNAMES( MXSPEC,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPCNAMES', PROGNAME )
        ALLOCATE( MOLUNITS( MXSPEC,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MOLUNITS', PROGNAME )

        SPCNAMES = ' '   ! array
        MOLUNITS = ' '   ! array

C...........   Sort master species names
        CALL SORTIC( ISP , INDX1A, SPECNMA )  ! sort on SPECNMA

C...........   Initialize species names table and counter per pollutant
        SPCNAMES = ' ' ! array

C...........   Cycle through count of all valid pollutants (NIPOL) and all 
C              species associated with these pollutants (ISP).  When the 
C              position of the pollutants in EINAM equals the pollutant 
C              position stored for the species, then write to SPCNAMES table.
        DO I = 1, NIPOL

            ICOUNT = 0
            DO J = 1, ISP
    
                K = INDX1A( J )
                IF ( ISPOL( K ) .EQ. I ) THEN 

                     ICOUNT = ICOUNT + 1
                     SPCNAMES( ICOUNT, I ) = SPECNMA( K )

C......................  When the pollutant does not have molar factors, store
C                        the molar units as mass units
                     IF( LMOLAR( I ) ) THEN
                         MOLUNITS( ICOUNT, I ) = SMOLUNIT
                     ELSE
                         MOLUNITS( ICOUNT, I ) = SMASUNIT
                     END IF

                END IF
       
            END DO

        END DO
 
C........  Rewind file

       REWIND( FDEV )
       
       RETURN
       
C......... Error message for reaching the end of file too soon
999    MESG = 'End of file reached unexpectedly. ' //
     &        'Check format of temporal' // CRLF() // BLANK5 //
     &        'cross reference file.'
       CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

       
C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93100   FORMAT( A5, 1X, A16, 1X, A16, 1X, G13.0, 1X, G13.0, 1X, G9.0 )

C...........   Internal buffering formats............ 94xxx

94010  FORMAT( 10( A, :, I8, :, 1X ) )
        
       END SUBROUTINE DSCSPROF                                                                            
