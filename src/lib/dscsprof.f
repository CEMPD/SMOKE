
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
C       Created ??/???? by ????
C       09/2025 by HT UNC-IE:  Use M3UTILIO; change  POLNAMA's size from 16 to IOVLEN3
C                             Change EINAM's size from (*) to (IOVLEN3) in consistent with MODINFO
C
C****************************************************************************/
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
C***************************************************************************
        USE M3UTILIO

C...........   Modules for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CRL

C...........   This module contains the speciation profile tables
        USE MODSPRO, ONLY: MXSPFUL, MXSPEC, SPCNAMES, MOLUNITS

        IMPLICIT NONE

C...........   Include files

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
c       INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        
C...........   EXTERNAL FUNCTIONS and their descriptions:
c       CHARACTER(2)    CRLF
c       INTEGER         GETFLINE
c       INTEGER         INDEX1
c       REAL            STR2REAL
c       LOGICAL         BLKORCMT

c       EXTERNAL        BLKORCMT, CRLF, GETFLINE, INDEX1, STR2REAL
        INTEGER, EXTERNAL :: GETFLINE
        LOGICAL, EXTERNAL :: BLKORCMT

C...........   Subroutine arguments (note- outputs MXSPFUL, MXSPEC, and SPCNAMES
C              passed via module MODSPRO)

        INTEGER     , INTENT  (IN) :: FDEV            ! file unit number
        INTEGER     , INTENT  (IN) :: NIPOL           ! number of pollutants
        CHARACTER(IOVLEN3), INTENT  (IN) :: EINAM( NIPOL )  ! pollutant names

C.........  Local parameters
        INTEGER, PARAMETER :: MXSEG = 6        ! # of potential line segments
        INTEGER, PARAMETER :: TMPNSPEC = 5000  ! # tmp number of species per pollutant

C...........   Arrays for getting pollutant-specific information from file
        INTEGER       NENTRA ( NIPOL )      ! number of table entries per pollutant
        INTEGER       NSPECA ( NIPOL )      ! number of species per pollutant
        CHARACTER(IOVLEN3) POLNAMA( NIPOL ) ! unsorted pollutant names; HT: change from 16 to IOVLEN3 for consistency

C...........   Arrays for getting species-specific information from file
        INTEGER,            ALLOCATABLE :: INDX1A ( : )    ! sorting index for SPECNMA
        CHARACTER(IOVLEN3), ALLOCATABLE :: SPECNMA ( : )   ! unsort spcs names
        CHARACTER(IOVLEN3), ALLOCATABLE :: TMPNAMES( :,: ) ! unsort names per pollutant
        LOGICAL,            ALLOCATABLE :: LMOLAR ( : )    ! true: moles conversion is not mass
                
        INTEGER        IPOS( 10 )       ! position in input pollutant list

C...........   Other arrays
        CHARACTER(32) SEGMENT( MXSEG )  ! Segments of parsed lines; HT: increase lenght to 32 to support variable name expansion 

C...........   Local variables

        INTEGER        I, J, K, L, M, N ! counters and indices
        INTEGER        ICOUNT     ! tmp counter while populating SPCNAMES
        INTEGER        INPRFTP    ! tmp. profile number
        INTEGER        IOS        ! i/o status
        INTEGER        IPOL       ! pollutant counter
        INTEGER        IREC       ! record counter
        INTEGER        ISP        ! species names counter
        INTEGER        NIPOS      ! number of pollutant matches
        INTEGER        NLINES     ! number of lines in data file
        INTEGER        PPOS       ! tmp position (from INDEX1) of pol in POLNAMA
        INTEGER        SPOS       ! tmp position (from INDEX1) of pol in SPECNMA

        REAL           FAC1, FAC2, FAC3 ! tmp speciation profile factors

        LOGICAL     :: EFLAG    = .FALSE.   ! true: error found
        LOGICAL     :: INHEADER = .FALSE.   ! true: in file header
        LOGICAL     :: BFLAG    = .FALSE.   ! true: running for biogenics

        CHARACTER(256) LINE       ! read buffer for a line
        CHARACTER(256) MESG       ! message buffer
        
        CHARACTER(SPNLEN3)  TMPPRF     ! tmp profile number
        CHARACTER(IOVLEN3)  POLNAM     ! tmp pollutant name
        CHARACTER(IOVLEN3)  SPECNM     ! tmp species name
        CHARACTER(SPNLEN3)  SPPRO      ! biogenics speciation profile to use

        CHARACTER(16) :: PROGNAME = 'DSCSPROF' ! program name

C***********************************************************************
C   Begin body of subroutine DSCSPROF
        
C...........  Make sure routine arguments are valid
        IF( FDEV .LE. 0 .OR. NIPOL .LE. 0 ) THEN
            MESG = 'INTERNAL ERROR: Invalid subroutine arguments'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C...........  Get source category, transfered by MODINFO module
        CALL GETCTGRY

C...........  For biogenics, evaluate speciation profile code to use
        IF( CRL == 'B' ) THEN
            BFLAG = .TRUE.
            MESG = 'Speciation profile to use for biogenics'
            CALL ENVSTR( 'BIOG_SPRO', MESG, ' ', SPPRO, IOS )

            IF( IOS .NE. 0 ) THEN
                MESG = 'ERROR: Variable BIOG_SPRO needs to be set ' //
     &                 'to run'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END IF

C...........   Determine length of input file and allocate memory for
C              a temporary species names array, an array that
C              associates a pollutant with each species name, an
C              index array, and an array to determine output units

        NLINES = GETFLINE( FDEV, 'Speciation profile file' )
        ALLOCATE( SPECNMA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPECNMA', PROGNAME )
        ALLOCATE( INDX1A( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDX1A', PROGNAME )
        ALLOCATE( TMPNAMES( TMPNSPEC,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPNAMES', PROGNAME )
        ALLOCATE( LMOLAR( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LMOLAR', PROGNAME )

C...........   Initialize species count per pollutant and flag for indicating
C              true molar conversions (NOTE - for some pollutants like PM10,
C              there is no mole-based factor and outputu should be in units
C              of gm/mole into the mole-base speciation matrix)
        NENTRA   = 0        ! array
        NSPECA   = 0.       ! array
        POLNAMA  = ' '      ! array
        TMPNAMES = ' '      ! array
        LMOLAR   = .FALSE.  ! array

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
        
            READ( FDEV,93000,END=999, IOSTAT=IOS ) LINE
     
            IREC = IREC + 1
             
            IF ( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &              'I/O error', IOS, 'reading speciation profile '//
     &              'file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Skip blank and comment lines

            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Separate the line of data into each part
            CALL PARSLINE( LINE, MXSEG, SEGMENT )

C.............  Left-justify character strings and convert factors to reals
            TMPPRF = ADJUSTL ( SEGMENT( 1 ) ) 

C.............  For biogenics, skip any entry that is not needed
            IF( BFLAG .AND. TMPPRF .NE. SPPRO ) CYCLE

            POLNAM = ADJUSTL ( SEGMENT( 2 ) )
            SPECNM = ADJUSTL ( SEGMENT( 3 ) )
            FAC1   = STR2REAL( SEGMENT( 4 ) )
            FAC2   = STR2REAL( SEGMENT( 5 ) )
            FAC3   = STR2REAL( SEGMENT( 6 ) )

C.............  Make sure divsor factor is not zero
            IF( FAC2 .EQ. 0. ) THEN
                WRITE( MESG,94010 ) 'WARNING: Zero divisor found at '//
     &                 'line ', IREC, '. Setting to 1.'
                CALL M3MESG( MESG )
                FAC2 = 1.
            END IF

C.............  Check width of character fields of fixed width
            L = LEN_TRIM( TMPPRF )
            IF( L .GT. SPNLEN3 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Speciation profile code ' //
     &                 'exceeds max width of', SPNLEN3, 'at line', IREC
                CALL M3MESG( MESG )
            END IF

            L = LEN_TRIM( POLNAM )
            IF( L .GT. IOVLEN3 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Pollutant name ' //
     &                 'exceeds max width of', IOVLEN3, 'at line', IREC
                CALL M3MESG( MESG )
            END IF

            L = LEN_TRIM( SPECNM )
            IF( L .GT. IOVLEN3 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Species name ' //
     &                 'exceeds max width of', IOVLEN3, 'at line', IREC
                CALL M3MESG( MESG )
            END IF

C.............  Search for pollutant in list of valid names, and go to the end
C               of the loop if none found (skip entry).  Record number
C               and position of all matches.
            M    = 0
            IPOS = 0   ! local array
            DO N = 1, NIPOL
                IF( POLNAM .EQ. EINAM( N ) ) THEN
                   M = M + 1
                   IF( M .LE. 10 ) THEN 
                       IPOS( M ) = N
                   ELSE
                       MESG = 'INTERNAL ERROR: IPOS array overflow '//
     &                        'prevented in ' // PROGNAME
                       CALL M3MSG2( MESG )
                       CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                   END IF 
                END IF
            END DO
            NIPOS = M
    
            IF ( MAXVAL( IPOS ) .EQ. 0 ) CYCLE
            
C.............  Search for pollutant unique list of all pollutants
            PPOS = INDEX1( POLNAM, IPOL, POLNAMA )
            
            IF ( PPOS .EQ. 0 ) THEN       ! if current POLNAM is not in
                                           !    POLNAMA, then
                IPOL = IPOL + 1
                POLNAMA( IPOL ) = POLNAM    ! add POLNAM to POLNAMA
                NENTRA ( IPOL ) = 1         ! init for first entry per pol

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

C.................  If mole-based = mass based, then use molar transform
                IF( FAC1/FAC2 .NE. FAC3 ) LMOLAR( ISP ) = .TRUE.
     
            END IF

C.............  Check if species is already stored for current pollutant, and
C               if not, increment species-per-pollutant counter and 
C               add species to list.
            DO M = 1, NIPOS

                K = NSPECA( IPOS( M ) )
                J = INDEX1( SPECNM, K, TMPNAMES( 1,IPOS( M ) ) )

                IF( J .LE. 0 ) THEN

                    K = K + 1

                    IF( K .LE. TMPNSPEC ) THEN
                        TMPNAMES( K, IPOS( M ) ) = SPECNM

                    ELSE

                       WRITE(MESG,94010)
     &                   'INTERNAL ERROR: The', TMPNSPEC, 'species '//
     &                   'per pollutant limit was exceeded by ' //
     &                   'pollutant '//TRIM(POLNAM)//' in ' // PROGNAME
                       CALL M3MSG2( MESG )
                       CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                    END IF

                    NSPECA( IPOS( M ) ) = K

                END IF
            END DO

        END DO              ! End loop over speciation profile input lines

        IF( IPOL .EQ. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: No pollutants found in speciation '//
     &             'profiles match the inventory!'
            CALL M3MSG2( MESG )
        END IF

        IF( ISP .EQ. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: No species found in speciation profile!'
            CALL M3MSG2( MESG )
        END IF

        IF( EFLAG ) THEN
            MESG = 'Problem(s) with speciation profiles file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C...........   Determine the max species per full table, and the max species
C              per pollutant
        MXSPFUL = MAXVAL( NENTRA )
        MXSPEC  = MAXVAL( NSPECA )
       
C...........   Allocate memory for species names array and units to use for
C              mole-based transformations. Also initialize.
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
C              species associated with these pollutants (ISP).  Check if species
C              is valid for the current pollutant, and if so, store in the
C              output species name list.
        DO I = 1, NIPOL

            ICOUNT = 0
            DO J = 1, ISP

C.................  Process species in sorted order    
                K = INDX1A( J )

C.................  Find species in list of valid species per pollutant
                N = INDEX1( SPECNMA( K ), NSPECA( I ),
     &                      TMPNAMES( 1,I )            ) 

                IF ( N .GT. 0 ) THEN 

                     ICOUNT = ICOUNT + 1
                     SPCNAMES( ICOUNT, I ) = SPECNMA( K )

C......................  When the species does not have molar factors, store
C                        the molar units as mass units
                     IF( LMOLAR( K ) ) THEN
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
     &        'Check format of speciation' // CRLF() // BLANK5 //
     &        'profile file.'
       CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

       
C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010  FORMAT( 10( A, :, I8, :, 1X ) )
        
       END SUBROUTINE DSCSPROF                                                                            
