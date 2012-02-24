
        SUBROUTINE RDCODNAM( FDEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C     Reads the inventory table to obtain a list of valid pollutants and
C     activities and their properties. The routine converts all names to 
C     uppercase and allocates the MODLISTS arrays needed for storing the
C     inventory table. It creates sorted lists of CAS numbers, both with
C     duplicates and without. The subroutine returns the names in their
C     original, unsorted order. 
C
C  PRECONDITIONS REQUIRED:
C     File unit FDEV already is opened
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C     Subroutines: Models-3 subroutines
C     Functions: Models-3 functions
C
C  REVISION  HISTORY:
C     Created 10/98 by M. Houyoux
C     Updated Aug 2006 by M. Houyoux for SMOKE 2.3 to add mode and change
C         internal code documentation from SAROAD to SPECIATE4 ID.
C     Updated Sep 2006 by M. Houyoux to automatically add _NOI pollutants
C         for entries that are no-integrate model species.
C
C****************************************************************************
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

C.........  MODULES for public variables
C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: NINVTBL, ITIDXA, ITIDXA2, ITLINNO, ITCODA,
     &                      ITNTIA, ITREAA, ITSTATA, ITKEEPA, ITMSPC,
     &                      ITEXPL, ITFACA, ITVTSA, ITNAMA, ITUNTA,
     &                      ITCASA, ITDSCA, ITCASDSCA, ITCASDNMA,
     &                      NINVKEEP, NINVDROP, NUNIQCAS, SCASIDX,
     &                      SORTCAS, UCASIDX, UCASNPOL, UCASNKEP,
     &                      UNIQCAS, EMISBYPOL, EMISBYCAS, RECSBYCAS,
     &                      MXIDAT, INVDCOD, INVSTAT, INVDCNV, INVDVTS,
     &                      INVDNAM, INVDUNT, INVDDSC, IDXCOD, SORTCOD

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL         CHKINT
        LOGICAL         CHKREAL
        LOGICAL         BLKORCMT
        CHARACTER(2)    CRLF
        INTEGER         GETFLINE
        INTEGER         INDEX1
        INTEGER         STR2INT
        REAL            STR2REAL

        EXTERNAL CHKINT, CHKREAL, CRLF, GETFLINE, INDEX1, STR2INT, 
     &           STR2REAL, BLKORCMT

C...........   SUBROUTINE ARGUMENTS
        INTEGER , INTENT (IN) :: FDEV   ! iventory table unit no.

C...........   Parameters
        INTEGER, PARAMETER :: NFIELDS = 15  ! no. input fields
        INTEGER, PARAMETER :: FBEG( NFIELDS ) = 
     &                      ( / 1 , 13, 17, 34, 40, 42, 44,
     &                          50, 52, 54, 56, 58, 62, 79, 119 / )
        INTEGER, PARAMETER :: FEND( NFIELDS ) = 
     &                      ( / 11, 15, 32, 38, 40, 42, 49,
     &                          50, 52, 54, 56, 60, 77, 118, 158 / )
        CHARACTER(4),       PARAMETER :: NOIEND = '_NOI'

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: LOCATIDX  ( : ) ! sorting index
        INTEGER, ALLOCATABLE :: INVDCODA( : ) ! 5-digit SPECIATE4 ID (if any)
        INTEGER, ALLOCATABLE :: INVSTATA( : ) ! Status (<0 activity; >0 pol)

        LOGICAL, ALLOCATABLE :: LADDNOI ( : ) ! true: add NOI entry for this pollutant

        REAL   , ALLOCATABLE :: INVPOSA ( : ) ! position in original list

        CHARACTER         , ALLOCATABLE :: INVDVTSA( : ) ! V=VOC, T=TOG, N=not
        CHARACTER(IOULEN3), ALLOCATABLE :: INVDUNTA( : ) ! units for SMOKE intmdt inventory
        CHARACTER(CASLEN3), ALLOCATABLE :: INVDCASA( : ) ! CAS nubmer
        CHARACTER(DDSLEN3), ALLOCATABLE :: INVDDSCA( : ) ! inventory data description

C...........   Local arrays
        INTEGER       FLEN   ( NFIELDS )
        CHARACTER(64) SEGMENT( NFIELDS )

C...........   Other local variables
        INTEGER         J, K, L, N     !  counters and indices

        INTEGER         CNT     !  count of pollutants per CAS
        INTEGER         CNTK    !  count of kept pollutants per CAS
        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  record number
        INTEGER         NDAT    !  number of data values
        INTEGER         NNOI    !  number of "no integrate" pollutants

        LOGICAL       :: AFLAG    = .FALSE.  ! true: duplicate activitiy found
        LOGICAL       :: CFLAG    = .FALSE.  ! true: duplicate CAS number found
        LOGICAL       :: NFLAG    = .FALSE.  ! true: duplicate unit found
        LOGICAL       :: SFLAG    = .FALSE.  ! true: duplicate SPECIATE 4 found
        LOGICAL       :: VFLAG    = .FALSE.  ! true: duplicate V/T found
        LOGICAL       :: EFLAG    = .FALSE.  ! true:  found
        LOGICAL          KEEPVAL             ! initialized keep status

        CHARACTER             BUF1    !  tmp 1-char buffer
        CHARACTER(16)      :: FNAME = 'INVTABLE' ! Logical name of input file
        CHARACTER(256)        MESG    !  message buffer
        CHARACTER(512)     :: LINE    !  Input line
        CHARACTER(CASLEN3) :: PCAS    !  previous CAS code
        CHARACTER(IOVLEN3) :: LNAM    !  previous pollutant name
        CHARACTER(IOVLEN3) :: DNAM    !  tmp for pollutant or activity name
        CHARACTER(DDSLEN3) :: DSCVAL  !  initialized CAS description value

        CHARACTER(16) :: PROGNAME = 'RDCODNAM' ! program name

C***********************************************************************
C   begin body of subroutine RDCODNAM

C.........  Compute field lengths
        DO N = 1, NFIELDS
            FLEN( N ) = FEND( N ) - FBEG( N ) + 1
        END DO

C.........  Get file sizes and allocate memory...

C.........  Get no. lines in pollutant codes & activities files for allocating
C           memory     
        IF( FDEV .GT. 0 ) THEN
            NINVTBL = GETFLINE( FDEV, 'Inventory data table file' )
        END IF
      
C.........  Allocate memory for storing raw unsorted inventory table
        ALLOCATE( ITIDXA( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ITIDXA', PROGNAME )
        ALLOCATE( ITIDXA2( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ITIDXA2', PROGNAME )
        ALLOCATE( ITLINNO( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ITLINNO', PROGNAME )
        ALLOCATE( ITCODA( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ITCODA', PROGNAME )
        ALLOCATE( ITNTIA( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ITNTIA', PROGNAME )
        ALLOCATE( ITREAA( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ITREAA', PROGNAME )
        ALLOCATE( ITSTATA( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ITSTATA', PROGNAME )
        ALLOCATE( ITKEEPA( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ITKEEPA', PROGNAME )
        ALLOCATE( ITEXPL( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ITEXPL', PROGNAME )
        ALLOCATE( ITMSPC( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ITMSPC', PROGNAME )
        ALLOCATE( ITFACA( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ITFACA', PROGNAME )
        ALLOCATE( ITVTSA( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ITVTSA', PROGNAME )
        ALLOCATE( ITNAMA( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ITNAMA', PROGNAME )
        ALLOCATE( ITUNTA( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ITUNTA', PROGNAME )
        ALLOCATE( ITCASA( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ITCASA', PROGNAME )
        ALLOCATE( ITDSCA( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ITDSCA', PROGNAME )
        ALLOCATE( ITCASDSCA( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ITCASDSCA', PROGNAME )
        ALLOCATE( ITCASDNMA( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ITCASDNMA', PROGNAME )

        ITIDXA  = 0
        ITIDXA2 = 0
        ITLINNO = 0
        ITCODA  = 0
        ITNTIA  = -9            ! default missing
        ITREAA  = -9            ! default missing
        ITSTATA = 1             ! array (expected by PROCINVEN)
        ITKEEPA = .FALSE.       ! all pollutants dropped unless otherwise set in INVTABLE
        ITMSPC  = .FALSE.       ! all pollutants no model species
        ITEXPL  = .FALSE.       ! all pollutants not explicit in mechanism
        ITFACA  = 1.            ! initialize to have no change
        ITVTSA  = 'N'           ! initialize as not part of VOC or TOG
        ITNAMA  = ' '           ! names are blank
        ITUNTA  = ' '           ! blank
        ITCASA  = ' '           ! blank
        ITDSCA  = ' '           ! pollutant description blank
        ITCASDSCA = ' '         ! CAS description blank
        ITCASDNMA = ' '         ! CAS number with data name

C.........  Read Inventory Table
        NDAT = 0
        IREC = 0
        DO

            READ( FDEV, 93000, END=22, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010 ) 
     &                 'ERROR: System I/O error', IOS, 'reading ' // 
     &                 'inventory table file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Skip comment and blank lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Check if line is a process/pollutant combination
C            IF( INDEX( LINE( 1:16 ), ETJOIN ) > 0 ) THEN
C                SEGMENT( 1 ) = ADJUSTL( LINE( 1:16 ) )

C.................  Skip CAS number, SPECIATE4 code, and reactivity, then
C                   store remaining fields
C  commented out by Marc Houyoux for SMOKE2.3 - unsure why they were being skipped, but
C  they are needed for the case when process-pollutant combinations are
C  provided in the emissions inventory.
C                DO N = 6, NFIELDS
C                    SEGMENT( N ) =
C     &                  ADJUSTL( LINE( FBEG( N ):FEND( N ) ) )
C                END DO
C            ELSE

C.................  Parse the line into its sections based on file format def'n
                DO N = 1, NFIELDS

                    SEGMENT( N ) = 
     &                  ADJUSTL( LINE( FBEG( N ):FEND( N ) ) )

                END DO
C            END IF

C.............  Get length of data name
            L = LEN_TRIM( SEGMENT( 1 ) )

C.............  Error if data name is not provided
            IF( L .LE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Data name is blank '//
     &                 'at line', IREC
                CALL M3MSG2( MESG )

            END IF

C.............  Perform checks on input fields...
C.............  Check that data name does not have illegal characters
            DO N = 1, L
                BUF1 = SEGMENT( 1 ) (L:L) 
                IF( ( BUF1 .GE. '0' .AND. BUF1 .LE. '9' ) .OR.
     &              ( BUF1 .GE. 'A' .AND. BUF1 .LE. 'Z' ) .OR.
     &              ( BUF1 .EQ. '_' ) ) THEN
                 ! all is well
                ELSE
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )'ERROR: Disallowed character "'//
     &                 BUF1 // '" in data name at line', IREC, '.' //
     &                 CRLF()// BLANK10// 'Allowed characters are '//
     &                 '"0" to "9", "A" to "Z", and "_".'
                    CALL M3MSG2( MESG )
                END IF
            END DO

C.............  Check that integer fields are integers
            IF( .NOT. CHKINT( SEGMENT( 4 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: SPECIATE4 ID is not an ' //
     &                 'integer at line', IREC
                CALL M3MSG2( MESG )
            END IF

            IF( .NOT. CHKINT( SEGMENT( 5 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Reactivity group code ' //
     &                 'is not an integer at line', IREC
                CALL M3MSG2( MESG )
            END IF

            IF( .NOT. CHKINT( SEGMENT( 12 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: NTI code is not an ' //
     &                 'integer at line', IREC
                CALL M3MSG2( MESG )
            END IF

C.............  Check that float fields are floats
            IF( .NOT. CHKREAL( SEGMENT( 7 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Factor value is not a '//
     &                 'floating point value at line', IREC
                CALL M3MSG2( MESG )
            END IF

C.............  Convert to upper case, where needed
            CALL UPCASE( SEGMENT( 1 )( 1:FLEN(1) ) )
            CALL UPCASE( SEGMENT( 2 )( 1:FLEN(2) ) )
            CALL UPCASE( SEGMENT( 3 )( 1:FLEN(3) ) )
            CALL UPCASE( SEGMENT( 6 )( 1:FLEN(6) ) )
            CALL UPCASE( SEGMENT( 7 )( 1:FLEN(7) ) )
            CALL UPCASE( SEGMENT( 8 )( 1:FLEN(8) ) )
            CALL UPCASE( SEGMENT( 9 )( 1:FLEN(9) ) )
            CALL UPCASE( SEGMENT( 10)( 1:FLEN(10) ) )
            CALL UPCASE( SEGMENT( 11)( 1:FLEN(11) ) )

C.............  Correct unknown entries
            IF( SEGMENT( 8 ) .NE. 'V' .AND.
     &          SEGMENT( 8 ) .NE. 'T'       ) SEGMENT( 8 ) = 'N'
            IF( SEGMENT( 9 ) .NE. 'Y'       ) SEGMENT( 9 ) = 'N'
            IF( SEGMENT( 10) .NE. 'Y'       ) SEGMENT( 10) = 'N'
            IF( SEGMENT( 11) .NE. 'Y'       ) SEGMENT( 11) = 'N'

C.............  Error for blank units
            IF( SEGMENT( 13 ) .EQ. ' ' ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &                 'ERROR: Inventory units not set in the ' //
     &                 FNAME // ' file at line', IREC
                CALL M3MSG2( MESG )

            END IF

C.............  Store unsorted variables
            NDAT = NDAT + 1
            IF( NDAT .LE. NINVTBL .AND. .NOT. EFLAG ) THEN

                ITIDXA   ( NDAT ) = NDAT
                ITIDXA2  ( NDAT ) = NDAT
                ITLINNO  ( NDAT ) = IREC
                ITNAMA   ( NDAT ) = TRIM( SEGMENT( 1 ) )
                ITCASA   ( NDAT ) = TRIM( SEGMENT( 3 ) )
                ITCODA   ( NDAT ) = STR2INT( SEGMENT( 4 ) )
                ITREAA   ( NDAT ) = STR2INT( SEGMENT( 5 ) )
                ITKEEPA  ( NDAT ) = ( SEGMENT( 6 ) .EQ. 'Y' )
                IF( SEGMENT( 7 ) == '' ) SEGMENT( 7 ) = '1.0'
                ITFACA   ( NDAT ) = STR2REAL( SEGMENT( 7 ) )
                ITVTSA   ( NDAT ) = TRIM( SEGMENT( 8 ) )
                ITMSPC   ( NDAT ) = ( SEGMENT( 9 ) .EQ. 'Y' )
                ITEXPL   ( NDAT ) = ( SEGMENT( 10 ) .EQ. 'Y' )
                ITNTIA   ( NDAT ) = STR2INT( SEGMENT( 12 ) )
                ITUNTA   ( NDAT ) = TRIM( SEGMENT( 13 ) )
                ITDSCA   ( NDAT ) = TRIM( SEGMENT( 14 ) )
                ITCASDSCA( NDAT ) = TRIM( SEGMENT( 15 ) )
                ITCASDNMA( NDAT ) = ITCASA( NDAT ) // ITNAMA( NDAT )

C.................  Set name depending on whether or not the mode is
C                   included in the format.
                L = LEN_TRIM( SEGMENT( 2 ) )

C.................  If no mode provided:
                IF ( L .LE. 0 ) THEN
                    ITNAMA( NDAT ) = TRIM( SEGMENT( 1 ) )

C.................  If the mode is provided
                ELSE
                    ITNAMA( NDAT ) = TRIM( SEGMENT( 2 ) ) // ETJOIN //
     &                               TRIM( SEGMENT( 1 ) )
                END IF

                IF( SEGMENT( 11 ) .EQ. 'Y' ) ITSTATA( NDAT ) = -1
                IF( ITKEEPA( NDAT ) ) NINVKEEP = NINVKEEP + 1

            END IF

        END DO          !  end read loop

22      CONTINUE        !  exit from loop reading FDEV

C.........  Check dimensions
        IF( NDAT .GT. NINVTBL ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 
     &             'ERROR: Number of pollutant records :', NDAT, 
     &             CRLF() // BLANK10 // 'Memory allocated :', NINVTBL
            CALL M3MSG2( MESG )

        ELSE IF( NDAT .EQ. 0 ) THEN
            EFLAG = .TRUE.
            MESG ='ERROR: No entries in inventory table file.'
            CALL M3MSG2( MESG )

        END IF

C.........  Abort if error found during read
        IF ( EFLAG ) THEN
            MESG = 'Problem with INVTABLE input file contents.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Exclude comment lines from count of inventory table entries and
C           calculate dropped entries
        NINVTBL = NDAT
        NINVDROP = NINVTBL - NINVKEEP

C.........  Sort by CAS number then pollutant name
        CALL SORTIC( NINVTBL, ITIDXA, ITCASDNMA )

C.........  Count the number of unique CAS numbers.  Blank CAS numbers
C           will be first in the sorted list, so they will not be counted.
        PCAS =  ' ' 
        DO N = 1, NINVTBL
            J = ITIDXA( N )
            IF ( ITCASA( J ) .NE. PCAS ) NUNIQCAS = NUNIQCAS + 1
            PCAS = ITCASA( J )
        END DO

C.........  Allocate memory for sorted CAS number arrays
        ALLOCATE( SCASIDX( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SCASIDX', PROGNAME )
        ALLOCATE( SORTCAS( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SORTCAS', PROGNAME )
        ALLOCATE( UCASIDX( NUNIQCAS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UCASIDX', PROGNAME )
        ALLOCATE( UCASNPOL( NUNIQCAS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UCASNPOL', PROGNAME )
        ALLOCATE( UCASNKEP( NUNIQCAS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UCASNKEP', PROGNAME )
        ALLOCATE( UNIQCAS( NUNIQCAS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UNIQCAS', PROGNAME )
        
C.........  Allocate memory for reporting of emissions (used in reader routines)
        ALLOCATE( EMISBYPOL( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMISBYPOL', PROGNAME )
        ALLOCATE( EMISBYCAS( NUNIQCAS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMISBYCAS', PROGNAME )
        ALLOCATE( RECSBYCAS( NUNIQCAS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RECSBYCAS', PROGNAME )        

C.........  Initialize new arrays
        SCASIDX = 0    ! array
        SORTCAS = ' '  ! array
        UCASIDX = 0    ! array
        UCASNPOL= 0    ! array
        UCASNKEP= 0    ! array
        UNIQCAS = ' '  ! array
        EMISBYPOL = 0. ! array
        EMISBYCAS = 0. ! array
        RECSBYCAS = 0  ! array

C.........  Store sorted CAS numbers in both sorted order arrays
        PCAS =  ' ' 
        NUNIQCAS = 0
        CNT      = 0
        CNTK     = 0
        DO N = 1, NINVTBL
            J    = ITIDXA( N )
            IREC = ITLINNO( J )

            IF ( ITCASA( J ) .EQ. ' ' ) CYCLE

C.............  Increment unique count if this CAS is different than previous
C.............  Also, reset count of records per CAS and keep status
            IF ( ITCASA( J ) .NE. PCAS ) THEN
                NUNIQCAS = NUNIQCAS + 1
                KEEPVAL  = ITKEEPA( J )
                DSCVAL   = ITCASDSCA( J )
                CNT      = 0
                CNTK     = 0
            END IF

C.............  Give warning if keep status for this iteration is not the same as
C               that for the current CAS code
            IF( ( ITKEEPA( J ) .AND. .NOT. KEEPVAL ) .OR. 
     &          ( .NOT. ITKEEPA( J ) .AND. KEEPVAL )      ) THEN
                BUF1 = 'N'
                IF( KEEPVAL ) BUF1 = 'Y'
                WRITE( MESG,94010 ) 'WARNING: Keep status at '//
     &                 'line', IREC, 'different from previously '//
     &                 'set status of "'// BUF1// '" for CAS "' //
     &                 TRIM( ITCASA( J ) ) // '".'
                CALL M3MSG2( MESG )

            END IF

C.............  Store the data units if not yet set, otherwise, check that its the same
            IF( ITCASDSCA( J ) .NE. DSCVAL ) THEN
                WRITE( MESG,94010 ) 'WARNING: Different '//
     &                 'descriptions for the same CAS number at line',
     &                 IREC, '.'
                CALL M3MSG2( MESG )
            END IF

C.............  Sorted arrays, with duplicates
            SCASIDX( N ) = J
            SORTCAS( N ) = ITCASA( J )

C.............  Sorted arrays, without duplicates
            CNT = CNT + 1
            IF( ITKEEPA( J ) ) CNTK = CNTK + 1

            IF( UCASIDX( NUNIQCAS ) .EQ. 0 ) UCASIDX( NUNIQCAS ) = N
            UCASNPOL( NUNIQCAS ) = CNT
            UCASNKEP( NUNIQCAS ) = CNTK
            UNIQCAS ( NUNIQCAS ) = ITCASA( J ) 

C.............  Set previous CAS for next iteration
            PCAS = ITCASA( J )

        END DO

C.........  Sort inventory table to create unique data list and associated arrays
        CALL SORTIC( NINVTBL, ITIDXA2, ITNAMA )

C.........  Count the number of unique data names that are kept
        LNAM = ' '
        NNOI = 0
        MXIDAT = 0
        DO N = 1, NINVTBL
            J = ITIDXA2( N ) 
            IF( .NOT. ITKEEPA( J ) ) CYCLE

C.............  Check if this pollutant  name is not the same as the previous one            
            IF( ITNAMA( J ) .NE. LNAM ) THEN 

C.................  Count the number of unique data names
                MXIDAT = MXIDAT + 1

C.................  If pollutant is a part of VOC or TOG, is "no-integrate" 
C                   and a model species, increase the count for inserting 
C                   the _NOI pollutants.
                IF ( ( ITVTSA(J) .EQ. 'V' .OR. ITVTSA(J) .EQ. 'T' ).AND.
     &               ITMSPC( J ) .AND. .NOT. ITEXPL( J ) ) THEN

                    NNOI = NNOI + 1

                ENDIF

            END IF

            LNAM = ITNAMA( J )
        END DO

        MXIDAT = MXIDAT + NNOI

C.........  Allocate memory for local inventory pollutants/activities
        ALLOCATE( LOCATIDX( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LOCATIDX', PROGNAME )
        ALLOCATE( INVPOSA( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVPOSA', PROGNAME )
        ALLOCATE( INVDCODA( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDCODA', PROGNAME )
        ALLOCATE( INVSTATA( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVSTATA', PROGNAME )
        ALLOCATE( LADDNOI( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LADDNOI', PROGNAME )
        ALLOCATE( INVDVTSA( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDVTSA', PROGNAME )
        ALLOCATE( INVDUNTA( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDUNTA', PROGNAME )
        ALLOCATE( INVDCASA( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDCASA', PROGNAME )
        ALLOCATE( INVDDSCA( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDDSCA', PROGNAME )

        LOCATIDX = 0    ! array
        INVPOSA  = 0.    ! array
        INVDCODA = 0    ! array
        INVSTATA = 0    ! array
        LADDNOI  = .FALSE. ! array
        INVDVTSA = ' '  ! array
        INVDUNTA = ' '  ! array
        INVDCASA = ' '  ! array
        INVDDSCA = ' '  ! array

C.........  Allocate memory for output inventory pollutants/activities
        ALLOCATE( INVDCOD( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDCOD', PROGNAME )
        ALLOCATE( INVSTAT( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVSTAT', PROGNAME )
        ALLOCATE( INVDCNV( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDCNV', PROGNAME )
        ALLOCATE( INVDVTS( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDVTS', PROGNAME )
        ALLOCATE( INVDNAM( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDNAM', PROGNAME )
        ALLOCATE( INVDUNT( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDUNT', PROGNAME )
        ALLOCATE( INVDDSC( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDDSC', PROGNAME )

        INVDCOD = 0    ! array
        INVSTAT = 0    ! array
        INVDCNV = 1.   ! array
        INVDVTS = ' '  ! array
        INVDNAM = ' '  ! array
        INVDUNT = ' '  ! array
        INVDDSC = ' '  ! array

C.........  Allocate memory for sorted SPECATE4 IDs and index
        ALLOCATE( IDXCOD( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDXCOD', PROGNAME )
        ALLOCATE( SORTCOD( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SORTCOD', PROGNAME )

        IDXCOD  = 0  ! array
        SORTCOD = 0  ! array
      
C.........  Store kept pollutants and activities in original (unsorted) order
        LNAM = ' '
        MXIDAT = 0
        DO N = 1, NINVTBL
            J = ITIDXA2( N ) 

C.............  If current data name is not kept, skip to next
            IF( .NOT. ITKEEPA( J ) ) CYCLE
            
            IREC = ITLINNO( J )

C.............  For the next data name...
            IF( ITNAMA( J ) .NE. LNAM ) THEN
                MXIDAT = MXIDAT + 1

C.................  Reset keep status for current inventory data value
                KEEPVAL = ITKEEPA( J )

            END IF

C.............  Store first position in original table
            IF( INVPOSA( MXIDAT ) .EQ. 0. ) THEN
                 LOCATIDX ( MXIDAT ) = MXIDAT
                 INVPOSA( MXIDAT ) = REAL( 
     &                    INDEX1( ITNAMA( J ), NINVTBL, ITNAMA ) )
            END IF
  
C.............  Store SPECIATE4 ID if not yet set, otherwise, check that its the same code
            IF( INVDCODA( MXIDAT ) .EQ. 0 ) THEN
                INVDCODA( MXIDAT ) = ITCODA( J )

            ELSE IF( ITCODA( J ) .NE. INVDCODA( MXIDAT ) ) THEN
                WRITE( MESG,94010 ) 'WARNING: Different SPECIATE4 ' //
     &                 'IDs for the same data name at line', IREC,
     &                 CRLF()// BLANK10// 'Using code', 
     &                 INVDCODA( MXIDAT ), 'and ignoring code',
     &                 ITCODA( J )
                CALL M3MSG2( MESG )

            ELSE IF( ITCODA( J ) .EQ. INVDCODA( MXIDAT ) ) THEN
                SFLAG = .TRUE.

            END IF

C.............  Store activity status if not yet set, otherwise, check that its the same
            IF( INVSTATA( MXIDAT ) .EQ. 0 ) THEN
                INVSTATA( MXIDAT ) = ITSTATA( J )

            ELSE IF( ITSTATA( J ) .NE. INVSTATA( MXIDAT ) ) THEN
                EFLAG = .TRUE.
                BUF1 = 'N'
                IF( INVSTATA( MXIDAT ) .LT. 0 ) BUF1 = 'Y'
                WRITE( MESG,94010 ) 'ERROR: Activity status at '//
     &                 'line', IREC, 'different from previously '//
     &                 'set status of "'// BUF1// '".' 
                CALL M3MSG2( MESG )

            ELSE IF( ITSTATA( J ) .EQ. INVSTATA( MXIDAT ) ) THEN
                AFLAG = .TRUE.

            END IF

C.............  Store VOC/TOG status if not yet set, otherwise, check that its the same
            IF( INVDVTSA( MXIDAT ) .EQ. ' ' ) THEN
                INVDVTSA( MXIDAT ) = ITVTSA( J )

            ELSE IF( ITVTSA( J ) .NE. INVDVTSA( MXIDAT ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: VOC/TOG status "'//
     &                 ITVTSA( J ) // '" at line', IREC, 
     &                 'is different from previously '//
     &                 'set status of "' // INVDVTSA( MXIDAT )// '".'
                CALL M3MSG2( MESG )

            ELSE IF( ITVTSA( J ) .EQ. INVDVTSA( MXIDAT ) ) THEN
                VFLAG = .TRUE.

            END IF

C.............  Store the data description if not yet set, otherwise, check that its the same
            IF( INVDUNTA( MXIDAT ) .EQ. ' ' ) THEN
                INVDUNTA( MXIDAT ) = ITUNTA( J )

            ELSE IF( ITUNTA( J ) .NE. INVDUNTA( MXIDAT ) ) THEN
                WRITE( MESG,94010 ) 'WARNING: Different units ' //
     &                 'for the same data name at line', IREC, '.'//
     &                 CRLF()// BLANK10// 'Using units "' // 
     &                 TRIM( INVDUNTA( MXIDAT ) ) // '" and ' //
     &                 'ignoring units "' // TRIM( ITUNTA( J ) ) //
     &                 '".'
                CALL M3MSG2( MESG )

            ELSE IF( ITUNTA( J ) .EQ. INVDUNTA( MXIDAT ) ) THEN
                NFLAG = .TRUE.

            END IF

C.............  Store the CAS number if not yet set, otherwise, check that its the same
            IF( INVDCASA( MXIDAT ) .EQ. ' ' ) THEN
                INVDCASA( MXIDAT ) = ITCASA( J )

            ELSE IF( ITCASA( J ) .EQ. INVDCASA( MXIDAT ) ) THEN
                CFLAG = .TRUE.

            END IF

C.............  Store the data units if not yet set, otherwise, check that its the same
            IF( INVDDSCA( MXIDAT ) .EQ. ' ' ) THEN
                INVDDSCA( MXIDAT ) = ITDSCA( J )

            ELSE IF( ITDSCA( J ) .NE. INVDDSCA( MXIDAT ) ) THEN
                WRITE( MESG,94010 ) 'WARNING: Different descriptions '//
     &                 'for the same data name at line', IREC, '.'//
     &                 CRLF()// BLANK10// 'Using description:' // 
     &                 CRLF()// BLANK16// '"'// 
     &                 TRIM( INVDDSCA(MXIDAT) ) // '" and ignoring:' //
     &                 CRLF()// BLANK16// '"'// 
     &                 TRIM( ITDSCA( J ) ) // '".'
                CALL M3MSG2( MESG )

            END IF

C.............  Give warning if keep status for this iteration is not the same as
C               that for the current data name
            IF( ( ITKEEPA( J ) .AND. .NOT. KEEPVAL ) .OR. 
     &          ( .NOT. ITKEEPA( J ) .AND. KEEPVAL )      ) THEN
                BUF1 = 'N'
                IF( KEEPVAL ) BUF1 = 'Y'
                WRITE( MESG,94010 ) 'WARNING: Keep status at '//
     &                 'line', IREC, 'different from previously '//
     &                 'set status of "'// BUF1// '" for data "' //
     &                 TRIM( ITNAMA( J ) ) // '".'
                CALL M3MSG2( MESG )

            END IF

C.............  Error msg for exact duplicate entries
C              : added to prevent doubling emssions by exact duplicate entries (ex: PM10)
            IF( AFLAG ) THEN    ! duplicate activity entries
            IF( CFLAG ) THEN    ! duplicate CAS#
            IF( NFLAG ) THEN    ! duplicate units 
            IF( SFLAG ) THEN    ! duplicate SPECAITE 4
            IF( VFLAG ) THEN    ! duplicate V or T
                WRITE( MESG,94010 ) 'ERROR: Duplicate entry of ' // 
     &                TRIM( ITNAMA( J ) ) // ' at line', IREC
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            END IF
            END IF
            END IF
            END IF
            END IF

C.............  Insert no-integrate entries (memory already increased with NNOI variable)
            IF ( ( ITVTSA(J) .EQ. 'V' .OR. ITVTSA(J) .EQ. 'T' ).AND.
     &             ITMSPC( J ) .AND. .NOT. ITEXPL( J ) ) THEN

                K = MXIDAT
                MXIDAT = MXIDAT + 1
                LOCATIDX  ( MXIDAT ) = MXIDAT
                INVPOSA ( MXIDAT ) = INVPOSA ( K ) + 0.1
                INVDCODA( MXIDAT ) = INVDCODA( K )
                INVSTATA( MXIDAT ) = INVSTATA( K )
                INVDVTSA( MXIDAT ) = INVDVTSA( K )
                INVDUNTA( MXIDAT ) = INVDUNTA( K )
                INVDDSCA( MXIDAT ) = INVDDSCA( K )
                LADDNOI ( MXIDAT ) = .TRUE.

            END IF

C.............  Set data previous name for next iteration
            LNAM = ITNAMA( J )

C.............  Initialize flags
            AFLAG = .FALSE.
            CFLAG = .FALSE.
            NFLAG = .FALSE.
            SFLAG = .FALSE.
            VFLAG = .FALSE.

        END DO       ! end loop to store unique data names

C.........  Abort if error found during read
        IF ( EFLAG ) THEN
            MESG = 'Problem with INVTABLE input file contents.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Resort inventory data names based on original position in inventory table
        CALL SORTR1( MXIDAT, LOCATIDX, INVPOSA )

        DO N = 1, MXIDAT
            J = LOCATIDX( N )

            IF( LADDNOI( J ) ) THEN

                INVDNAM( N ) = TRIM( ITNAMA( INT( INVPOSA( J ) ) ) ) // 
     &                         NOIEND
                INVDCOD( N ) = INVDCODA( J )
                INVSTAT( N ) = INVSTATA( J )
                INVDVTS( N ) = 'N' 
                INVDUNT( N ) = INVDUNTA( J )
                INVDDSC( N ) = 'No-integrate '// INVDDSCA( J )

            ELSE

                INVDNAM( N ) = ITNAMA( INT( INVPOSA( J ) ) )
                INVDCOD( N ) = INVDCODA( J )
                INVSTAT( N ) = INVSTATA( J )
                INVDVTS( N ) = INVDVTSA( J )
                INVDUNT( N ) = INVDUNTA( J )
                INVDDSC( N ) = INVDDSCA( J )

            END IF

C.............  Reset LOCATIDX for use in next sorting step
            LOCATIDX( N ) = N

        END DO

C.........  Sort data list by SPECIATE IDs
        CALL SORTI1( MXIDAT, LOCATIDX, INVDCOD )

C.........  Store SPECIATE IDs in sorted order
        DO N = 1, MXIDAT
            J = LOCATIDX( N )
            IDXCOD( N )  = J
            SORTCOD( N ) = INVDCOD( J )
        END DO

C.........  Rewind inventory table file
        IF ( FDEV .NE. 0 ) REWIND( FDEV )

C.........  Deallocate local memory
        DEALLOCATE( LOCATIDX, INVPOSA, INVDCODA, INVSTATA, INVDVTSA,
     &              INVDUNTA, INVDDSCA, LADDNOI )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDCODNAM
