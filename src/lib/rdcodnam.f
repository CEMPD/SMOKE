
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
C
C****************************************************************************
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

C.........  MODULES for public variables
C.........  This module contains the lists of unique inventory information
        USE MODLISTS

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL         CHKINT
        LOGICAL         CHKREAL
        CHARACTER*2     CRLF
        INTEGER         GETFLINE
        INTEGER         INDEX1
        INTEGER         STR2INT
        REAL            STR2REAL

        EXTERNAL CHKINT, CHKREAL, CRLF, GETFLINE, INDEX1, STR2INT, 
     &           STR2REAL

C...........   SUBROUTINE ARGUMENTS
        INTEGER , INTENT (IN) :: FDEV   ! iventory table unit no.

C...........   Parameters
        INTEGER, PARAMETER :: NFIELDS = 14  ! no. input fields
        INTEGER, PARAMETER :: FBEG( NFIELDS ) = 
     &                      ( / 1 , 13, 24, 30, 32, 34,
     &                          40, 42, 44, 46, 48, 52, 69, 110 / )
        INTEGER, PARAMETER :: FEND( NFIELDS ) = 
     &                      ( / 11, 22, 28, 30, 32, 38,
     &                          40, 42, 44, 46, 50, 67, 108, 153 / )

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: LOCIDX  ( : ) ! sorting index
        INTEGER, ALLOCATABLE :: INVPOSA ( : ) ! position in original list
        INTEGER, ALLOCATABLE :: INVDCODA( : ) ! 5-digit SAROAD code (if any)
        INTEGER, ALLOCATABLE :: INVSTATA( : ) ! Status (<0 activity; >0 pol)

        CHARACTER(LEN=1)      , ALLOCATABLE :: INVDVTSA( : ) ! V=VOC, T=TOG, N=not
        CHARACTER(LEN=IOULEN3), ALLOCATABLE :: INVDUNTA( : ) ! units for SMOKE intmdt inventory
        CHARACTER(LEN=DDSLEN3), ALLOCATABLE :: INVDDSCA( : ) ! inventory data description

C...........   Local arrays
        INTEGER      FLEN   ( NFIELDS )
        CHARACTER*64 SEGMENT( NFIELDS )

C...........   Other local variables
        INTEGER         J, L, N     !  counters and indices

        INTEGER         CNT     !  count of pollutants per CAS
        INTEGER         CNTK    !  count of kept pollutants per CAS
        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  record number
        INTEGER         NDAT    !  number of data values

        LOGICAL       :: EFLAG    = .FALSE.  ! true: error found
        LOGICAL          KEEPVAL             ! initialized keep status

        CHARACTER*1               BUF1    !  tmp 1-char buffer
        CHARACTER*16           :: FNAME = 'INVTABLE' ! Logical name of input file
        CHARACTER*256             MESG    !  message buffer
        CHARACTER*512          :: LINE    !  Input line
        CHARACTER(LEN=CASLEN3) :: PCAS    !  previous CAS code
        CHARACTER(LEN=IOVLEN3) :: LNAM    !  previous pollutant name
        CHARACTER(LEN=IOVLEN3) :: DNAM    !  tmp for pollutant or activity name
        CHARACTER(LEN=DDSLEN3) :: DSCVAL  !  initialized CAS description value

        CHARACTER*16 :: PROGNAME = 'RDCODNAM' ! program name

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

C.............  Skip comment lines
            IF( LINE( 1:1 ) .EQ. CINVHDR ) CYCLE

C.............  Check if line is a process/pollutant combination
            IF( INDEX( LINE( 1:16 ), ETJOIN ) > 0 ) THEN
                SEGMENT( 1 ) = ADJUSTL( LINE( 1:16 ) )

C.................  Skip CAS number, SAROAD code, and reactivity, then
C                   store remaining fields               
                DO N = 5, NFIELDS
                    SEGMENT( N ) = 
     &                  ADJUSTL( LINE( FBEG( N ):FEND( N ) ) )
                END DO
            ELSE

C.................  Parse the line into its sections based on file format def'n
                DO N = 1, NFIELDS

                    SEGMENT( N ) = 
     &                  ADJUSTL( LINE( FBEG( N ):FEND( N ) ) )

                END DO
            END IF

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
            IF( .NOT. CHKINT( SEGMENT( 3 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: SAROAD code is not an ' //
     &                 'integer at line', IREC
                CALL M3MSG2( MESG )
            END IF

            IF( .NOT. CHKINT( SEGMENT( 4 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Reactivity group code ' //
     &                 'is not an integer at line', IREC
                CALL M3MSG2( MESG )
            END IF

            IF( .NOT. CHKINT( SEGMENT( 11 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: NTI code is not an ' //
     &                 'integer at line', IREC
                CALL M3MSG2( MESG )
            END IF

C.............  Check that float fields are floats
            IF( .NOT. CHKREAL( SEGMENT( 6 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Factor value is not an '//
     &                 'floating point value at line', IREC
                CALL M3MSG2( MESG )
            END IF

C.............  Convert to upper case, where needed
            CALL UPCASE( SEGMENT( 1 )( 1:FLEN(1) ) )
            CALL UPCASE( SEGMENT( 2 )( 1:FLEN(2) ) )
            CALL UPCASE( SEGMENT( 5 )( 1:FLEN(5) ) )
            CALL UPCASE( SEGMENT( 7 )( 1:FLEN(7) ) )
            CALL UPCASE( SEGMENT( 8 )( 1:FLEN(8) ) )
            CALL UPCASE( SEGMENT( 9 )( 1:FLEN(9) ) )
            CALL UPCASE( SEGMENT( 10)( 1:FLEN(10) ) )

C.............  Correct unknown entries
            IF( SEGMENT( 7 ) .NE. 'V' .AND.
     &          SEGMENT( 7 ) .NE. 'T'       ) SEGMENT( 7 ) = 'N'
            IF( SEGMENT( 8 ) .NE. 'Y'       ) SEGMENT( 8 ) = 'N'
            IF( SEGMENT( 9 ) .NE. 'Y'       ) SEGMENT( 9 ) = 'N'
            IF( SEGMENT( 10) .NE. 'Y'       ) SEGMENT( 10) = 'N'

C.............  Error for blank units
            IF( SEGMENT( 12 ) .EQ. ' ' ) THEN

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
                ITCASA   ( NDAT ) = TRIM( SEGMENT( 2 ) )
                ITCODA   ( NDAT ) = STR2INT( SEGMENT( 3 ) )
                ITREAA   ( NDAT ) = STR2INT( SEGMENT( 4 ) )
                ITKEEPA  ( NDAT ) = ( SEGMENT( 5 ) .EQ. 'Y' )
                ITFACA   ( NDAT ) = STR2REAL( SEGMENT( 6 ) )
                ITVTSA   ( NDAT ) = TRIM( SEGMENT( 7 ) )
                ITMSPC   ( NDAT ) = ( SEGMENT( 8 ) .EQ. 'Y' )
                ITEXPL   ( NDAT ) = ( SEGMENT( 9 ) .EQ. 'Y' )
                ITNTIA   ( NDAT ) = STR2INT( SEGMENT( 11 ) )
                ITUNTA   ( NDAT ) = TRIM( SEGMENT( 12 ) )
                ITDSCA   ( NDAT ) = TRIM( SEGMENT( 13 ) )
                ITCASDSCA( NDAT ) = TRIM( SEGMENT( 14 ) )
                ITCASDNMA( NDAT ) = ITCASA( NDAT ) // ITNAMA( NDAT )

                IF( SEGMENT( 10 ) .EQ. 'Y' ) ITSTATA( NDAT ) = -1
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
            MESG = 'Problem with input file contents.'
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
        DO N = 1, NINVTBL
            J = ITIDXA2( N ) 
            IF( .NOT. ITKEEPA( J ) ) CYCLE
            
            IF( ITNAMA( J ) .NE. LNAM ) MXIDAT = MXIDAT + 1
            LNAM = ITNAMA( J )
        END DO

C.........  Allocate memory for local inventory pollutants/activities
        ALLOCATE( LOCIDX( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LOCIDX', PROGNAME )
        ALLOCATE( INVPOSA( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVPOSA', PROGNAME )
        ALLOCATE( INVDCODA( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDCODA', PROGNAME )
        ALLOCATE( INVSTATA( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVSTATA', PROGNAME )
        ALLOCATE( INVDVTSA( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDVTSA', PROGNAME )
        ALLOCATE( INVDUNTA( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDUNTA', PROGNAME )
        ALLOCATE( INVDDSCA( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDDSCA', PROGNAME )

        LOCIDX   = 0    ! array
        INVPOSA  = 0    ! array
        INVDCODA = 0    ! array
        INVSTATA = 0    ! array
        INVDVTSA = ' '  ! array
        INVDUNTA = ' '  ! array
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

C.........  Allocate memory for sorted SAROAD codes and index
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
            IF( INVPOSA( MXIDAT ) .EQ. 0 ) THEN
                 LOCIDX ( MXIDAT ) = MXIDAT
                 INVPOSA( MXIDAT ) = 
     &                    INDEX1( ITNAMA( J ), NINVTBL, ITNAMA )
            END IF
  
C.............  Store SAROAD code if not yet set, otherwise, check that its the same code
            IF( INVDCODA( MXIDAT ) .EQ. 0 ) THEN
                INVDCODA( MXIDAT ) = ITCODA( J )

            ELSE IF( ITCODA( J ) .NE. INVDCODA( MXIDAT ) ) THEN
                WRITE( MESG,94010 ) 'WARNING: Different SAROAD ' //
     &                 'code for the same data name at line', IREC,
     &                 CRLF()// BLANK10// 'Using code', 
     &                 INVDCODA( MXIDAT ), 'and ignoring code',
     &                 ITCODA( J )
                CALL M3MSG2( MESG )
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

C.............  Set data previous name for next iteration
            LNAM = ITNAMA( J )

        END DO       ! end loop to store unique data names

C.........  Resort inventory data names based on original position in inventory table
        CALL SORTI1( MXIDAT, LOCIDX, INVPOSA )

        DO N = 1, MXIDAT
            J = LOCIDX( N )
            INVDNAM( N ) = ITNAMA( INVPOSA( J ) )
            INVDCOD( N ) = INVDCODA( J )
            INVSTAT( N ) = INVSTATA( J )
            INVDVTS( N ) = INVDVTSA( J )
            INVDUNT( N ) = INVDUNTA( J )
            INVDDSC( N ) = INVDDSCA( J )

C.............  Reset LOCIDX for use in next sorting step
            LOCIDX( N ) = N

        END DO


C.........  Sort data list by SAROAD codes
        CALL SORTI1( MXIDAT, LOCIDX, INVDCOD )

C.........  Store SAROAD in sorted order
        DO N = 1, MXIDAT
            J = LOCIDX( N )
            IDXCOD( N )  = J
            SORTCOD( N ) = INVDCOD( J )
        END DO

C.........  Rewind inventory table file
        IF ( FDEV .NE. 0 ) REWIND( FDEV )

C.........  Deallocate local memory
        DEALLOCATE( LOCIDX, INVPOSA, INVDCODA, INVSTATA, INVDVTSA,
     &              INVDUNTA, INVDDSCA )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDCODNAM
