
        PROGRAM ADJMOBIL

C***********************************************************************
C  program body starts at line 
C
C  DESCRIPTION:
C       Adjusts variables in the MOBL or MTMP file and outputs another file
C       NOTE- Initial implementation limited to MTMP
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; PROMPTMFILE
C
C  REVISION  HISTORY:
C       Prototype 1/99 by M Houyoux
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
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
C****************************************************************************

      IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'MBDIMS3.EXT'   !  mobil-source dimensioning parameters
        INCLUDE 'TMDIMS3.EXT'   !  temporal dimensioning parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2     CRLF
        INTEGER         FIND1
        INTEGER         FIND3
        INTEGER         FIND4
        INTEGER         GETFLINE
        LOGICAL         GETYN
        CHARACTER*10    HHMMSS
        INTEGER         INDEX1
        CHARACTER*14    MMDDYY
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         STR2INT
        REAL            STR2REAL
        INTEGER         WKDAY           !  day of week (1...7)

        EXTERNAL        CRLF, FIND1, FIND3, FIND4, GETFLINE, GETYN, 
     &                  HHMMSS, INDEX1, MMDDYY, PROMPTFFILE, 
     &                  PROMPTMFILE, STR2INT, STR2REAL, WKDAY

C...........   Local parameters
        INTEGER, PARAMETER :: MXTABL1 = NRCLAS * MXEMIS * NVTYPE
        INTEGER, PARAMETER :: MXTABL2 = NRCLAS * MXSTA * MXEMIS * NVTYPE
        INTEGER, PARAMETER :: MXTABL3 = NRCLAS * MXFIP * MXEMIS * NVTYPE

C...........   LOCAL VARIABLES and their descriptions:
C...........   Mobil Sources input and output arrays (one variable at a time)
        
        REAL            EMISV( NMSRC, NVTYPE )  !  input emissions.
        
C...........   Inventory arrays
        INTEGER         IFIP  ( NMSRC )   ! state and county FIPS code
        INTEGER         IRCLAS( NMSRC )   ! road class
        INTEGER         IFIDX ( NMSRC )   ! index to county report
        INTEGER         IRIDX ( NMSRC )   ! index to road class report
        INTEGER         ISIDX ( NMSRC )   ! index to state report
        INTEGER         SRCIDX( NMSRC, NVTYPE, MXEMIS )  ! table position

C...........   Unsorted allocatable adjustment factor table
        INTEGER, ALLOCATABLE :: ADJIDXA( : )
        INTEGER, ALLOCATABLE :: ADJFIPA( : )
        INTEGER, ALLOCATABLE :: ADJRDTA( : )
        INTEGER, ALLOCATABLE :: ADJVIDA( : )
        INTEGER, ALLOCATABLE :: ADJCODA( : )
        REAL   , ALLOCATABLE :: ADJFACA( : )

C...........   Sorted, grouped adjustment factors tables
        INTEGER     N1, N2, N3    ! table counters

        INTEGER     ADJID0( NVTYPE, MXEMIS )           ! fallback

        INTEGER     RDT1  ( MXTABL1 )                  ! roadclass only
        INTEGER     VID1  ( MXTABL1 )
        INTEGER     COD1  ( MXTABL1 )
        INTEGER     ADJID1( MXTABL1 )

        INTEGER     FIP2  ( MXTABL2 )                  ! state & roadclass
        INTEGER     RDT2  ( MXTABL2 )
        INTEGER     VID2  ( MXTABL2 )
        INTEGER     COD2  ( MXTABL2 )
        INTEGER     ADJID2( MXTABL2 )

        INTEGER     FIP3  ( MXTABL3 )                   ! county & roadclass
        INTEGER     RDT3  ( MXTABL3 )
        INTEGER     VID3  ( MXTABL3 )
        INTEGER     COD3  ( MXTABL3 )
        INTEGER     ADJID3( MXTABL3 )

C...........   Emissions report arrays
        INTEGER     FIPREP ( MXFIP )
        INTEGER     STAREP ( MXFIP )

        REAL        OLDEMIS ( NRCLAS, NVTYPE, MXEMIS )
        REAL        NEWEMIS ( NRCLAS, NVTYPE, MXEMIS )
        REAL        OLDEMSTA( MXSTA , NVTYPE, MXEMIS )
        REAL        NEWEMSTA( MXSTA , NVTYPE, MXEMIS )
        REAL        OLDEMFIP( MXFIP , NVTYPE, MXEMIS )
        REAL        NEWEMFIP( MXFIP , NVTYPE, MXEMIS )

C...........   Variable names arrays 
        INTEGER      NVAR
        LOGICAL      LADJUST( MXVARS3 )
        CHARACTER*16 VNAMES ( MXVARS3 )

C...........   Parsing buffer
        CHARACTER*20 SEGMENT( 5 )

C...........   File names and unit numbers
        INTEGER         LDEV
        INTEGER         MDEV               ! input factors unit number
        CHARACTER*16    ANAME   !  logical name for mobile-source input file
        CHARACTER*16    INAME   !  logical name for input inventory name
        CHARACTER*16    ONAME   !  logical name for output file name

C...........   Other local variables
        INTEGER         F, I, J, K, L, N, S, T, V

        INTEGER         FIP                ! tmp state & county code
        INTEGER         ICOD               ! tmp emis-type position
        INTEGER         IOS                ! i/o status
        INTEGER         JDATE, JTIME, TSTEP
        INTEGER         LCOD               ! previous emis-type pos'n
        INTEGER         LDATE, LTIME
        INTEGER         LFIP               ! previous state & county code
        INTEGER         LFIPREP            ! prev. reporting state & county code
        INTEGER         LRDT               ! previous road class code
        INTEGER         LSTAREP            ! prev. reporting state code
        INTEGER         LVID               ! previous vehicle type ID
        INTEGER         NADJ               ! no. entries in adjustments file
        INTEGER         NETYPE             ! no. emission types
        INTEGER         NFIPREP            ! no. counties for reporting
        INTEGER         NLINES             ! no. lines in input file
        INTEGER         NSTAREP            ! no. states for reporting
        INTEGER         NSTEPS
        INTEGER         RDT                ! temporary road class code
        INTEGER         VID                ! tmp vehicle type ID

        REAL            FAC  ! temporary adjustment factor

        LOGICAL      :: EFLAG = .FALSE.   ! true: error found

        CHARACTER*3     CVID    !  tmp char vehicle type ID code
        CHARACTER*6     CFIP    !  tmp char state & county code
        CHARACTER*6     CRDT    !  tmp char road class code
        CHARACTER*7     OUTNAM  !  contains default  logical output name
        CHARACTER*16    POL     !  tmp emission type name

        CHARACTER*300   LINE    !  file input line
        CHARACTER*300   MESG    !  scratch message area

        CHARACTER*16 :: PROGNAME = 'ADJMOBIL'

C***********************************************************************
C   begin body of program ADJMOBIL
        
        LDEV = INIT3()

        CALL INITEM( LDEV )
        
        WRITE( *,92000 )
     &  ' ',
     &  'Program ADJMOBIL to adjust mobile source emissions by state',
     &  'and county code, road class, vehicle type, and emission type.',
     &  'It uses an ASCII input file to assign adjustment factors and',
     &  'modify the MTMP file, and outputs a new MTMP file.',
     &  ' '
        WRITE( *,92000 ) 
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").  Use "NONE" as the',
     &  'name for the control-matrix file if you wish to omit control',
     &  'from the operations performed, or as the name of the output',
     &  'file if you want only to time program performance without the',
     &  'overhead of additional I/O.',
     &  ' '
        WRITE( *,92000 ) 
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry.  Default responses are indicated',
     &  'in brackets [LIKE THIS].',
     &  ' '
        
        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, 'Ending Program', 2 )
        END IF
 
C.......   Get file names; open input mobile source, mods table, and
C.......   output files.
        INAME = PROMPTMFILE( 
     &          'Enter logical name for MOBILE INVENTORY ' // 
     &          'file', FSREAD3, 'MOBL', PROGNAME )

        ANAME = PROMPTMFILE( 
     &          'Enter logical name for TIME-STEPPED ' // 
     &          'MOBILE EMIS file', FSREAD3, 'MTMP', PROGNAME )

C.......   Open table for modifying input 
        MDEV = PROMPTFFILE( 
     &           'Enter logical name for adjustments file',
     &           .TRUE., .TRUE., 'MAJFACS', PROGNAME )

C.......   Build description of output file, and optionally open it:

        IF ( .NOT. DESC3( INAME ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, 
     &              'Could not get description of file ' // INAME, 2 )

        ELSEIF( NROWS3D .NE. NMSRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  MOBILE INVEN file:', NROWS3D,
     &      'program:', NMSRC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
 
        END IF

        IF ( .NOT. DESC3( ANAME ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, 
     &              'Could not get description of file ' // ANAME, 2 )

        ELSEIF( NROWS3D .NE. NMSRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  MOBILE SOURCES file:', NROWS3D,
     &      'program:', NMSRC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
 
        END IF

        JDATE   = SDATE3D
        JTIME   = STIME3D
        TSTEP   = TSTEP3D
        NSTEPS  = MXREC3D
        NETYPE  = NVARS3D
        NVAR    = NVARS3D
        VNAMES  = VNAME3D   ! Array

C.........  Annual inventory VMT file as input
        IF ( TSTEP .EQ. 0 ) THEN
            MESG = 'Cannot run using MOBL file!'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Get the number of lines in the input file
        NLINES = GETFLINE( MDEV, 'Mobile adjustments' )

C.........  Allocate memory for unsorted records
        ALLOCATE( ADJIDXA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ADJIDXA', PROGNAME )
        ALLOCATE( ADJFIPA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ADJFIPA', PROGNAME )
        ALLOCATE( ADJRDTA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ADJRDTA', PROGNAME )
        ALLOCATE( ADJVIDA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ADJVIDA', PROGNAME )
        ALLOCATE( ADJCODA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ADJCODA', PROGNAME )
        ALLOCATE( ADJFACA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ADJFACA', PROGNAME )
 
C.......   Read adjustment factors
        N = 0
        DO I = 1, NLINES

            READ( MDEV, 93000, END=999, IOSTAT=IOS ) LINE

            IF( IOS .GT. 0 ) THEN
                WRITE( MESG,94010 ) 'I/O error', IOS,
     &                 'reading factors at line', I
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Skip blank or header lines
            IF( LINE .EQ. ' ' .OR. LINE( 1:1 ) .EQ. '#' ) CYCLE

C.............  Convert line to fields.
            CALL PARSLINE( LINE, 5, SEGMENT )

            CFIP =  ADJUSTL( SEGMENT( 1 ) )
            CRDT =  ADJUSTL( SEGMENT( 2 ) )
            CVID =  ADJUSTL( SEGMENT( 3 ) )
            POL  =  ADJUSTL( SEGMENT( 4 ) )
            FAC  = STR2REAL( SEGMENT( 5 ) )

            CALL PADZERO( CFIP )
            CALL PADZERO( CRDT )
            CALL PADZERO( CVID )

            FIP = STR2INT( CFIP )
            RDT = STR2INT( CRDT )
            VID = STR2INT( CVID )

C.............  Check that road class and vehicle type ID are valid
            J = FIND1( RDT, NRCLAS, MROADS3 ) 
            IF( RDT .NE. 0 .AND. J .LE. 0 ) THEN
                WRITE( MESG,94010 ) 'WARNING: Skipping entry for '//
     &                 'unrecognized road class', RDT, 'at line', I
                CALL M3MESG( MESG ) 
                CYCLE

            END IF

            IF( VID .NE. 0. .AND. VID .GT. NVTYPE ) THEN
                WRITE( MESG,94010 ) 'WARNING: Skipping entry for '//
     &                 'unrecognized vehicle type code', VID,
     &                 'at line', I
                CALL M3MESG( MESG ) 
                CYCLE

            END IF

C.............  Check that emission types are valid
            ICOD = INDEX1( POL, NVAR, VNAMES )

            IF( ICOD .LE. 0 ) THEN
                 L = LEN_TRIM( POL )
                 WRITE( MESG,94010 ) 'WARNING: Skipping entry for '//
     &                 'unrecognized emission type "' // POL( 1:L ) //
     &                 '" at line', I
                CALL M3MESG( MESG ) 
                CYCLE

            END IF

            N = N + 1

            ADJIDXA( N ) = N
            ADJFIPA( N ) = FIP
            ADJRDTA( N ) = RDT
            ADJVIDA( N ) = VID
            ADJCODA( N ) = ICOD
            ADJFACA( N ) = FAC

        END DO

        NADJ = N

        MESG = 'Processing mobile adjustments file...'
        CALL M3MSG2( MESG )

C.........  Sort entries
        CALL SORTI4( NADJ, ADJIDXA, ADJFIPA, ADJRDTA, ADJVIDA, ADJCODA )

C.........  Initialize category-specific cross-reference tables
        ADJID0 = IMISS3   ! array
        ADJID1 = IMISS3   ! array
        ADJID2 = IMISS3   ! array
        ADJID3 = IMISS3   ! array

C.........  Group cross-reference profiles by how specifically they are 
C           referenced
        N1 = 0
        N2 = 0
        N3 = 0
        LFIP = -9
        LRDT = -9
        LVID = -9
        LCOD = -9
        LFIPREP = -9
        LSTAREP = -9
        DO K = 1, NADJ

            J    = ADJIDXA( K )
            FIP  = ADJFIPA( J )
            RDT  = ADJRDTA( J )
            VID  = ADJVIDA( J )
            ICOD = ADJCODA( J )

C.............  Flag emission-types that are referenced in the input
C               file
            LADJUST( ICOD ) = .TRUE.

C.............  Give warning for duplicate records
            IF( FIP  .EQ. LFIP .AND.
     &          RDT  .EQ. LRDT .AND.
     &          VID  .EQ. LVID .AND.
     &          ICOD .EQ. LCOD       ) THEN

                WRITE( MESG,94010 ) 'WARNING: Duplicate record near ' //
     &                 'line', J, 'of adjustments file.'
                CALL M3MESG( MESG )

            ELSE
                LFIP = FIP
                LRDT = RDT
                LVID = VID
                LCOD = ICOD

            END IF

C.............  Assign group and store
            IF ( FIP .EQ. 0 ) THEN !  FIP-independent fallback profiles
 
                IF ( RDT .EQ. 0 ) THEN  !  ultimate fallback profiles
                    
        	    IF ( VID .EQ. 0 ) THEN
                	ADJID0( 1:NVTYPE, ICOD ) = J
        	    ELSE
                	ADJID0( VID, ICOD ) = J
        	    END IF

                ELSE                    !  RDT?-dependent-only profiles
 
                    IF( N1 .EQ. 0 ) THEN
                        N1 = N1 + 1

                    ELSE IF( RDT  .NE. RDT1( N1 ) .OR.
     &                       VID  .NE. VID1( N1 ) .OR.
     &                       ICOD .NE. COD1( N1 )      ) THEN
                        N1 = N1 + 1

                    END IF

                    IF ( N1 .LE. MXTABL1 ) THEN

                        RDT1  ( N1 ) = RDT
                        VID1  ( N1 ) = VID
                        COD1  ( N1 ) = ICOD
                        ADJID1( N1 ) = J

                    END IF

                END IF

            ELSE IF ( MOD( FIP,1000 ) .EQ. 0 ) THEN     ! State dependent
 
                FIP = FIP/1000
 
                IF( N2 .EQ. 0 ) THEN
                    N2 = N2 + 1

                ELSE IF( FIP  .NE. FIP2( N2 ) .OR.
     &                   RDT  .NE. RDT2( N2 ) .OR.
     &                   VID  .NE. VID2( N2 ) .OR.
     &                   ICOD .NE. COD2( N2 )      ) THEN
                    N2 = N2 + 1

                END IF

                IF ( N2 .LE. MXTABL2 ) THEN

                    FIP2  ( N2 ) = FIP
                    RDT2  ( N2 ) = RDT
                    VID2  ( N2 ) = VID
                    COD2  ( N2 ) = ICOD
                    ADJID2( N2 ) = J

                END IF

C.................  Update list of states for reporting
                IF( FIP .NE. LSTAREP ) THEN
                   NSTAREP = NSTAREP + 1
                   IF( NSTAREP .LE. MXSTA ) THEN
                       STAREP( NSTAREP ) = FIP
                   END IF
                END IF

                LSTAREP = FIP
 
            ELSE                              !  FIP-RDT dependent profiles
 
                IF( N3   .EQ. 0 ) THEN
                    N3 = N3 + 1

                ELSE IF( FIP  .NE. FIP3( N3 ) .OR.
     &                   RDT  .NE. RDT3( N3 ) .OR.
     &                   VID  .NE. VID3( N3 ) .OR.
     &                   ICOD .NE. COD3( N3 )      ) THEN
                     N3 = N3 + 1

                END IF

                IF ( N3 .LE. MXTABL3 ) THEN

                    FIP3  ( N3 ) = FIP
                    RDT3  ( N3 ) = RDT
                    VID3  ( N3 ) = VID
                    COD3  ( N3 ) = ICOD
                    ADJID3( N3 ) = J

                END IF

C.................  Update list of counties for reporting
                IF( FIP .NE. LFIPREP ) THEN
                   NFIPREP = NFIPREP + 1
                   IF( NFIPREP .LE. MXFIP ) THEN
                       FIPREP( NFIPREP ) = FIP
                   END IF
                END IF

                LFIPREP = FIP

            END IF

        END DO   ! End of x-ref parsing loop

C.........  Check table dimensions and print errors
        IF( N1 .GT. NRCLAS ) THEN

            EFLAG = .TRUE.
            WRITE( MESG,94010 )
     &             'Max roadclass adjustments NRCLAS=', NRCLAS,
     &             'exceeded in input file: count=', N1
            CALL M3MSG2( MESG )

        END IF

        IF( N2 .GT. NRCLAS ) THEN

            EFLAG = .TRUE.
            WRITE( MESG,94010 )
     &             'Max state-specific adjustments MXSTRC=', MXSTRC,
     &             'exceeded in input file: count=', N2
            CALL M3MSG2( MESG )

        END IF

        IF( N3 .GT. NRCLAS ) THEN

            EFLAG = .TRUE.
            WRITE( MESG,94010 )
     &             'Max county-specific adjustments MXFRC=', MXFRC,
     &             'exceeded in input file: count=', N3
            CALL M3MSG2( MESG )

        END IF

        IF ( EFLAG ) THEN
            MESG = 'Problem processing adjustments input file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Read source information
        IF( .NOT. READ3( INAME, 'IFIP', ALLAYS3, 
     &                   0, 0, IFIP )         ) THEN

            MESG = 'Could not read "IFIP" from ' // INAME
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

        END IF          !  if read3() succeeds or not

        IF( .NOT. READ3( INAME, 'IRCLAS', ALLAYS3, 
     &                   0, 0, IRCLAS )         ) THEN

            MESG = 'Could not read "IRCLAS" from ' // INAME
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

        END IF          !  if read3() succeeds or not

C.........  Initialize the indices to the reporting tables
        IRIDX = 0             !  array
        ISIDX = NSTAREP + 1   !  array
        IFIDX = NFIPREP + 1   !  array

C.........  Populate the indices to the reporting tables
        DO S = 1, NMSRC

C.............  Road class index
            J = FIND1( IRCLAS( S ), NRCLAS, MROADS3 )

            IF( J .LE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Road class', IRCLAS( S ),
     &                 'from inventory not found in valid list' //
     &                 CRLF() // '          for source', S
                CALL M3MESG( MESG )

            ELSE
                IRIDX( S ) = J

            END IF

C.............  County index
            FIP = IFIP( S )
            J = FIND1( FIP, NFIPREP, FIPREP )
            IF( J .GT. 0 ) THEN
                IFIDX( S ) = J
            END IF

C.............  State index
            FIP = FIP / 1000
            J = FIND1( FIP, NSTAREP, STAREP ) 
            IF( J .GT. 0 ) THEN
                ISIDX( S ) = J
            END IF

        END DO

C.........  Abort if error found
        IF( EFLAG ) THEN
            MESG = 'Problem with inventory file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Assign table position to each source, vehicle type, and 
C           emission type (or IMISS3 for no adjustments)
        DO V = 1, NVAR
            DO K = 1, NVTYPE 
        	DO S = 1, NMSRC

        	    FIP = IFIP  ( S )
        	    RDT = IRCLAS( S )

C.....................  FIP-RDT-VTYPE matches:
        	    F= FIND4(FIP, RDT, K, V, N3, FIP3, RDT3, VID3, COD3)
        	    IF ( F .GT. 0 ) THEN
                	SRCIDX( S, K, V ) = ADJID3( F )
                	CYCLE
        	    END IF
           
C.....................  FIP-RDT matches:
        	    F= FIND4(FIP, RDT, 0, V, N3, FIP3, RDT3, VID3, COD3)
        	    IF ( F .GT. 0 ) THEN
                	SRCIDX( S, K, V ) = ADJID3( F )
                	CYCLE
        	    END IF
           
C.....................  FIP-VTYPE matches:
        	    F= FIND4(FIP, 0, K, V, N3, FIP3, RDT3, VID3, COD3)
        	    IF ( F .GT. 0 ) THEN
                	SRCIDX( S, K, V ) = ADJID3( F )
                	CYCLE
        	    END IF
           
C.....................  FIP-only matches:
        	    F= FIND4(FIP, 0, 0, V, N3, FIP3, RDT3, VID3, COD3)
        	    IF ( F .GT. 0 ) THEN
                	SRCIDX( S, K, V ) = ADJID3( F )
                	CYCLE
        	    END IF
           
C.....................  State-RDT-VTYPE matches:
        	    FIP = FIP / 1000
        	    F= FIND4(FIP, RDT, K, V, N2, FIP2, RDT2, VID2, COD2)
        	    IF ( F .GT. 0 ) THEN
                	SRCIDX( S, K, V ) = ADJID2( F )
                	CYCLE
        	    END IF
           
C.....................  State-RDT matches:
        	    F= FIND4(FIP, RDT, 0, V, N2, FIP2, RDT2, VID2, COD2)
        	    IF ( F .GT. 0 ) THEN
                	SRCIDX( S, K, V ) = ADJID2( F )
                	CYCLE
        	    END IF
           
C.....................  State-VTYPE matches:
        	    F= FIND4(FIP, 0, K, V, N2, FIP2, RDT2, VID2, COD2)
        	    IF ( F .GT. 0 ) THEN
                	SRCIDX( S, K, V ) = ADJID2( F )
                	CYCLE
        	    END IF
           
C.....................  State-only matches:
        	    F= FIND4(FIP, 0, 0, V, N2, FIP2, RDT2, VID2, COD2)
        	    IF ( F .GT. 0 ) THEN
                	SRCIDX( S, K, V ) = ADJID2( F )
                	CYCLE
        	    END IF
           
C.....................  RDT-VTYPE matches:
        	    F = FIND3( RDT, K, V, N1, RDT1, VID1, COD1 )
        	    IF ( F .GT. 0 ) THEN
                	SRCIDX( S, K, V ) = ADJID1( F )
                	CYCLE
        	    END IF
           
C.....................  RDT-only matches:
        	    F = FIND3( RDT, 0, V, N1, RDT1, VID1, COD1 )
        	    IF ( F .GT. 0 ) THEN
                	SRCIDX( S, K, V ) = ADJID1( F )
                	CYCLE
        	    END IF
           
C.....................  Fallback:
        	    SRCIDX( S, K, V ) = ADJID0( K, V )

                END DO  ! end of source loop
            END DO      ! end of vehicle type loop
        END DO          ! end of emission type loop

C.........  Set up to open output file. Note that most header variables are
C           being set by the DESC3 call on the input file, earlier.
        L = LEN_TRIM( FDESC3D( 1 ) )
        FDESC3D( 1 ) = FDESC3D( 1 )( 1:L ) // ' adjusted by ADJMOBIL'
        IF ( TSTEP .EQ. 0 ) THEN
            OUTNAM = 'MOBLOUT'

            LADJUST = .FALSE.  ! array
  
            K = INDEX1( 'VMT', NVARS3D, VNAME3D )
            LADJUST( K ) = .TRUE.
            MESG = 'NOTE: Adjusting variable VMT'
            CALL M3MSG2( MESG )

        ELSE  
            OUTNAM  = 'MTMPOUT'
            NVARS3D = NVAR

C.............  Scan the input variable names and flag the ones that
C               contain the name of the emission type to adjust.
            DO V = 1, NVAR

                IF( LADJUST( V ) ) THEN

                    MESG = 'NOTE: Adjusting variable ' // VNAMES( V )
                    CALL M3MSG2( MESG )

                END IF
            END DO

        END IF

        ONAME = PROMPTMFILE( 'Enter logical name for output file',
     &                       FSUNKN3, OUTNAM, PROGNAME )

C.........  Initialize report arrays
        OLDEMIS = 0.    ! array
        NEWEMIS = 0.    ! array

C.........  Transform and write out mobile source emissions values:
        LDATE = 0
        DO T = 1, NSTEPS

C.............  If this is a new month, or new day, write message
            IF ( JDATE .GT. 0 .AND. LDATE. NE. JDATE ) THEN
 
                MESG = 'Processing ' //
     &                 DAYS( WKDAY( JDATE ) ) // MMDDYY( JDATE )
                CALL M3MSG2( MESG( 1:LEN_TRIM( MESG ) ) )
 
            END IF

C.............  Write to screen because WRITE3 only writes to LDEV
            IF( JDATE .GT. 0 ) WRITE( *, 93020 ) HHMMSS( JTIME )

C.............   Read in entire emissions field
            DO V = 1, NVAR

                IF( .NOT. READ3( ANAME, VNAMES( V ), ALLAYS3, 
     &                           JDATE, JTIME, EMISV )         ) THEN
                    L = LEN_TRIM( VNAMES( V ) )
                    MESG = 'Could not read "' // VNAMES( V )( 1:L ) //
     &                     '" from ' // ANAME
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                END IF          !  if read3() succeeds or not

C.................  Adjust emissions, if needed
                IF( LADJUST( V ) ) THEN

                    DO  K = 1, NVTYPE
                        DO  S = 1, NMSRC

C.............................  Retreive report indices
                            N = IRIDX( S )
                            F = IFIDX( S )
                            L = ISIDX( S )

C.............................  Store old emissions for report arrays
                            OLDEMIS ( N, K, V ) = OLDEMIS ( N, K, V ) +
     &                                            EMISV( S,K )

                            OLDEMSTA( L, K, V ) = OLDEMSTA( L, K, V ) +
     &                                            EMISV( S,K )

                            OLDEMFIP( F, K, V ) = OLDEMFIP( F, K, V ) +
     &                                            EMISV( S,K )

C.............................  Retrieve index to adjustments table and skip
C                               if no adjustments are to be made
                            J = SRCIDX( S, K, V )

                            IF( J .GT. 0 ) THEN
                                EMISV( S,K )= EMISV( S,K )* ADJFACA( J )
                            END IF

C.............................  Store new emissions for report arrays
                            NEWEMIS ( N, K, V ) = NEWEMIS ( N, K, V ) +
     &                                            EMISV( S,K )

                            NEWEMSTA( L, K, V ) = NEWEMSTA( L, K, V ) +
     &                                            EMISV( S,K )

                            NEWEMFIP( F, K, V ) = NEWEMFIP( F, K, V ) +
     &                                            EMISV( S,K )

                        END DO
                    END DO

                END IF

                IF ( .NOT. WRITE3( ONAME, VNAMES( V ),
     &                             JDATE, JTIME, EMISV  ) ) THEN

                    L = LEN_TRIM( VNAMES( V ) )
                    MESG = 'Could not write "' // VNAMES( V )( 1:L ) //
     &                     '" to ' // ONAME
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                END IF          ! if write3() failed

            END DO              ! end loop on input variables

            LDATE = JDATE
            LTIME = JTIME

            CALL NEXTIME( JDATE, JTIME, TSTEP )

        END DO          !  end loop on time steps

C.........  Report summary emissions to log file

        MESG = 'EMISSIONS REPORT'
        CALL M3MESG( MESG )

        DO V = 1, NVAR

            IF( LADJUST( V ) ) THEN

C.................  Write emission types header
                L = LEN_TRIM( VNAMES( V ) )
                MESG = '    Emissions type ' // 
     &                 VNAMES( V )( 1:L ) // ':'
                CALL M3MESG( MESG )

C.................  Write columns header for road class report
                WRITE( LDEV, 93030 )

C.................  Write emissions for road class report
                DO I = 1, NRCLAS 
                    DO K = 1, NVTYPE

                        WRITE( LDEV, 93040 ) MROADS3( I ), VTYPE3( K ),
     &                         OLDEMIS( I, K, V ), NEWEMIS( I, K, V )
                    END DO
                END DO 

C.................  Write columns header for state report
                WRITE( LDEV, 93050 )

C.................  Write emissions for state report
                DO I = 1, NSTAREP 
                    DO K = 1, NVTYPE

                        WRITE( LDEV, 93060 ) STAREP( I ), VTYPE3( K ),
     &                         OLDEMSTA( I, K, V ), NEWEMSTA( I, K, V )
                    END DO
                END DO 

C.................  Write columns header for county report
                WRITE( LDEV, 93070 )

C.................  Write emissions for county report
                DO I = 1, NFIPREP 
                    DO K = 1, NVTYPE

                        WRITE( LDEV, 93080 ) FIPREP( I ), VTYPE3( K ),
     &                         OLDEMFIP( I, K, V ), NEWEMFIP( I, K, V )
                    END DO
                END DO 

            END IF

        END DO

C.........  Normal completion
      
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C.........  Abnormal program ending
999     MESG = 'Unexpected end of file reached.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

C...........   Formatted file I/O formats............ 93xxx
 
93000   FORMAT( A )

93020   FORMAT( 8X, 'at time ', A8 )

93030   FORMAT( 7X, 'Road ', 1X, 'Vehicle', 4X, '   Old   ', 3X, 
     &          '   New   ', /, 7X, 'Class', 1X, ' Type  ', 4X, 
     &          'Emissions', 3X, 'Emissions' )

93040   FORMAT( 8X, I3.3, 3X, A5, 1X, F11.3, 1X, F11.3 )

93050   FORMAT( /, 6X, '     ', 1X, 'Vehicle', 4X, '   Old   ', 3X, 
     &          '   New   ', /, 6X, 'State', 1X, ' Type  ', 4X, 
     &          'Emissions', 3X, 'Emissions' )

93060   FORMAT( 8X, I2.2, 3X, A5, 1X, F11.3, 1X, F11.3 )

93070   FORMAT( /, 5X, '      ', 1X, 'Vehicle', 4X, '   Old   ', 3X, 
     &          '   New   ', /, 5X, 'County', 1X, ' Type  ', 4X, 
     &          'Emissions', 3X, 'Emissions' )

93080   FORMAT( 5X, I5.5, 3X, A5, 1X, F11.3, 1X, F11.3 )

C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT ( 10 ( A, :, I10, :, 2X ) )

        END PROGRAM ADJMOBIL

