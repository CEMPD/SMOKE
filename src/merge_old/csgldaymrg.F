
C Version "@(#)$Id$ $Source$ $Date$ 

        PROGRAM CSGLDAYMRG

C***********************************************************************
C  program     CSGLDAYMRG body starts at line  262
C  subroutine  DAYFILE         starts at line 1442
C  subroutine  CTLSPC          starts at line 1555
C
C  DESCRIPTION:
C       Use gridding-transform, model-speciation, and control matrices 
C       from programs GRDAMAT, GRDPMAT, SPCAMAT, SPCPMAT, and optionally
C       CTLAMAT and CTLPMAT to merge, grid, speciate, and control 
C       source-level area and point source emissions files.
C
C  PRECONDITIONS REQUIRED:  
C       M3IO source-level source emissions in canonical order produced
C       by RAWAREA and RAWPOINT
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; GETYN, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  8/95 by CJC.
C
C       Version    2/96 by CJC for combinations of area, biogenic, 
C       mobile, and point
C
C       Version   12/96 by CJC:  day-specific input files for each
C       category of emissions, each with its own clock.
C
C       Version    1/97 by  SL:  Modified mobile speciation and gridding 
C       modules to make csgldaymrg.F consistent with csgmobil.F. 
C       
C       Version    2/97 by  SL:  A fully validated version that can run
C       each of the four processes (mb, pt, bo, ar) separately or in any 
C       combination of them.
C
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  0-based I/O API file desc. data structures.
        INCLUDE 'ARDIMS3.EXT'   !  area -source dimensioning parameters
        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'VGRD3.EXT'     !  (AQM) vertical dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'MBDIMS3.EXT'   !  mobile source parms and tables


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        LOGICAL         DSCGRID
        INTEGER         GETDATE
        INTEGER         GETNUM
        LOGICAL         GETYN
        INTEGER         INDEX1
        CHARACTER*16    PROMPTMFILE
        INTEGER         PROMPTFFILE
        INTEGER         TRIMLEN
        INTEGER         CSVNAME
        INTEGER         WKDAY
        CHARACTER*14    MMDDYY
        CHARACTER*10    HHMMSS


        
        EXTERNAL    DSCGRID, GETDATE, GETNUM, GETYN, INDEX1, 
     &              PROMPTMFILE, PROMPTFFILE, TRIMLEN, CSVNAME,
     &              WKDAY, MMDDYY, HHMMSS

C...........   PARAMETER

        CHARACTER*16    NONE
        PARAMETER     ( NONE = 'NONE' )

cmh temp
        INTEGER     NETYPE
c sl        PARAMETER (NETYPE=6)
        PARAMETER (NETYPE=8)
cmh temp done

        INTEGER     ANETYPE      ! No. of scenario-spec MV emis types
        INTEGER     ANMPOL       ! Model species derived from CSVNAME
        CHARACTER*16   AVNAME3D( MXVARS3 ) ! variable names (length MXDLEN3=80)

C...........   LOCAL VARIABLES and their descriptions:
C...........   Emissions Sources input and output arrays (one variable at a time)

        REAL            AREAV ( NASRC + 
     &                          NPSRC + NMSRC*NVTYPE*MXEMIS )  ! input area emis
        REAL            AREAC ( NASRC )  ! controlled speciated area emissions

        REAL            POINTV( NPSRC )  ! input point emissions.
        REAL            POINTC( NPSRC )  ! controlled speciated point emissions

        REAL            MOBILV( NMSRC, NVTYPE, MXEMIS ) !input mobile emissions.
        REAL            MOBILC( NMSRC )  ! controlled speciated mobile emissions

        REAL            EMISG ( NGRID, EMLAYS )  !  gridded output emissions

C...........   Gridding Matrices
      
        INTEGER         NA( NGRID )
        INTEGER         IA( NMATX )
        REAL            MA( NMATX )

        COMMON  / GAMAT / NA, IA, MA

        INTEGER         NM( NGRID )
        INTEGER         IM( NMATX )
        REAL            MM( NMATX )

        COMMON  / GMMAT / NM, IM, MM

        INTEGER         NP( NGRID )
        INTEGER         IP( NPSRC )

        COMMON  / GPMAT / NP, IP

C...........   Vertical layering fractions for point source plume rise
C...........   Note that LFRAC(K,L) and FRACL(L,K ) refer to point 
C...........   source IP(K).  (Note that transpose FRACL gives sequential
C...........   memory access in MRGLVEC() ).

        REAL            LFRAC( NPSRC, EMLAYS )
        REAL            FRACL( EMLAYS, NPSRC )

C...........   Control and projection matrices (re-used for both purposes)

        REAL            CA( NASRC, NIPOL )
        REAL            CP( NPSRC, NIPOL )

cmh temp
        INTEGER  NEMODE
        PARAMETER( NEMODE = 4 )
        CHARACTER*16 MEMODE( NEMODE )  ! Emissions modes for given scenario
        DATA MEMODE  / 'EXH', 'EVP', 'RNL', 'RST' /
cmh temp done

        REAL            CM( NMSRC, NVTYPE, NIPOL + NEMODE - 1 )
        CHARACTER*16    MCNAMES( NIPOL + NEMODE - 1 )
        CHARACTER*16    MPNAMES( NIPOL + NEMODE - 1 )

C...........   Speciation Matrices, and description of how it maps input 
C...........   variables to output variables:  SS( S, SPCREF( V,W ) )
C...........   is used to produce species # OUTREF in the output file
C...........   from species V in the input file
C...........   (NOTE that names in the speciation-matrix file are of the 
C...........   form "<INV-SPC-NAME>to<MODEL-SPC-NAME>")

        REAL            SA( NASRC, NMPOL )
        REAL            SP( NPSRC, NMPOL )
        REAL            SM( NMSRC, NVTYPE, MXMPOL*MXEPRC )
        INTEGER         MBLMODE( MXMPOL*MXEPRC )
        INTEGER         MBLIPOL( MXMPOL*MXEPRC )
        INTEGER         MBLSDEX( NMPOL )
        INTEGER         MBLSCNT( NMPOL )

        INTEGER         NSPCOUT( NIPOL )      !  # of active rows in SPCREF
        INTEGER         SPCREF ( NIPOL, NMPOL )
        INTEGER         NPPM( MXMPOL )        !  No. of spcs facs per model spcs
        
        CHARACTER*16    INVNAME( NMPOL )      !  inv-spc part of name
        CHARACTER*16    MSNAMES( MXMPOL*MXEPRC )!  mobile spc-names
        CHARACTER*16    SPCFNM( MXMPOL*MXEPRC ) !  Names of speciation factors
        CHARACTER*16    NMATCH                  !  Name for matching
        CHARACTER*16    EMTNAM( MXEMIS )        ! Names of emission types
        CHARACTER*16    LNAME            ! Var name from previous loop iteration

C...........   Gridding Matrix, Control Matrix, Speciation Matrix

        INTEGER         NX( NGRID )
        INTEGER         IX( NMATX )
        REAL            CX( NMATX )

C...........   Other local variables
        
        INTEGER         JDATE, JTIME, TSTEP   !  output-file        date&time
        INTEGER         NSTEPS                !  number of output time steps
        INTEGER         LDATE                 !  last output-step's date
        INTEGER         ADATE, ATIME          !  area-source   file date&time
        INTEGER         BDATE, BTIME          !  bio-source    file date&time
        INTEGER         MDATE, MTIME          !  mobile-source file date&time
        INTEGER         PDATE, PTIME          !  point-source  file date&time
        INTEGER         SIDX                  !  Speciation matrix index no.
        INTEGER         TDEV                  !  Speciation list file unit
        INTEGER         NSNAME                ! No. of speciation factor names
        INTEGER         ICNT                  ! Counter
        INTEGER         LDEV
        INTEGER         I, J, K, L, M, N, F, S, T, V
        INTEGER         KK, JJ, III
        
        LOGICAL         AFLAG   !  area     source emissions go into the merge
        LOGICAL         BFLAG   !  biogenic source emissions go into the merge
        LOGICAL         MFLAG   !  mobile   source emissions go into the merge
        LOGICAL         PFLAG   !  point    source emissions go into the merge
        LOGICAL         ACFLAG  !  use area   source control matrix
        LOGICAL         MCFLAG  !  use mobile source control matrix
        LOGICAL         PCFLAG  !  use point  source control matrix
        LOGICAL         APFLAG  !  use area   source projection matrix
        LOGICAL         MPFLAG  !  use mobile source projection matrix
        LOGICAL         PPFLAG  !  use point  source projection matrix
        LOGICAL         METFLAG !  check consistency of input files' met headers
        LOGICAL         BSPCFLAG( NMPOL )       !  is spc in bio input file?
        LOGICAL         OUTFLAG !  generating file output

        CHARACTER*16    ANAME   !  logical name for area-source       input file
        CHARACTER*16    ACNAME  !  logical name for control    matrix input file
        CHARACTER*16    AGNAME  !  logical name for grid       matrix input file
        CHARACTER*16    APNAME  !  logical name for projection matrix input file
        CHARACTER*16    ASNAME  !  logical name for speciation matrix input file

        CHARACTER*16    MCNAME  !  logical name for projection matrix input file
        CHARACTER*16    MGNAME  !  logical name for grid       matrix input file
        CHARACTER*16    MPNAME  !  logical name for projection matrix input file
        CHARACTER*16    MSNAME  !  logical name for speciation matrix input file
        CHARACTER*16    MTNAME  !  logical name for time-stepped mb   input file

        CHARACTER*16    PNAME   !  logical name for point-source      input file
        CHARACTER*16    PCNAME  !  logical name for control    matrix input file
        CHARACTER*16    PGNAME  !  logical name for grid       matrix input file
        CHARACTER*16    PLNAME  !  logical name for layer coefficient input file
        CHARACTER*16    PPNAME  !  logical name for projection matrix input file
        CHARACTER*16    PSNAME  !  logical name for speciation matrix input file

        CHARACTER*16    BTNAME  !  logical name for biogenic temporal input file
        CHARACTER*16    ENAME   !  logical name for emission         output file

        CHARACTER*16    SCRBUF, TMPBUF, ABUF, BBUF  !  scratch vbles for names
        CHARACTER*80    BMETHEAD( 2 )  ! met headers from BTMP
        CHARACTER*80    MMETHEAD( 2 )  ! met headers from MTMP
        CHARACTER*80    PMETHEAD( 2 )  ! met headers from PLAY

        CHARACTER*10    DAYS( 7 )
        DATA            DAYS
     &          /
     &          'Monday', 'Tuesday', 'Wednesday', 'Thursday',
     &          'Friday', 'Saturday', 'Sunday'
     &          /

        CHARACTER*256   MESG
 

C.......   EQUIVALENCES/OVERLAYS to reduce executable size:
C...........   Need only one control matrix at a time; don't need
C...........   control matrices at the same time as input emissions
C...........   or layer fractions.

        EQUIVALENCE( CA, CP )
        EQUIVALENCE( CA, CM )
        EQUIVALENCE( CA, AREAV )
        EQUIVALENCE( AREAV( NASRC+1       ), POINTV )
        EQUIVALENCE( AREAV( NASRC+NPSRC+1 ), MOBILV )

C...........   STATEMENT FUNCTION:  floating point "definitely unequal"
        
        LOGICAL         FLTERR
        REAL            P, Q

        FLTERR( P, Q ) = 
     &      ( (P - Q)**2  .GT.  1.0E-10*( P*P + Q*Q + 1.0E-5 ) )


C***********************************************************************
C   begin body of program CSGLDAYMRG
        
        LDEV = INIT3()
        
        WRITE( *,92000 ) ' ', ' ', ' ',
     &'Program CSGLDAYMRG to take combinations of the sorted source    ',
     &'level area, bio, mobile and point source emissions files and   ',
     &'gridding-transform, speciation-transform, point-source layer   ',
     &'fractions if appropriate, and optionally control and projection',
     &'matrices, and produce layered gridded model emissions.  The    ',
     &'output time step structure may be either hourly time stepped   ',
     &'for a time period with selected starting date, time, and length',
     &'or time-independent with a single output record.  For each     ',
     &'input emissions type and each simulation day, the program will ',
     &'prompt you for an input file name and the date-and-time which  ',
     &'corresponds to the current simulation date-and-time.           ',
     &' '
        WRITE( *,92000 ) 
     &'You will need to have set the logical names for the input and ',
     &'output files prior to program launch, using                   ',
     &'   "setenv <logicalname> <pathname>".                         ',
     &'You will be prompted for these logical names during program   ',
     &'execution.  Use "NONE" as the name for the control matrix     ',
     &'files if you wish to omit control from the operations         ',
     &'performed, or as the name of the output file if you want only ',
     &'to time program performance without the overhead of output I/O',
     &' '
        WRITE( *,92000 ) 
     &'You may use END_OF-FILE (control-D) to quit the program       ',
     &'during logical-name entry.  Default responses are indicated   ',
     &'in brackets [LIKE THIS].                                      ',
     &' '
        
        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'CSGLDAYMRG', 0, 0, 
     &                   'Ending program at user request', 0 )
        END IF


C.......   Set appropriate flags and parameters for this run:

        AFLAG = GETYN( 'Do you want HOURLY output?', .TRUE. )
        
        IF ( AFLAG ) THEN

            JDATE  = GETDATE( JDATE,
     &           'Enter simulation starting date (YYYYDDD)|(YYYYMMDD)' )

            JTIME   = GETNUM( 0,  235959,       0, 
     &                        'Enter starting time  (HHMMSS)' )
            TSTEP   = 10000
            NSTEPS  = GETNUM( 1, 24, 999999,
     &                        'Enter run duration    (hours)' )

        ELSE

            JDATE   = 0
            JTIME   = 0
            TSTEP   = 0
            NSTEPS  = 1

        END IF                  !  if hourly output desired, or not


        ADATE = 0
        ATIME = 0
        AFLAG = GETYN( 'Include AREA SOURCES in merge?', .TRUE. )

        BDATE = 0
        BTIME = 0
        BFLAG = GETYN( 'Include BIOGENIC SOURCES in merge?', .TRUE. )

        MDATE = 0
        MTIME = 0
        MFLAG = GETYN( 'Include MOBILE SOURCES in merge?', .TRUE. )

        PDATE = 0
        PTIME = 0
        PFLAG = GETYN( 'Include POINT SOURCES in merge?', .TRUE. )

        IF ( .NOT.( AFLAG .OR. BFLAG .OR. MFLAG .OR. PFLAG ) ) THEN
            CALL M3EXIT( 'CSGLDAYMRG', 0, 0, 'No input files', 2 )
        END IF

        METFLAG = GETYN( 'Check consistency of met headers?', .TRUE. )

C.......   As appropriate, open and check dimensionality of matrix files:
C...........   NOTE that CLTSPC() and DAYFILE() called later perform
C...........   the dimension checks for control and time-stepped emis
C...........   files, respectively.  DAYFILE() also _opens_ emis files.

C...........   Area source files:

        IF ( AFLAG ) THEN

            AGNAME = PROMPTMFILE( 
     &          'Enter name for  AREA GRIDDING MATRIX file',
     &          FSREAD3, 'AGMAT', 'CSGLDAYMRG' )

            ASNAME = PROMPTMFILE( 
     &          'Enter name for  AREA SPECIATION MATRIX file',
     &          FSREAD3, 'ASMAT', 'CSGLDAYMRG' )

            ACNAME = PROMPTMFILE( 
     &          'Enter name for  AREA CONTROL MATRIX file ' //
     &          'or "NONE"',
     &          FSREAD3, 'ACMAT', 'CSGLDAYMRG' )

            APNAME = PROMPTMFILE( 
     &          'Enter name for  AREA PROJECTION MATRIX file ' 
     &          // 'or "NONE"',
     &          FSREAD3, 'APMAT', 'CSGLDAYMRG' )

            ACFLAG = ( ACNAME .NE. NONE )
            APFLAG = ( APNAME .NE. NONE )

            IF ( .NOT. DESC3( AGNAME ) ) THEN
                MESG = 'Could not get description of file "' //
     &                 AGNAME( 1:TRIMLEN( AGNAME ) ) // '"'
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )
            ELSE IF ( NROWS3D .NE. NGRID ) THEN
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  AGMAT file:', NROWS3D,
     &          'program:', NGRID
                CALL M3EXIT( 'CSGMERGE', 0, 0, MESG, 2 )
            ELSE IF ( NCOLS3D .NE. NMATX ) THEN
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  NMATX in AGMAT file:', NCOLS3D,
     &          'and program:', NMATX, 'do not match.'
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )
            ELSE IF ( NTHIK3D .NE. NASRC ) THEN
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  NASRC in AGMAT file:', NTHIK3D,
     &          'and program:', NASRC, 'do not match.'
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )
            END IF 	!  if agmat file's dimensionality matches, or not

            IF ( .NOT. DESC3( ASNAME ) ) THEN
                MESG = 'Could not get description of file "' // 
     &                 ASNAME( 1:TRIMLEN( ASNAME ) ) // '"'
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )
            ELSE IF ( NROWS3D .NE. NASRC ) THEN
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  ASMAT file:', NROWS3D,
     &          'program:', NASRC
                CALL M3EXIT( 'CSGMERGE', 0, 0, MESG, 2 )
            END IF      ! if asmat file's dimensionality matches, or not

        ELSE	!  area not used.  set area-related flags to false

            ACFLAG = .FALSE.
            APFLAG = .FALSE.

        END IF          !  if area emissions used, or not


C.......   Mobile source matrix files:

        IF ( MFLAG ) THEN

            MGNAME = PROMPTMFILE( 
     &          'Enter name for MOBILE GRIDDING MATRIX file',
     &          FSREAD3, 'MGMAT', 'CSGLDAYMRG' )

            MSNAME = PROMPTMFILE( 
     &          'Enter name for MOBILE SPECIATION MATRIX file',
     &          FSREAD3, 'MSMAT', 'CSGLDAYMRG' )

            TDEV = PROMPTFFILE(
     &           'Enter logical name for SPECIATION LIST file',
     &           .TRUE., .TRUE., 'SPCS', 'CSGLDAYMRG' )

            MTNAME = PROMPTMFILE(
     &          'Enter logical name for TIME-STEPPED MOBILE EMIS file',
     &          FSREAD3, 'MTMP', 'CSGLDAYMRG' )
            IF ( .NOT. DESC3( MTNAME ) ) THEN
                MESG =  'Could not get description for file "' //
     &                  MTNAME( 1:TRIMLEN( MTNAME ) ) // '"'
                CALL M3EXIT( 'CSGLDAYMRG', JDATE, JTIME, MESG, 2 )
            END IF          !  if DESC3() failed
            IF ( NROWS3D .NE. NMSRC ) THEN
                WRITE( MESG, 94010 )
     &              'Dimension mismatch in file "',
     &              MTNAME( 1:TRIMLEN( MTNAME ) ), '" -- file:', 
     &              NROWS3D, 'caller:', NMSRC
                CALL M3EXIT( 'CSGLDAYMRG', JDATE, JTIME, MESG, 2 )
            END IF          !  if dimension mismatch
            MMETHEAD( 1 ) = FDESC3D( 2 )
            MMETHEAD( 2 ) = FDESC3D( 3 )

            MCNAME = PROMPTMFILE( 
     &          'Enter name for MOBILE CONTROL MATRIX file ' //
     &          'or "NONE"',
     &          FSREAD3, 'MCMAT', 'CSGLDAYMRG' )

            MPNAME = PROMPTMFILE( 
     &          'Enter name for MOBILE PROJECTION MATRIX file ' 
     &          // 'or "NONE"',
     &          FSREAD3, 'MPMAT', 'CSGLDAYMRG' )

            MCFLAG = ( MCNAME .NE. NONE )
            MPFLAG = ( MPNAME .NE. NONE )

            IF ( .NOT. DESC3( MGNAME ) ) THEN
                MESG = 'Could not get description of file "' //
     &                 MGNAME( 1:TRIMLEN( MGNAME ) ) // '"'
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )
            ELSE IF ( NROWS3D .NE. NGRID ) THEN
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  MGMAT file:', NROWS3D,
     &          'program:', NGRID
                CALL M3EXIT( 'CSGMERGE', 0, 0, MESG, 2 )
            ELSE IF ( NCOLS3D .NE. NMATX ) THEN
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  NMATX in MGMAT file:', NCOLS3D,
     &          'and program:', NMATX, 'do not match.'
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )
            ELSE IF ( NTHIK3D .NE. NMSRC ) THEN
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  NMSRC in AGMAT file:', NTHIK3D,
     &          'and program:', NMSRC, 'do not match.'
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )
            END IF 	!  if mgmat file's dimensionality matches, or not

            IF ( .NOT. DESC3( MSNAME ) ) THEN
                MESG = 'Could not get description of file ' // MSNAME
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )
            ELSE IF ( NROWS3D .NE. NMSRC ) THEN
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  MSMAT file:', NROWS3D,
     &          'program:', NMSRC
                CALL M3EXIT( 'CSGMERGE', 0, 0, MESG, 2 )
            END IF

            IF( MCFLAG ) THEN
                IF ( .NOT. DESC3( MCNAME ) ) THEN
                    MESG = 'Could not get description of file '//MCNAME
                    CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )
                ELSE IF ( NROWS3D .NE. NMSRC ) THEN
                    WRITE( MESG, 94010 )
     &              'Dimension mismatch.  MCMAT file:', NROWS3D,
     &              'program:', NMSRC
                    CALL M3EXIT( 'CSGMERGE', 0, 0, MESG, 2 )
                END IF
                DO  5  I = 1, NIPOL + NEMODE - 1

5               CONTINUE
            END IF

            IF( MPFLAG ) THEN
                IF ( .NOT. DESC3( MPNAME ) ) THEN
                    MESG = 'Could not get description of file '//MPNAME
                    CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )
                ELSE IF ( NROWS3D .NE. NMSRC ) THEN
                    WRITE( MESG, 94010 )
     &              'Dimension mismatch.  MPMAT file:', NROWS3D,
     &              'program:', NMSRC
                    CALL M3EXIT( 'CSGMERGE', 0, 0, MESG, 2 )
                END IF
                DO  6  I = 1, NIPOL + NEMODE - 1
                    MCNAMES( I ) = VNAME3D( I )( 6:16 )
6               CONTINUE
            END IF

        ELSE		!  no mobile.  set mobile-related flags to false

            MCFLAG = .FALSE.
            MPFLAG = .FALSE.

        END IF          !  if mobile emissions used, or not


C.......   Point source matrix and layer fraction files:

        IF ( PFLAG ) THEN

            PGNAME = PROMPTMFILE( 
     &          'Enter name for POINT GRIDDING MATRIX file',
     &          FSREAD3, 'PGMAT', 'CSGLDAYMRG' )
     
            PSNAME = PROMPTMFILE( 
     &          'Enter name for POINT SPECIATION MATRIX file',
     &          FSREAD3, 'PSMAT', 'CSGLDAYMRG' )
     
            PLNAME = PROMPTMFILE( 
     &              'Enter name for POINT LAYER COEFF file',
     &              FSREAD3, 'PLAY', 'CSGLDAYMRG' )

            PCNAME = PROMPTMFILE( 
     &          'Enter name for POINT CONTROL MATRIX file ' //
     &          'or "NONE"',
     &          FSREAD3, 'PCMAT', 'CSGLDAYMRG' )

            PPNAME = PROMPTMFILE( 
     &      'Enter name for POINT PROJECTION MATRIX file ' //
     &      'or "NONE"',
     &      FSREAD3, 'PPMAT', 'CSGLDAYMRG' )

            PCFLAG = ( PCNAME .NE. NONE )
            PPFLAG = ( PPNAME .NE. NONE )

            IF ( .NOT. DESC3( PLNAME ) ) THEN
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0,
     &                  'Could not get description of file "'
     &                  // PLNAME( 1:TRIMLEN( PLNAME ) ) // '"', 2 )
            ELSE IF ( NROWS3D .NE. NPSRC ) THEN
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  LAYFRAC file:', NROWS3D,
     &          'program:', NPSRC
                CALL M3EXIT( 'CSGMERGE', 0, 0, MESG, 2 )
            ELSE IF ( NLAYS3D .NE. EMLAYS ) THEN
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  NLAYS in PLAY file:', NLAYS3D,
     &          'and program:', EMLAYS, 'do not match.'
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )
            END IF !  if file's dimensionality matches, or niot
            PMETHEAD( 1 ) = FDESC3D( 2 ) ! save met header to an array
            PMETHEAD( 2 ) = FDESC3D( 3 ) 

            DO  11  L = 0, EMLAYS
                IF ( FLTERR( SIGMAF( L ), VGLVS3D( L ) ) ) THEN
                    WRITE( MESG,94020 )
     &              'Layer', L, 'mismatch.  ' //
     &              'in file:',      VGLVS3D( L ),
     &              'in program:', SIGMAF ( L )
                    CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )
                END IF
11          CONTINUE

            IF ( .NOT. DESC3( PGNAME ) ) THEN
                MESG = 'Could not get description of file "' //
     &                 PGNAME( 1:TRIMLEN( PGNAME ) ) // '"'
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )
            ELSE IF ( NROWS3D .NE. NGRID ) THEN
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  NGRID in PGMAT file:', NROWS3D,
     &          'program:', NGRID
                CALL M3EXIT( 'CSGMERGE', 0, 0, MESG, 2 )
            ELSE IF ( NCOLS3D .NE. NPSRC ) THEN
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  NPSRC in PGMAT file:', NCOLS3D,
     &          'and program:', NPSRC, 'do not match.'
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )
            END IF !  if pgmat file's dimensionality matches, or not

            IF ( .NOT. DESC3( PSNAME ) ) THEN
                MESG = 'Could not get description of file ' // PSNAME
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )
            ELSE IF ( NROWS3D .NE. NPSRC ) THEN
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  PSMAT file:', NROWS3D,
     &          'program:', NPSRC
                CALL M3EXIT( 'CSGMERGE', 0, 0, MESG, 2 )
            END IF !  if psmat file's dimensionality matches, or not

        ELSE            !  no point. set default layer parameters, pt-flags

            PCFLAG = .FALSE.
            PPFLAG = .FALSE.

        END IF          !  if point sources used, or not


C.......   Biogenic source file

C........  Adding DESC3(BTMP) to avoid failure while opening output file  
C........  when bio is running alone 

        IF ( BFLAG .AND. (.NOT. (MFLAG.OR.PFLAG.OR.AFLAG) ) ) THEN

            BTNAME = PROMPTMFILE(
     &         'Enter logical name for TIME-STEPPED BIOGENIC EMIS file',
     &          FSREAD3, 'BTMP', 'CSGLDAYMRG' )

            IF ( .NOT. DESC3( BTNAME ) ) THEN
                MESG = 'Could not get description of file ' // BTNAME
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )
            END IF

            IF ( NROWS3D .NE. NROWS ) THEN
                WRITE( MESG, 94010 )
     &              'Dimension mismatch in file "',
     &              BTNAME( 1:TRIMLEN( BTNAME ) ), '" -- file:', 
     &              NROWS3D, 'caller:', NROWS
                CALL M3EXIT( 'CSGLDAYMRG', JDATE, JTIME, MESG, 2 )
            END IF          !  if dimension mismatch

            BMETHEAD( 1 ) = FDESC3D( 2 ) 
            BMETHEAD( 2 ) = FDESC3D( 3 ) 
            DO  101  V = 1, NMPOL	! Just set vname3d(*)
                BSPCFLAG( V ) = 
     &          ( INDEX1( EMNAM(V), NVARS3D, VNAME3D ) .GT. 0 )
101         CONTINUE

        END IF          !  if biogenic sources used, or not

C...........   Set up header for output file

        IF ( .NOT. DSCGRID( GRDNM, SCRBUF,
     &                      GDTYP3D, P_ALP3D, P_BET3D, P_GAM3D, 
     &                      XCENT3D, YCENT3D, XORIG3D, YORIG3D, 
     &                      XCELL3D, YCELL3D, 
     &                      NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            CALL M3EXIT( 'CSGLDAYMRG', 0, 0,
     &              'Could not get description of grid ' // GRDNM, 2 )

        END IF

        VGTYP3D = VGTYP
        VGTOP3D = VGTOP

        DO  12  L = 0, EMLAYS
            VGLVS3D( L ) = SIGMAF( L )
12      CONTINUE

        IF ( .NOT. PFLAG ) THEN
           NLAYS3D = 1
        ELSE
           NLAYS3D = EMLAYS
        END IF

        NVARS3D = NMPOL
        FTYPE3D = GRDDED3
        TSTEP3D = TSTEP
        SDATE3D = JDATE
        STIME3D = JTIME
        GDNAM3D = GRDNM
        
        IF ( TSTEP .EQ. 0 ) THEN
            MESG = 'gridded speciated emissions'
            ABUF = 'moles'
        ELSE
            MESG = 'hourly gridded speciated emissions'
            ABUF = 'moles/hour'
        END IF

        DO  22 V = 1, NMPOL
            VNAME3D( V ) = EMNAM( V )
            UNITS3D( V ) = ABUF
            VDESC3D( V ) = MESG
            VTYPE3D( V ) = M3REAL
22      CONTINUE
        FDESC3D( 1 ) = 'area source speciated gridded emissions.'
                 
        IF ( ACFLAG .OR. MCFLAG .OR. PCFLAG ) THEN
            ENAME = 'LGTSC'
        ELSE
            ENAME = 'LGTS'
        END IF

        ENAME = PROMPTMFILE(  
     &          'Enter name for OUTPUT EMIS file, or "NONE"',
     &          FSUNKN3, ENAME, 'CSGLDAYMRG' )

        OUTFLAG = ( ENAME .NE. NONE )


C...........   Build description of how to use speciation matrices:

        DO  31  V = 1, NIPOL
            NSPCOUT( V ) = 0
31      CONTINUE

        DO  32  V = 1, NMPOL
            MBLSDEX( V ) = NMPOL*NMPOL*NEMODE   !  "huge" initialization
            MBLSCNT( V ) = 0
32      CONTINUE

        DO  33  V = 1, NMPOL

            K = INDEX( SFNAM( V ), '_to_' )
            IF( K .LT. 1 ) THEN
                CALL M3EXIT( 'CSGAREA', 0, 0,
     &              'Bad conversion factor name "' // 
     &              SFNAM( V )( 1:TRIMLEN( SFNAM( V ) ) ) //
     &              '" in ' // ASNAME, 2 )
            END IF
            INVNAME( V ) = SFNAM( V )( 1 : K-1 )

            I = INDEX1( INVNAME( V ), NIPOL, EINAM )
            IF ( I .LT. 0 ) THEN
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, 
     &                       'Bad input species "' // 
     &                       INVNAME( V )( 1:TRIMLEN( INVNAME( V ))) //
     &                       '" in ' // ASNAME, 2 )
            END IF

            SCRBUF = SFNAM( V )( K+4:16 )
            J = INDEX1( SCRBUF, NMPOL, EMNAM )
            IF ( J .LT. 0 ) THEN
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, 
     &                       'Bad output species "' // 
     &                       SCRBUF( 1:TRIMLEN( SCRBUF ) ) // 
     &                       '" in ' // ASNAME, 2 )
            END IF

            F = NSPCOUT( I ) + 1
            SPCREF ( I,F ) = J
            NSPCOUT( I )   = F

33      CONTINUE


C.......   Read relevant gridding and speciation matrices:

        IF ( AFLAG ) THEN

            WRITE( *,92000 ) 
     &      ' ', 'Reading area transform matrices.', ' '

            IF ( .NOT. READ3( AGNAME, 'ALL', 1, 0, 0, NA ) ) THEN

                MESG = 'Could not read gridding matrix from file "' //
     &                 AGNAME( 1 : TRIMLEN( AGNAME ) ) // '".'
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )

            END IF      !  if read3() failed for gridding matrix

            IF ( .NOT. READ3( ASNAME, 'ALL', 1, 0, 0, SA ) ) THEN

                MESG = 'Could not read SPECIATION matrix from file "'//
     &                 ASNAME( 1 : TRIMLEN( ASNAME ) ) // '".'
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )

            END IF      ! if desc3(), dim check, or read3() failed for spc mtx

        END IF      !  if aflag


        IF ( MFLAG ) THEN

            WRITE( *,92000 ) 
     &      ' ', 'Reading mobile transform matrices.', ' '

            IF ( .NOT. READ3( MGNAME, 'ALL', 1, 0, 0, NM ) ) THEN

                MESG = 'Could not read gridding matrix from file "' //
     &                 MGNAME( 1 : TRIMLEN( MGNAME ) ) // '".'
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )

            END IF      !  if read3() failed for gridding matrix

            IF ( .NOT. READ3( MSNAME, 'ALL', ALLAYS3, 0, 0, SM ) ) THEN

                MESG = 'Could not read speciation matrix from file "'//
     &                 MSNAME( 1 : TRIMLEN( MSNAME ) ) // '".'
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )

            END IF      !  if read3() failed for speciation matrix

C...........   Parse names of the variables in msmat:

C.......   Start description of output file with ANAME and grid description

        IF ( .NOT. DESC3( MTNAME ) ) THEN
            CALL M3EXIT( 'CSGMOBIL', 0, 0,
     &              'Could not get description of file ' // MTNAME, 2 )
        END IF

C.......   Get speciation naming information from list file TDEV
C..........  Get number of scenario-specific MV emission types
C..........     ( Use VNAME3D as dummy so don't have to declare new array)

        ANETYPE  = CSVNAME( TDEV, .TRUE., .FALSE., .TRUE., EMTNAM )
        NSNAME  = CSVNAME( TDEV, .TRUE., .TRUE. , .TRUE., SPCFNM )

C..........  Get number and names of scenario-specific model species

        ANMPOL   = CSVNAME(TDEV, .FALSE., .TRUE. , .FALSE., AVNAME3D)

C.........  Initialize "last name" holder "LNAME"
        L     = TRIMLEN( SPCFNM( 1 ) )
        J     = INDEX  ( SPCFNM( 1 )(   1:L ), '_' )
        K     = INDEX  ( SPCFNM( 1 )( J+1:L ), '_' ) + J
        LNAME = SPCFNM( 1 )( K+1:L )

        V     = 0
        ICNT  = 0
C.........  Count tail part of speciation factor variable names
        DO 42 I = 1, NSNAME

C.............  Find second "_" in variable names
            L   = TRIMLEN( SPCFNM( I ) )
            J   = INDEX  ( SPCFNM( I )(   1:L ), '_' )
            K   = INDEX  ( SPCFNM( I )( J+1:L ), '_' ) + J

C.............  If current name and last name don't match or on last
            IF ( SPCFNM( I )( K+1:L ) .NE. LNAME ) THEN

C.................  Store count of names
                V = V + 1
                NPPM( V ) = ICNT

C.................  Reset last name and counter
                LNAME = SPCFNM( I )( K+1:L )
                ICNT  = 1

            ELSE
                ICNT = ICNT + 1

            ENDIF

   42   CONTINUE  ! End counting model species in speciation factor names

C.........  Store last counter value
        V = V + 1
        NPPM( V ) = ICNT

        END IF      !  if mflag


        IF ( PFLAG ) THEN

            WRITE( *,92000 ) 
     &      ' ', 'Reading point transform matrices.', ' '

            IF ( .NOT. READ3( PGNAME, 'ALL', 1, 0, 0, NP ) ) THEN

                MESG = 'Could not read gridding matrix from file "' //
     &                 PGNAME( 1 : TRIMLEN( PGNAME ) ) // '".'
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )

            END IF      !  if read3() failed for gridding matrix

            IF ( .NOT. READ3( PSNAME, 'ALL', 1, 0, 0, SP ) ) THEN

                MESG = 'Could not read SPECIATION matrix from file "'//
     &                 PSNAME( 1 : TRIMLEN( PSNAME ) ) // '".'
                CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )

            END IF      !  if desc3(), dim check, or read3() failed for spc mtx

        END IF      !  if pflag


C...........   Apply control and/or projection matrices.
C.......   Area source controls and projections:

        IF( ACFLAG ) THEN
            CALL CTLSPC( ACNAME, NASRC, SA, CA, INVNAME )
        END IF          !  if acflag:  control area


        IF( APFLAG ) THEN
            CALL CTLSPC( APNAME, NASRC, SA, CA, INVNAME )
        END IF          !  if acflag:  control area


C.......   Point source controls and projections:

        IF( PCFLAG ) THEN
            CALL CTLSPC( PCNAME, NPSRC, SP, CP, INVNAME )
        END IF          !  if pcflag:  control point sources

        IF( PPFLAG ) THEN
            CALL CTLSPC( PPNAME, NPSRC, SP, CP, INVNAME )
        END IF          !  if ppflag:  project point sources


C.......   Mobile source controls and projections:

        IF( MCFLAG ) THEN

C...............   Parse the variable-names from the speciation matrix and 
C...............   map to inventory species and mobile emissions modes:

            DO  53  I = 1, MXMPOL*MXEPRC

                L = MBLIPOL( I )
                M = MBLMODE( I )

                TMPBUF = 'CTL_' // MEMODE( M )( 1:3 ) // 
     &                   '_' // EINAM( L )
                M = INDEX1( TMPBUF, NIPOL + NEMODE - 1, MCNAMES )
                IF ( M .LE. 0 ) THEN                             
                    MESG = 'Missing control variable "' //
     &                     TMPBUF( 1:TRIMLEN( TMPBUF ) ) //
     &                     '" FOR spc matrix vble name "' //
     &                     VNAME3D( I )( 1:TRIMLEN( VNAME3D( I ) ) ) //
     &                     '"'
                    CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )       
                END IF

                DO  52  L = 1, NVTYPE
                DO  51  S = 1, NMSRC
                    SM( S,L,I ) = SM( S,L,I ) * CM( S,L,M )  
51              CONTINUE
52              CONTINUE

53          CONTINUE

        END IF          !  if mcflag:  control mobile


        IF( MPFLAG ) THEN

C...............   Parse the variable-names from the speciation matrix and 
C...............   map to inventory species and mobile emissions modes:

            DO  57  I = 1, MXMPOL*MXEPRC

                L = MBLIPOL( I )
                M = MBLMODE( I )

                TMPBUF = 'PRJ_' // MEMODE( M )( 1:3 ) // 
     &                   '_' // EINAM( L )
                M = INDEX1( TMPBUF, NIPOL + NEMODE - 1, MPNAMES )
                IF ( M .LE. 0 ) THEN                             
                    MESG = 'Missing projection variable "' //
     &                     TMPBUF( 1:TRIMLEN( TMPBUF ) ) //
     &                     '" FOR spc matrix vble name "' //
     &                     VNAME3D( I )( 1:TRIMLEN( VNAME3D( I ) ) ) //
     &                     '"'
                    CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )       
                END IF

                DO  56  L = 1, NVTYPE
                DO  55  S = 1, NMSRC
                    SM( S,L,I ) = SM( S,L,I ) * CM( S,L,M )
55              CONTINUE
56              CONTINUE

57          CONTINUE

        END IF          !  if mpflag:  project mobile



C.......   Transform and write out area source emissions values:
C...........   Use one time-stepped biogenics file and one time-stepped
C...........   mobile file for this time period.  Use a (possibly-new)
C...........   area and point soruce time-stepped file for each new day
C...........   encountered in the run.
C...........   Note:  DAYFILE() prompts the user for a file name,
C...........   opens the file and checks its dimensionality against
C...........   what is expected, and then prompts the user for the 
C...........   file-date&time to match with the current simulation
C...........   date&time.

        IF ( OUTFLAG ) THEN
            MESG = 'Calculating and writing speciated gridded emissions'
        ELSE 
             MESG = 'Calculating speciated gridded emissions...'
        END IF
        CALL M3MESG( MESG )

C.........  Compare the consistency of met headers from BTMP, MTMP, PLAY to
C           make sure they are using the same met file

        IF ( METFLAG ) THEN

            IF ( (BMETHEAD(1).EQ.MMETHEAD(1)) .AND.   ! All three (bo, mb, pt)
     &           (BMETHEAD(2).EQ.MMETHEAD(2)) .AND.   ! are consistent
     &           (MMETHEAD(1).EQ.PMETHEAD(1)) .AND.
     &           (MMETHEAD(2).EQ.PMETHEAD(2))
     &         ) THEN

               WRITE( MESG,92000 ) 'Met file consistency is verified '
     &                          // 'among BTMP, MTMP, and PLAY'
               CALL M3MESG( MESG )
               WRITE( MESG,92020 ) BMETHEAD( 1 )
               CALL M3MESG( MESG )
               WRITE( MESG,92020 ) BMETHEAD( 2 )
               CALL M3MESG( MESG )

            ELSE IF ( (BMETHEAD(1).EQ.MMETHEAD(1)) .AND.
     &                (BMETHEAD(2).EQ.MMETHEAD(2)) ) THEN ! bo & mb are consist.
               IF ( PMETHEAD(1).NE.BMETHEAD(1) ) THEN ! Pt's scen. inconsistent
                  WRITE( MESG,92000 )
     &               'Point source has different met scenario ' //
     &               'than Biogenics and Mobile'
                  CALL M3WARN( 'CSGLDAYMRG', 0, 0, MESG )
               END IF
               IF ( PMETHEAD(2).NE.BMETHEAD(2) ) THEN ! Pt's scheme incons.
                  WRITE( MESG,92000 )
     &               'Point source has different cloud scheme ' //
     &               'than Biogenics and Mobile'
                  CALL M3WARN( 'CSGLDAYMRG', 0, 0, MESG )
               END IF

               WRITE( MESG,92000 ) 'For Point Sources:'
               CALL M3MESG( MESG )
               WRITE( MESG,92020 ) PMETHEAD( 1 )
               CALL M3MESG( MESG )
               WRITE( MESG,92020 ) PMETHEAD( 2 )
               CALL M3MESG( MESG )
               WRITE( MESG,92000 ) 'For Biogenics and Mobile:'
               CALL M3MESG( MESG )
               WRITE( MESG,92020 ) BMETHEAD( 1 )
               CALL M3MESG( MESG )
               WRITE( MESG,92020 ) BMETHEAD( 2 )
               CALL M3MESG( MESG )

            ELSE IF ( (MMETHEAD(1).EQ.PMETHEAD(1)) .AND.
     &                (MMETHEAD(2).EQ.PMETHEAD(2)) ) THEN ! Try mb & pt
               IF ( BMETHEAD(1).NE.MMETHEAD(1) ) THEN ! Bo's scen. inconsistent
                  WRITE( MESG,92000 )
     &               'Biogenics has different met scenario ' //
     &               'than Mobile and Point'
                  CALL M3WARN( 'CSGLDAYMRG', 0, 0, MESG )
               END IF
               IF ( BMETHEAD(2).NE.MMETHEAD(2) ) THEN ! Bo's scheme incons.
                  WRITE( MESG,92000 )
     &               'Biogenics has different cloud scheme ' //
     &               'than Mobile and Point'
                  CALL M3WARN( 'CSGLDAYMRG', 0, 0, MESG )
               END IF

               WRITE( MESG,92000 ) 'For Biogenics:'
               CALL M3MESG( MESG )
               WRITE( MESG,92020 ) BMETHEAD( 1 )
               CALL M3MESG( MESG )
               WRITE( MESG,92020 ) BMETHEAD( 2 )
               CALL M3MESG( MESG )
               WRITE( MESG,92000 ) 'For Mobile and Point:'
               CALL M3MESG( MESG )
               WRITE( MESG,92020 ) MMETHEAD( 1 )
               CALL M3MESG( MESG )
               WRITE( MESG,92020 ) MMETHEAD( 2 )
               CALL M3MESG( MESG )

            ELSE IF ( (PMETHEAD(1).EQ.BMETHEAD(1)) .AND.
     &                (PMETHEAD(2).EQ.BMETHEAD(2)) ) THEN ! Try pt & bo
               IF ( MMETHEAD(1).NE.PMETHEAD(1) ) THEN ! Mb's scen. inconsistent
                  WRITE( MESG,92000 )
     &               'Mobile has different met scenario ' //
     &               'than Point and Biogenics'
                  CALL M3WARN( 'CSGLDAYMRG', 0, 0, MESG )
               END IF
               IF ( MMETHEAD(2).NE.PMETHEAD(2) ) THEN ! Mb's scheme incons.
                  WRITE( MESG,92000 )
     &               'Mobile has different cloud scheme ' //
     &               'than Point and Biogenics'
                  CALL M3WARN( 'CSGLDAYMRG', 0, 0, MESG )
               END IF

               WRITE( MESG,92000 ) 'For Mobile:'
               CALL M3MESG( MESG )
               WRITE( MESG,92020 ) MMETHEAD( 1 )
               CALL M3MESG( MESG )
               WRITE( MESG,92020 ) MMETHEAD( 2 )
               CALL M3MESG( MESG )
               WRITE( MESG,92000 ) 'For Point and Biogenics:'
               CALL M3MESG( MESG )
               WRITE( MESG,92020 ) PMETHEAD( 1 )
               CALL M3MESG( MESG )
               WRITE( MESG,92020 ) PMETHEAD( 2 )
               CALL M3MESG( MESG )

            ELSE     ! Neither two inputs use the same met data
               WRITE( MESG,92000 ) 'Biogenic, Mobile and Point all '
     &                          // 'used different met inputs '
               CALL M3WARN( 'CSGLDAYMRG', 0, 0, MESG )

               WRITE( MESG,92000 )
     &            'For Biogenics:'
               CALL M3MESG( MESG )
               WRITE( MESG,92020 ) BMETHEAD( 1 )
               CALL M3MESG( MESG )
               WRITE( MESG,92020 ) BMETHEAD( 2 )
               CALL M3MESG( MESG )

               WRITE( MESG,92000 )
     &            'For Mobile:'
               CALL M3MESG( MESG )
               WRITE( MESG,92020 ) MMETHEAD( 1 )
               CALL M3MESG( MESG )
               WRITE( MESG,92020 ) MMETHEAD( 2 )
               CALL M3MESG( MESG )

               WRITE( MESG,92000 )
     &            'For Point Source:'
               CALL M3MESG( MESG )
               WRITE( MESG,92020 ) PMETHEAD( 1 )
               CALL M3MESG( MESG )
               WRITE( MESG,92020 ) PMETHEAD( 2 )
               CALL M3MESG( MESG )

            END IF   ! Met header comparison

        END IF       ! if METFLAG

        LDATE = -1

        DO  229  T = 1, NSTEPS

            IF ( MOD( T,24 ) .EQ. 1 ) THEN
 
                MESG = 'Processing ' //
     &                 DAYS( WKDAY( JDATE ) ) // MMDDYY( JDATE )
                CALL M3MSG2( MESG( 1:TRIMLEN( MESG ) ) )
 
            END IF
 
C.............  Write to screen because WRITE3 only writes to LDEV
            WRITE( *, 93020 ) HHMMSS( JTIME )
        
            IF ( LDATE .NE. JDATE ) THEN  !  get new day's input files, etc

                LDATE = JDATE

                IF ( AFLAG ) THEN
                    CALL DAYFILE( ANAME, ADATE, ATIME, JDATE, JTIME,
     &                            NASRC, 'Enter  AREA  SOURCE file' )
                END IF

                IF ( PFLAG ) THEN
                    CALL DAYFILE( PNAME, PDATE, PTIME, JDATE, JTIME,
     &                            NPSRC, 'Enter  POINT SOURCE file' )
                END IF

            END IF  !  if ldate different from jdate (so get new input files)

C...........   Read mobile source emissions:

            IF ( MFLAG ) THEN

C.............   Read in entire emissions field
 
               IF( .NOT. READ3( MTNAME, ALLVAR3, ALLAYS3,
     &                    JDATE, JTIME, MOBILV ) ) THEN

                   MESG = 'Could not read timestep from ' // MTNAME
                   CALL M3EXIT( 'CSGLDAYMRG', JDATE, JTIME, MESG, 2 )
                   GO TO 229

               END IF          !  if read3() succeeds or not

            END IF      !  if mflag


C...........   Read point source layer fractions:

            IF ( PFLAG ) THEN

                IF( .NOT. READ3( PLNAME, 'LFRAC', ALLAYS3, 
     &                     JDATE, JTIME, LFRAC ) ) THEN

                    CALL M3EXIT( 'CSGLDAYMRG', JDATE, JTIME, 
     &              'Could not read LFRAC from '// PLNAME, 2 )

                END IF          !  if read3() failed

                DO  112  K = 1, NPSRC  !  transpose layer fraction coeffs:
                DO  112  L = 1, EMLAYS
                    FRACL( L,K ) = LFRAC( K,L )
111             CONTINUE            ! end loop on input layers
112             CONTINUE            ! end loop on input variables

            END IF          !  if pflag


C...........   For each inventory species, read area and mobile emissons
C...........   as appropriate, and for each model species dependent upon
C...........   that species, apply gridding matrices as appropriate and
C...........   write to output file:

            SIDX = 0  ! for MOBILE

            DO  218  V = 1, NIPOL

                IF( AFLAG ) THEN
                    IF( .NOT. READ3( ANAME, EINAM( V ), ALLAYS3, 
     &                         JDATE, JTIME, AREAV ) ) THEN

                        MESG = 'Could not read "' //
     &                  EINAM( V )( 1:TRIMLEN( EINAM( V ))) //
     &                  '" from file ' // ANAME
                        CALL M3EXIT('CSGLDAYMRG', JDATE,JTIME, MESG, 2)
                        GO TO  218

                    END IF          !  if read3() failed
                END IF          !  if aflag:  merging area emissions

                IF ( PFLAG ) THEN

                    IF( .NOT. READ3( PNAME, EINAM( V ), ALLAYS3, 
     &                         JDATE, JTIME, POINTV ) ) THEN

                        MESG = 'Could not read "' //
     &                  EINAM( V )( 1:TRIMLEN( EINAM( V ))) //
     &                  '" from ' // PNAME
                        CALL M3EXIT('CSGLDAYMRG', JDATE,JTIME, MESG, 2)
                        GO TO  218

                    END IF          !  if read3() failed

                END IF          !  if pflag


C...............   Process all output species for this input species:

                DO  207  J = 1,  NSPCOUT( V )

                    DO 197 S = 1, NMSRC        ! Reinit. each src for each 
                       MOBILC( S ) = 0.        ! new model species
197                 CONTINUE

                    K = SPCREF( V,J )
C...................   Initialization of emisg[]:  either read from bio file or 
C...................   set to zero:

                    IF ( BSPCFLAG( K ) ) THEN

                        IF( .NOT. READ3( BTNAME, EMNAM( K ), ALLAYS3,
     &                                   JDATE, JTIME, EMISG ) ) THEN
                            MESG = 'Could not read "' //
     &                      EMNAM( K )( 1:TRIMLEN( EMNAM( K ))) //
     &                      '" from ' // BTNAME
                            CALL M3EXIT( 'CSGLDAYMRG', JDATE,JTIME, 
     &                                   MESG, 2 )
                        END IF          !  if read3() failed
                        DO  192  L = 2, EMLAYS
                        DO  191  I = 1, NGRID
                            EMISG( I,L ) = 0.0
191                     CONTINUE
192                     CONTINUE

                    ELSE

                        DO  194  L = 1, EMLAYS
                        DO  193  I = 1, NGRID
                            EMISG( I,L ) = 0.0
193                     CONTINUE
194                     CONTINUE

                    END IF              !  if this is a biogenic species


                    IF ( AFLAG ) THEN

                        CALL DMATVEC( NASRC, SA( 1,K ), AREAV, AREAC )
                        CALL AMRGVEC( NASRC, AREAC, NMATX, NA, IA, MA,
     &                                NGRID, EMISG )

                    END IF      !  if aflag


                    IF ( MFLAG ) THEN   !  combine mobile modes, vehicle
                                        !  types for this model species:

C........  Loop through processes associated w/ model species (within model
C........       species loop)

               DO 199 JJ = 1, NPPM( K )

C........  Set index of SS for model species and process
                    SIDX = SIDX + 1

C........  Get index of EMISV for model species and process
C........  First find the location of the second "_"
                    L = TRIMLEN( SPCFNM ( SIDX ) )
                    N = INDEX  ( SPCFNM ( SIDX )(   1:L ), '_' )
                    M = INDEX  ( SPCFNM ( SIDX )( N+1:L ), '_' ) + N

C........  Then match the first 2 name sections to EMTNAM
                    NMATCH = SPCFNM( SIDX )( 1:M-1 )
                    KK     = INDEX1( NMATCH, ANETYPE, EMTNAM )

C........  Apply SM to MOBILV and add result to current MOBILC

                    CALL MSMATVEC( NMSRC, NVTYPE,
     &                             SM   ( 1,1,SIDX ),
     &                             MOBILV( 1,1,KK   ), MOBILC )

 199           CONTINUE

C........  Apply gridding matrix

                       CALL AMRGVEC ( NMSRC, MOBILC, NMATX, 
     &                                NM, IM, MM,
     &                                NGRID, EMISG )

                    END IF      !  if mflag


                    IF ( PFLAG ) THEN

                        CALL DMATVEC( NPSRC, SP( 1,K ), POINTV, POINTC )
                        CALL PLMVEC ( FRACL, POINTC, NP, IP, EMISG )

                    END IF      !  if pflag

                    IF ( OUTFLAG ) THEN

                        IF ( .NOT. WRITE3( ENAME, EMNAM( K ),
     &                                     JDATE, JTIME, 
     &                                     EMISG ) ) THEN

                            MESG = 'Could not write "' // 
     &                      EMNAM( K )( 1:TRIMLEN( EMNAM( K ) ) ) //
     &                      '" to ' // ENAME
                            CALL M3WARN( 'CSGLDAYMRG', JDATE, JTIME,
     &                                   MESG )

                        END IF          !  if write3() failed

                    END IF          !  if outputflag

207             CONTINUE    !  end loop on output vbles for this input vble

218         CONTINUE            ! end loop on input variables

            CALL NEXTIME( JDATE, JTIME, TSTEP )         !  update model clock
            CALL NEXTIME( ADATE, ATIME, TSTEP )         !  update  area clock
            CALL NEXTIME( BDATE, BTIME, TSTEP )         !  update  bio  clock
            CALL NEXTIME( MDATE, MTIME, TSTEP )         !  update mobil clock
            CALL NEXTIME( PDATE, PTIME, TSTEP )         !  update point clock

229     CONTINUE          !  end loop on time steps

999     CONTINUE          !  exit program
      
        CALL M3EXIT( 'CSGLDAYMRG', 0, 0, 
     &               'Normal completion  of PROGRAM CSGLDAYMRG', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Error and warning message formats..... 91xxx

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92020   FORMAT ( 10X, A )

C...........   Formatted file I/O formats............ 93xxx

93020   FORMAT( 8X, 'at time ', A8 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I10, :, 2X ) )

94020   FORMAT( A, I4, 2X, 10 ( A, :, 1PG14.6, :, 2X ) )


C...........   Miscellaneous formats................. 95xxx

        END


C****************   end program CSGLDAYMRG   *****************************
C***********************************************************************

        SUBROUTINE  DAYFILE( FNAME, FDATE, FTIME, JDATE, JTIME, 
     &                       NSRCS, PROMPT )

C***********************************************************************
C  subroutine body starts at line 1056
C
C  DESCRIPTION:
C       Prompts the user for a day-specific emissions file FNAME,
C       and for the date&time to be used in FNAME to correspond
C       with model date&time JDATE:JTIME.  Opens the requested file
C       and checks consistency of the file's dimensioning with the 
C       dimensioning expected by the caller.
C       Sets file description into the FDESC commons.
C
C  PRECONDITIONS REQUIRED:  
C       setenv's for appropriate logical names
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; GETYN, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  12/96 by CJC for new CSGLDAYMRG
C
C***********************************************************************
        
C...........   INCLUDES:

      INCLUDE 'PARMS3.EXT'      ! I/O API constants
      INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
      INCLUDE 'IODECL3.EXT'     ! I/O API function declarations


C...........   ARGUMENTS and their descriptions:

        CHARACTER*16    FNAME   !  [default] name of the file opened
        INTEGER         FDATE   !  this day's file-starting date
        INTEGER         FTIME   !  this day's file-starting time
        INTEGER         JDATE   !  current model date
        INTEGER         JTIME   !  current model time
        INTEGER         NSRCS   !  number of sources expected
        CHARACTER*(*)   PROMPT  !  prompt for user


C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         GETNUM
        LOGICAL         GETYN
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN

        EXTERNAL        GETNUM, GETYN, PROMPTMFILE, TRIMLEN


C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        CHARACTER*240   MESG


C***********************************************************************
C   begin body of subroutine  DAYFILE

        FNAME = PROMPTMFILE( PROMPT, FSREAD3, FNAME, 
     &                       'CSGLDAYMRG/DAYFILE' )
        
        IF ( .NOT. DESC3( FNAME ) ) THEN
            MESG =  'Could not get description for file "' //
     &              FNAME( 1:TRIMLEN( FNAME ) ) // '"'
            CALL M3EXIT( 'CSGLDAYMRG/DAYFILE', JDATE, JTIME, MESG, 2 )
        END IF          !  if DESC3() failed

        IF ( NROWS3D .NE. NSRCS ) THEN
            WRITE( MESG, 94010 )
     &          'Dimension mismatch in file "',
     &          FNAME( 1:TRIMLEN( FNAME ) ), '" -- file:', NROWS3D,
     &          'caller:', NSRCS
            CALL M3EXIT( 'CSGLDAYMRG/DAYFILE', JDATE, JTIME, MESG, 2 )
        END IF          !  if dimension mismatch
        
        IF ( TSTEP3D .GT. 0 ) THEN
            WRITE( *,92000 ) ' '
            WRITE( *,92010 )
     &          'The current model   date&time is', JDATE,   JTIME,
     &          'This file starting  date&time is', SDATE3D, STIME3D,
     &          'Enter this-file date&time to be used with ' //
     &          'current model date&time.'
            FDATE = GETNUM( SDATE3D, 9999999, SDATE3D, 
     &                      'Enter date (YYYYDDD)' )
            FTIME = GETNUM( 0, 235959, JTIME, 
     &                      'Enter time  (HHMMSS)' )
        ELSE
            FDATE = 0
            FTIME = 0
        END IF          !  if FNAME time-stepped, or not.

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT ( 5X , A )    !  generic string log format


92010   FORMAT ( 5X , A, :, I9, ':', I6.6 )



C...........   Internal buffering formats............ 94xxx

94010   FORMAT ( 3A, I9, 2X, A, I9 )    !  dimension mismatch format


        END


C****************   end subroutine DAYFILE  *****************************
C***********************************************************************


        SUBROUTINE  CTLSPC( FNAME, NSRCS, SPCMAT, CTLMAT, INVNAME )

C***********************************************************************
C  subroutine body starts at line 1166
C
C  DESCRIPTION:
C       Reads the indicated control matrix from FNAME, and applies
C       it to the indicated speciation matrix.
C
C  PRECONDITIONS REQUIRED:  
C       setenv's for appropriate logical names
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C
C  REVISION  HISTORY:
C       Prototype  12/96 by CJC for new CSGLDAYMRG
C
C***********************************************************************
        
C...........   INCLUDES:

        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   ARGUMENTS and their descriptions:

        CHARACTER*16    FNAME   !  of the control-matrix file
        INTEGER         NSRCS   !  number of sources expected
        REAL            SPCMAT( NSRCS, NMPOL )
        REAL            CTLMAT( NSRCS, NIPOL )


C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         TRIMLEN
        EXTERNAL        TRIMLEN


C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        CHARACTER*16    SCRBUF
        CHARACTER*80    LINE
        CHARACTER*240   MESG
        
        INTEGER         S, I, V

        CHARACTER*16    INVNAME( NMPOL )

C***********************************************************************
C   begin body of subroutine  CTLSPC
        
        IF ( .NOT. DESC3( FNAME ) ) THEN

            MESG = 'Could not get description of file "' //
     &             FNAME( 1:TRIMLEN( FNAME ) ) // '"'
            CALL M3EXIT( 'CSGLDAYMRG', 0, 0, MESG, 2 )

        ELSE IF ( NROWS3D .NE. NSRCS ) THEN

            LINE = 'Dimension mismatch.  "' // 
     &             FNAME( 1:TRIMLEN( FNAME ) ) // '" file:' 
            WRITE( MESG, 94010 )
     &              LINE( 1:TRIMLEN( LINE ) ), NROWS3D,
     &              'program:', NSRCS
            CALL M3EXIT( 'CSGMERGE', 0, 0, MESG, 2 )

        END IF      !  check dimensionality of FNAME

        LINE = 'Reading and applying "' // 
     &         FNAME( 1:TRIMLEN( FNAME ) ) // '"'
        CALL M3MESG( LINE )
        
        DO  33  S = 1, NIPOL
            
            SCRBUF = 'CONTROL_' // EINAM( S )( 1:8 )
            IF ( .NOT. READ3( FNAME, SCRBUF, 1, 0, 0, 
     &                        CTLMAT( 1,S ) ) ) THEN
                MESG = 'Could not read "' //
     &                  SCRBUF( 1:TRIMLEN( SCRBUF ) ) //
     &                  '" from CONTROL FILE "' //
     &                  FNAME( 1:TRIMLEN( FNAME ) ) // '"'
                CALL M3EXIT( 'CSGLDAYMRG/CTLSPC', 0, 0, MESG, 2 )
            END IF              !  if read3() failed

C...........   If this inventory species name matches the leading portion
C...........   of the spc-matrix-vble name, multiply this diagonal
C...........   spc-conversion matrix by this inventory-spc control matrix

            DO  22  V = 1, NMPOL
                IF ( EINAM( S ) .EQ. INVNAME( V ) ) THEN
                    DO  11  I = 1, NSRCS
                        SPCMAT( I,V ) = CTLMAT( I,S ) * SPCMAT( I,V )
11                  CONTINUE
                END IF
22          CONTINUE

33      CONTINUE

        RETURN

C...........   Internal buffering formats............ 94xxx

94010   FORMAT ( A, :, I9, :, 2X )    !  dimension mismatch format

        END     !  end subroutine CTLSPC

