
C Version "@(#)$Id$ $Source$ $Date$ 

        PROGRAM CSGLMERGE

C***********************************************************************
C  program body starts at line 179
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
C       Version    1/97 by SL for met header implementation and header
C       consistency checking scheme.
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'ARDIMS3.EXT'   !  area -source dimensioning parameters
        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'LADIMS3.EXT'   !  vertical dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'MBDIMS3.EXT'   !  mobile source parms and tables
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        LOGICAL         DSCGRID
        LOGICAL         GETYN
        INTEGER         INDEX1
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN
        
        EXTERNAL        DSCGRID, GETYN, INDEX1, PROMPTMFILE, TRIMLEN


C...........   PARAMETER

        CHARACTER*16    NONE
        PARAMETER     ( NONE = 'NONE' )


C...........   LOCAL VARIABLES and their descriptions:
C...........   Emissions Sources input and output arrays (one variable at a time)

        REAL            AREAV ( NASRC )  !  input area emissions.
        REAL            AREAC ( NASRC )  !  controlled speciated area emissions

        REAL            POINTV( NPSRC )  !  input point emissions.
        REAL            POINTC( NPSRC )  !  controlled speciated point emissions

cmh temp
        INTEGER NETYPE
        PARAMETER (NETYPE=6)
cmh temp done

        REAL            MOBILV( NMSRC, NVTYPE, NETYPE )  !  input mobile emissions.
        REAL            MOBILC( NMSRC )  !  controlled speciated mobile emissions

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

C...........   Speciation Matrices, and description of how it maps input 
C...........   variables to output variables:  SS( S, SPCREF( V,W ) )
C...........   is used to produce species # OUTREF in the output file
C...........   from species V in the input file
C...........   (NOTE that names in the speciation-matrix file are of the 
C...........   form "<INV-SPC-NAME>to<MODEL-SPC-NAME>")

        REAL            SA( NASRC, NMPOL )
        REAL            SP( NPSRC, NMPOL )

        REAL            SM( NMSRC, NVTYPE, MXMPOL*MXEPRC )
        CHARACTER*16    MSNAMES( MXMPOL*MXEPRC )
        INTEGER         MBLMODE( MXMPOL*MXEPRC )
        INTEGER         MBLIPOL( MXMPOL*MXEPRC )
        INTEGER         MBLSDEX( NMPOL )
        INTEGER         MBLSCNT( NMPOL )

        INTEGER         NSPCOUT( NIPOL )        !  # of active rows in SPCREF
        INTEGER         SPCREF ( NIPOL, NMPOL )
        
        CHARACTER*16    INVNAME( NMPOL )        !  inv-spc part of name

C...........   Other local variables
        
        INTEGER         JDATE, JTIME, TSTEP
        INTEGER         NSTEPS
        INTEGER         LDEV
        INTEGER         I, J, K, L, M, N, F, S, T, V
        
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
        LOGICAL         BSPCFLAG( NMPOL )	!  is spc in bio input file?
        LOGICAL         OUTFLAG !  generating file output

        CHARACTER*16    ANAME   !  logical name for  area-source      input file
        CHARACTER*16    BNAME   !  logical name for  biogenic-source  input file
        CHARACTER*16    MNAME   !  logical name for  mobile-source    input file
        CHARACTER*16    PNAME   !  logical name for point-source      input file
        CHARACTER*16    ACNAME  !  logical name for control    matrix input file
        CHARACTER*16    AGNAME  !  logical name for grid       matrix input file
        CHARACTER*16    APNAME  !  logical name for projection matrix input file
        CHARACTER*16    ASNAME  !  logical name for speciation matrix input file
        CHARACTER*16    MCNAME  !  logical name for projection matrix input file
        CHARACTER*16    MGNAME  !  logical name for grid       matrix input file
        CHARACTER*16    MPNAME  !  logical name for projection matrix input file
        CHARACTER*16    MSNAME  !  logical name for speciation matrix input file
        CHARACTER*16    PCNAME  !  logical name for control    matrix input file
        CHARACTER*16    PGNAME  !  logical name for grid       matrix input file
        CHARACTER*16    PLNAME  !  logical name for layer coefficient input file
        CHARACTER*16    PPNAME  !  logical name for projection matrix input file
        CHARACTER*16    PSNAME  !  logical name for speciation matrix input file
        CHARACTER*16    ENAME   !  logical name for emission output file
            
        CHARACTER*16    SCRBUF, TMPBUF, ABUF, BBUF  !  scratch vbles for names
        CHARACTER*80    BMETHEAD( 2 )  ! met headers from BTMP
        CHARACTER*80    MMETHEAD( 2 )  ! met headers from MTMP
        CHARACTER*80    PMETHEAD( 2 )  ! met headers from PLAY
        CHARACTER*256   MESG


C***********************************************************************
C   begin body of program CSGLMERGE
        
        LDEV = INIT3()
        
        WRITE( *,92000 ) ' ', ' ', ' ',
     &'Program CSGLMERGE to take combinations of the sorted source    ',
     &'level area, bio, mobile and point source emissions files and  ',
     &'gridding-transform, speciation-transform, point-source layer  ',
     &'fractions if appropriate, and optionally control and projection',
     &'matrices, and produce layered gridded model emissions.  The    ',
     &'output time step structure will be the same asthe time step   ',
     &'structure of the input emissions files (which must also match ',
     &'with each other)',
     &  ' '
        WRITE( *,92000 ) 
     &'You will need to have set the logical names for the input and ',
     &'output files prior to program launch, using                   ',
     &'   "setenv <logicalname> <pathname>".                         ',
     &'You will be prompted for these logical names during program   ',
     &'execution.  Use "NONE" as the name for the control-matrix     ',
     &'file if you wish to omit control from the operations          ',
     &'performed, or as the name of the output file if you want only ',
     &'to time program performance without the overhead of output I/O',
     &' '
        WRITE( *,92000 ) 
     &'You may use END_OF-FILE (control-D) to quit the program       ',
     &'during logical-name entry.  Default responses are indicated   ',
     &'in brackets [LIKE THIS].                                      ',
     &' '
        
        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            WRITE( *,92000 ) ' ', 'Ending program CSGLMERGE', ' '
            CALL EXIT( 0 )
        END IF


C.......   Open area-related files; set appropriate flags:

        JDATE   = -1
        JTIME   = -1
        TSTEP   = -1
        NSTEPS  = -1

        ANAME = PROMPTMFILE( 
     &          'Enter logical name for  AREA SOURCE file, or "NONE"',
     &          FSREAD3, 'ATMP', 'CSGLMERGE' )

        IF ( ANAME .NE. NONE ) THEN

            AGNAME = PROMPTMFILE( 
     &          'Enter logical name for  AREA GRID MATRIX file',
     &          FSREAD3, 'AGMAT', 'CSGLMERGE' )

            ASNAME = PROMPTMFILE( 
     &          'Enter logical name for  AREA SPECIATION MATRIX file',
     &          FSREAD3, 'ASMAT', 'CSGLMERGE' )

            ACNAME = PROMPTMFILE( 
     &          'Enter logical name for  AREA CONTROL MATRIX file ' //
     &          'or "NONE"',
     &          FSREAD3, 'ACMAT', 'CSGLMERGE' )

            APNAME = PROMPTMFILE( 
     &          'Enter logical name for  AREA PROJECTION MATRIX file ' 
     &          // 'or "NONE"',
     &          FSREAD3, 'APMAT', 'CSGLMERGE' )

            AFLAG  = .TRUE.
            ACFLAG = ( ACNAME .NE. NONE )
            APFLAG = ( APNAME .NE. NONE )

            IF ( .NOT. DESC3( ANAME ) ) THEN
                MESG = 'Could not get description of file ' // ANAME
                CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )
            END IF

            JDATE = SDATE3D
            JTIME = STIME3D
            TSTEP = TSTEP3D
            IF ( TSTEP .EQ. 0 ) THEN
                NSTEPS = 1
            ELSE
                NSTEPS = MXREC3D
            END IF

        ELSE

            AFLAG  = .FALSE.
            ACFLAG = .FALSE.
            APFLAG = .FALSE.

        END IF          !  if aname is "none " or not.


C.......   Open biogenic file; set appropriate flags:

        BNAME = PROMPTMFILE( 
     &          'Enter logical name for BIOGENIC SOURCE file, ' //
     &          'or "NONE"',
     &          FSREAD3, 'BTMP', 'CSGLMERGE' )

        IF ( BNAME .NE. NONE ) THEN

            BFLAG = .TRUE.

            IF ( .NOT. DESC3( BNAME ) ) THEN
                MESG = 'Could not get description of file ' // BNAME
                CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )
            END IF

            BMETHEAD( 1 ) = FDESC3D ( 2 )  ! record met header from BTMP
            BMETHEAD( 2 ) = FDESC3D ( 3 )

            DO  11  V = 1, NMPOL
                BSPCFLAG( V ) = ( INDEX1( EMNAM( V ), 
     &                                    NVARS3D, VNAME3D ) .GT. 0 )
11          CONTINUE

            IF ( JDATE .GE. 0 ) THEN    !  dates and times defined
                IF ( ( JDATE .LT. SDATE3D ) .OR.
     &               ( JDATE .EQ. SDATE3D .AND. 
     &                 JTIME .LT. STIME3D ) ) THEN
                    JDATE   = SDATE3D
                    JTIME   = STIME3D
                    NSTEPS  = MIN( NSTEPS, MXREC3D )
                END IF
                IF ( TSTEP .NE. TSTEP3D ) THEN
                    WRITE( MESG,94010 ) 
     &              'Time step mismatch:', TSTEP, 
     &              'from area,', TSTEP3D,
     &              'from bio'
                    CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )
                END IF
            ELSE                        !  dates and times not yet defined
                JDATE   = SDATE3D
                JTIME   = STIME3D
                TSTEP   = TSTEP3D
                IF ( TSTEP .EQ. 0 ) THEN
                    NSTEPS  = 1
                ELSE
                    NSTEPS  = MXREC3D
                END IF
            END IF                      !  if dates and times yet defined or not

        ELSE

            BFLAG = .FALSE.
            DO  12  V = 1, NMPOL
                BSPCFLAG( V ) = .FALSE.
12          CONTINUE

        END IF


C.......   Open mobile-related files; set appropriate flags:

        MNAME = PROMPTMFILE( 
     &          'Enter logical name for MOBILE SOURCE file, or "NONE"',
     &          FSREAD3, 'MTMP', 'CSGLMERGE' )

        IF ( MNAME .NE. NONE ) THEN

            MGNAME = PROMPTMFILE( 
     &          'Enter logical name for MOBILE GRID MATRIX file',
     &          FSREAD3, 'MGMAT', 'CSGLMERGE' )

            MSNAME = PROMPTMFILE( 
     &          'Enter logical name for MOBILE SPECIATION MATRIX file',
     &          FSREAD3, 'MSMAT', 'CSGLMERGE' )

            MCNAME = PROMPTMFILE( 
     &          'Enter logical name for MOBILE CONTROL MATRIX file ' //
     &          'or "NONE"',
     &          FSREAD3, 'MCMAT', 'CSGLMERGE' )

            MPNAME = PROMPTMFILE( 
     &          'Enter logical name for MOBILE PROJECTION MATRIX file ' 
     &          // 'or "NONE"',
     &          FSREAD3, 'MPMAT', 'CSGLMERGE' )

            IF ( .NOT. DESC3( MNAME ) ) THEN
                MESG = 'Could not get description of file ' // MNAME
                CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )
            END IF

            MMETHEAD( 1 ) = FDESC3D ( 2 )  ! record met header from MTMP
            MMETHEAD( 2 ) = FDESC3D ( 3 )

            IF ( JDATE .GE. 0 ) THEN    !  dates and times defined
                IF ( ( JDATE .LT. SDATE3D ) .OR.
     &               ( JDATE .EQ. SDATE3D .AND. 
     &                 JTIME .LT. STIME3D ) ) THEN
                    JDATE   = SDATE3D
                    JTIME   = STIME3D
                    NSTEPS  = MIN( NSTEPS, MXREC3D )
                END IF
                IF ( TSTEP .NE. TSTEP3D ) THEN
                    WRITE( MESG,94010 ) 
     &              'Time step mismatch:', TSTEP, 
     &              'from area/bio,', TSTEP3D,
     &              'from mobile'
                    CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )
                END IF
            ELSE                        !  dates and times not yet defined
                JDATE   = SDATE3D
                JTIME   = STIME3D
                TSTEP   = TSTEP3D
                IF ( TSTEP .EQ. 0 ) THEN
                    NSTEPS  = 1
                ELSE
                    NSTEPS  = MXREC3D
                END IF
            END IF                      !  if dates and times yet defined or not

            MFLAG  = .TRUE.
            MCFLAG = ( MCNAME .NE. NONE )
            MPFLAG = ( MPNAME .NE. NONE )

        ELSE

            MFLAG  = .FALSE.
            MCFLAG = .FALSE.
            MPFLAG = .FALSE.

        END IF          !  if mname is "none " or not.


C.......   Open point-related files; set appropriate flags:

        PNAME = PROMPTMFILE( 
     &          'Enter logical name for required POINT SOURCE file"',
     &          FSREAD3, 'PTMP', 'CSGLMERGE' )

        IF ( PNAME .NE. NONE ) THEN

            PGNAME = PROMPTMFILE( 
     &          'Enter logical name for POINT GRID MATRIX file',
     &          FSREAD3, 'PGMAT', 'CSGLMERGE' )
     
            PSNAME = PROMPTMFILE( 
     &          'Enter logical name for POINT SPECIATION MATRIX file',
     &          FSREAD3, 'PSMAT', 'CSGLMERGE' )
     
            PLNAME = PROMPTMFILE( 
     &              'Enter logical name for POINT LAYER COEFF file',
     &              FSREAD3, 'PLAY', 'CSGLMERGE' )

            PCNAME = PROMPTMFILE( 
     &          'Enter logical name for POINT CONTROL MATRIX file ' //
     &          'or "NONE"',
     &          FSREAD3, 'PCMAT', 'CSGLMERGE' )

            PPNAME = PROMPTMFILE( 
     &      'Enter logical name for POINT PROJECTION MATRIX file ' //
     &      'or "NONE"',
     &      FSREAD3, 'PCMAT', 'CSGLMERGE' )

            IF ( .NOT. DESC3( PNAME ) ) THEN
                MESG = 'Could not get description of file ' // PNAME
                CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )
            END IF

            IF ( JDATE .GE. 0 ) THEN    !  dates and times defined
                IF ( ( JDATE .LT. SDATE3D ) .OR.
     &               ( JDATE .EQ. SDATE3D .AND. 
     &                 JTIME .LT. STIME3D ) ) THEN
                    JDATE   = SDATE3D
                    JTIME   = STIME3D
                    NSTEPS  = MIN( NSTEPS, MXREC3D )
                END IF
                IF ( TSTEP .NE. TSTEP3D ) THEN
                    WRITE( MESG,94010 ) 
     &              'Time step mismatch:', TSTEP, 
     &              'from area/bio/mobile,', TSTEP3D,
     &              'from point'
                    CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )
                END IF
            ELSE                        !  dates and times not yet defined
                JDATE   = SDATE3D
                JTIME   = STIME3D
                TSTEP   = TSTEP3D
                IF ( TSTEP .EQ. 0 ) THEN
                    NSTEPS  = 1
                ELSE
                    NSTEPS  = MXREC3D
                END IF
            END IF                      !  if dates and times yet defined or not

            PFLAG  = .TRUE.
            PCFLAG = ( PCNAME .NE. NONE )
            PPFLAG = ( PPNAME .NE. NONE )

        ELSE            !  pname is "none"

            PFLAG  = .FALSE.
            PCFLAG = .FALSE.
            PPFLAG = .FALSE.

        END IF          !  if pname is "none " or not.

        IF ( .NOT.( AFLAG .OR. BFLAG .OR. MFLAG .OR. PFLAG ) ) THEN
            CALL M3EXIT( 'CSGLMERGE', 0, 0, 'No input files', 2 )
        END IF

C.........  NOTE: have horizontal grid description from DESC3 calls above

C.........  Get vertical description fro PLNAME file
C.........  Create 3D file if have point sources
        IF( PFLAG ) THEN

            IF ( .NOT. DESC3( PLNAME ) ) THEN
                CALL M3EXIT( 'CSGLMERGE', 0, 0,
     &                  'Could not get description of file "'
     &                  // PLNAME( 1:TRIMLEN( PLNAME ) ) // '"', 2 )
            END IF

            PMETHEAD( 1 ) = FDESC3D ( 2 )  ! record met header from PLAY
            PMETHEAD( 2 ) = FDESC3D ( 3 )
            
C.........  If no point sources, someone has changed the code, so give error
        ELSE

            CALL M3EXIT( 'CSGLMERGE', 0, 0, 
     &         'ERROR: point sources required for vertical parameters',
     &         2 )

        ENDIF

C.........  Compare the consistency of met headers from BTMP, MTMP, PLAY to
C           make sure they are using the same met file

            IF ( (BMETHEAD(1).EQ.MMETHEAD(1)) .AND.   ! All three (bo, mb, pt)
     &           (BMETHEAD(2).EQ.MMETHEAD(2)) .AND.   ! are consistent
     &           (MMETHEAD(1).EQ.PMETHEAD(1)) .AND.
     &           (MMETHEAD(2).EQ.PMETHEAD(2)) 
     &         ) THEN

               WRITE( MESG,92000 ) 'Met file consistency is verified '
     &                          // 'among BTMP, MTMP, and PLAY'
               CALL M3MESG( MESG )
               WRITE( MESG,92010 ) BMETHEAD( 1 )
               CALL M3MESG( MESG )
               WRITE( MESG,92010 ) BMETHEAD( 2 )
               CALL M3MESG( MESG )

            ELSE IF ( (BMETHEAD(1).EQ.MMETHEAD(1)) .AND.
     &                (BMETHEAD(2).EQ.MMETHEAD(2)) ) THEN ! bo & mb are consist.
               IF ( PMETHEAD(1).NE.BMETHEAD(1) ) THEN ! Pt's scen. inconsistent
                  WRITE( MESG,92000 ) 
     &               'Point source has different met scenario ' //
     &               'than Biogenics and Mobile' 
                  CALL M3WARN( 'CSGLMERGE', 0, 0, MESG )
               END IF
               IF ( PMETHEAD(2).NE.BMETHEAD(2) ) THEN ! Pt's scheme incons.
                  WRITE( MESG,92000 ) 
     &               'Point source has different cloud scheme ' //
     &               'than Biogenics and Mobile' 
                  CALL M3WARN( 'CSGLMERGE', 0, 0, MESG )
               END IF 
 
               WRITE( MESG,92000 ) 'For Point Sources:'
               CALL M3MESG( MESG )
               WRITE( MESG,92010 ) PMETHEAD( 1 )
               CALL M3MESG( MESG )
               WRITE( MESG,92010 ) PMETHEAD( 2 )
               CALL M3MESG( MESG )
               WRITE( MESG,92000 ) 'For Biogenics and Mobile:'
               CALL M3MESG( MESG )
               WRITE( MESG,92010 ) BMETHEAD( 1 )
               CALL M3MESG( MESG )
               WRITE( MESG,92010 ) BMETHEAD( 2 )
               CALL M3MESG( MESG )
 
            ELSE IF ( (MMETHEAD(1).EQ.PMETHEAD(1)) .AND.
     &                (MMETHEAD(2).EQ.PMETHEAD(2)) ) THEN ! Try mb & pt
               IF ( BMETHEAD(1).NE.MMETHEAD(1) ) THEN ! Bo's scen. inconsistent
                  WRITE( MESG,92000 ) 
     &               'Biogenics has different met scenario ' //
     &               'than Mobile and Point' 
                  CALL M3WARN( 'CSGLMERGE', 0, 0, MESG )
               END IF
               IF ( BMETHEAD(2).NE.MMETHEAD(2) ) THEN ! Bo's scheme incons.
                  WRITE( MESG,92000 ) 
     &               'Biogenics has different cloud scheme ' //
     &               'than Mobile and Point' 
                  CALL M3WARN( 'CSGLMERGE', 0, 0, MESG )
               END IF 
 
               WRITE( MESG,92000 ) 'For Biogenics:'
               CALL M3MESG( MESG )
               WRITE( MESG,92010 ) BMETHEAD( 1 )
               CALL M3MESG( MESG )
               WRITE( MESG,92010 ) BMETHEAD( 2 )
               CALL M3MESG( MESG )
               WRITE( MESG,92000 ) 'For Mobile and Point:'
               CALL M3MESG( MESG )
               WRITE( MESG,92010 ) MMETHEAD( 1 )
               CALL M3MESG( MESG )
               WRITE( MESG,92010 ) MMETHEAD( 2 )
               CALL M3MESG( MESG )
 
            ELSE IF ( (PMETHEAD(1).EQ.BMETHEAD(1)) .AND.
     &                (PMETHEAD(2).EQ.BMETHEAD(2)) ) THEN ! Try pt & bo
               IF ( MMETHEAD(1).NE.PMETHEAD(1) ) THEN ! Mb's scen. inconsistent
                  WRITE( MESG,92000 ) 
     &               'Mobile has different met scenario ' //
     &               'than Point and Biogenics' 
                  CALL M3WARN( 'CSGLMERGE', 0, 0, MESG )
               END IF
               IF ( MMETHEAD(2).NE.PMETHEAD(2) ) THEN ! Mb's scheme incons.
                  WRITE( MESG,92000 ) 
     &               'Mobile has different cloud scheme ' //
     &               'than Point and Biogenics' 
                  CALL M3WARN( 'CSGLMERGE', 0, 0, MESG )
               END IF 
 
               WRITE( MESG,92000 ) 'For Mobile:'
               CALL M3MESG( MESG )
               WRITE( MESG,92010 ) MMETHEAD( 1 )
               CALL M3MESG( MESG )
               WRITE( MESG,92010 ) MMETHEAD( 2 )
               CALL M3MESG( MESG )
               WRITE( MESG,92000 ) 'For Point and Biogenics:'
               CALL M3MESG( MESG )
               WRITE( MESG,92010 ) PMETHEAD( 1 )
               CALL M3MESG( MESG )
               WRITE( MESG,92010 ) PMETHEAD( 2 )
               CALL M3MESG( MESG )
 
            ELSE     ! Neither two inputs use the same met data 
               WRITE( MESG,92000 ) 'Biogenic, Mobile and Point all '
     &                          // 'used different met inputs '
               CALL M3WARN( 'CSGLMERGE', 0, 0, MESG )

               WRITE( MESG,92000 ) 
     &            'For Biogenics:'
               CALL M3MESG( MESG )
               WRITE( MESG,92010 ) BMETHEAD( 1 ) 
               CALL M3MESG( MESG )
               WRITE( MESG,92010 ) BMETHEAD( 2 ) 
               CALL M3MESG( MESG )

               WRITE( MESG,92000 ) 
     &            'For Mobile:'
               CALL M3MESG( MESG )
               WRITE( MESG,92010 ) MMETHEAD( 1 ) 
               CALL M3MESG( MESG )
               WRITE( MESG,92010 ) MMETHEAD( 2 ) 
               CALL M3MESG( MESG )

               WRITE( MESG,92000 ) 
     &            'For Point Source:'
               CALL M3MESG( MESG )
               WRITE( MESG,92010 ) PMETHEAD( 1 ) 
               CALL M3MESG( MESG )
               WRITE( MESG,92010 ) PMETHEAD( 2 ) 
               CALL M3MESG( MESG )

            END IF   ! Met header comparison             

        IF ( .NOT. DSCGRID( GRDNM, SCRBUF,
     &                      GDTYP3D, P_ALP3D, P_BET3D, P_GAM3D, 
     &                      XCENT3D, YCENT3D, XORIG3D, YORIG3D, 
     &                      XCELL3D, YCELL3D, 
     &                      NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            CALL M3EXIT( 'CSGLMERGE', 0, 0,
     &              'Could not get description of grid ' // GRDNM, 2 )

        END IF

        NVARS3D = NMPOL
        FTYPE3D = GRDDED3   !  shares most of file-description with input file.
        TSTEP3D = TSTEP     !  vertical description from PLNAME
        SDATE3D = JDATE
        STIME3D = JTIME
        GDNAM3D = GRDNM
        NLAYS3D = EMLAYS    !  Should be EMLAYS already, but just in case

        DO  22 V = 1, NMPOL
            VNAME3D( V ) = EMNAM( V )
            UNITS3D( V ) = 'moles/hour'
            VDESC3D( V ) = 'gridded speciated emissions'
            VTYPE3D( V ) = m3REAL
22      CONTINUE
        FDESC3D( 1 ) = 'area source speciated gridded emissions.'
                 
        IF ( ACFLAG .OR. MCFLAG .OR. PCFLAG ) THEN
            ENAME = 'LGTSC'
        ELSE
            ENAME = 'LGTS'
        END IF

        ENAME = PROMPTMFILE( 
     &  'Enter logical name for OUTPUT EMIS file, or "NONE"',
     &          FSUNKN3, ENAME, 'CSGLMERGE' )

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
                CALL M3EXIT( 'CSGLMERGE', 0, 0, 
     &                       'Bad input species "' // 
     &                       INVNAME( V )( 1:TRIMLEN( INVNAME( V ))) //
     &                       '" in ' // ASNAME, 2 )
            END IF

            SCRBUF = SFNAM( V )( K+4:16 )
            J = INDEX1( SCRBUF, NMPOL, EMNAM )
            IF ( J .LT. 0 ) THEN
                CALL M3EXIT( 'CSGLMERGE', 0, 0, 
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
                CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )

            END IF      !  if read3() failed for gridding matrix

            IF ( .NOT. READ3( ASNAME, 'ALL', 1, 0, 0, SA ) ) THEN

                MESG = 'Could not read SPECIATION matrix from file "'//
     &                 ASNAME( 1 : TRIMLEN( ASNAME ) ) // '".'
                CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )

            END IF      !  if read3() failed for speciation matrix

        END IF      !  if aflag

        IF ( MFLAG ) THEN

            WRITE( *,92000 ) 
     &      ' ', 'Reading mobile transform matrices.', ' '

            IF ( .NOT. READ3( MGNAME, 'ALL', 1, 0, 0, NM ) ) THEN

                MESG = 'Could not read gridding matrix from file "' //
     &                 MGNAME( 1 : TRIMLEN( MGNAME ) ) // '".'
                CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )

            END IF      !  if read3() failed for gridding matrix

            IF ( .NOT. READ3( MSNAME, 'ALL', 1, 0, 0, SM ) ) THEN

                MESG = 'Could not read speciation matrix from file "'//
     &                 MSNAME( 1 : TRIMLEN( MSNAME ) ) // '".'
                CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )

            END IF      !  if read3() failed for speciation matrix

C...........   Parse names of the variables in msmat:

            IF ( .NOT. DESC3( MSNAME ) ) THEN
                MESG = 'Could not get description of file ' // MSNAME
                CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )
            END IF

            DO  34  V = 1, NVARS3D

                MSNAMES( V ) = VNAME3D( V )
                SCRBUF       = VNAME3D( V )
                MBLMODE( V ) = INDEX1( SCRBUF, NEMODE, MEMODE )
                IF ( MBLMODE( V ) .LE. 0 ) THEN
                    MESG = 'Bad Mobile mode in vble name "' //
     &                     SCRBUF( 1:TRIMLEN( SCRBUF ) ) // 
     &                     '" from file ' // MSNAME
                    CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )
                END IF

                TMPBUF = SCRBUF( 8:16 )
                I      = INDEX1( TMPBUF, NMPOL, EMNAM )
                IF ( I .LE. 0 ) THEN
                    MESG = 'Bad model spc in spc matrix vble name "' //
     &                     SCRBUF( 1:TRIMLEN( SCRBUF ) ) // 
     &                     '" from file '// MSNAME
                    CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )
                END IF

                MBLSDEX( I ) = MIN( V, MBLSDEX( I ) )
                MBLSCNT( I ) = MBLSCNT( I ) + 1

                K    = INDEX( SFNAM( J ), '_to_' )
                BBUF = SFNAM( I )( 1:K-1 )
                L    = INDEX1( BBUF, NIPOL, EINAM )
                IF ( L .LE. 0 ) THEN
                    MESG = 'No inventory species "' //
     &                     ABUF( 1:TRIMLEN( ABUF ) ) //
     &                     '" for spc matrix vble name "' //
     &                     SCRBUF( 1:TRIMLEN( SCRBUF ) ) // 
     &                     '" from file '// MSNAME
                    CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )
                END IF
                MBLIPOL( V ) = L

34          CONTINUE            !  end loop on v, parsing msmat vble-names

        END IF      !  if mflag


        IF ( PFLAG ) THEN

            WRITE( *,92000 ) 
     &      ' ', 'Reading point transform matrices.', ' '

            IF ( .NOT. READ3( PGNAME, 'ALL', 1, 0, 0, NP ) ) THEN

                MESG = 'Could not read gridding matrix from file "' //
     &                 PGNAME( 1 : TRIMLEN( PGNAME ) ) // '".'
                CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )

            END IF      !  if read3() failed for gridding matrix

            IF ( .NOT. READ3( PSNAME, 'ALL', 1, 0, 0, SP ) ) THEN

                MESG = 'Could not read SPECIATION matrix from file "'//
     &                 PSNAME( 1 : TRIMLEN( PSNAME ) ) // '".'
                CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )

            END IF      !  if read3() failed for speciation matrix

        END IF      !  if pflag


C...........   Apply control and/or projection matrices.
C.......   Area source controls and projections:

        IF( ACFLAG ) THEN

            WRITE( *,92000 ) 
     &      ' ', 'Reading  AREA CONTROL matrix.', ' '

            DO  43  S = 1, NIPOL

                SCRBUF = 'CONTROL_' // EINAM( S )( 1:8 )
                IF ( .NOT. READ3( ACNAME, SCRBUF, 1, 0, 0, 
     &                            CA( 1,S ) ) ) THEN

                    MESG = 'Could not read "' // SCRBUF( 1:12 ) //
     &                     '" from CONTROL FILE "' //
     &                     ACNAME( 1 : TRIMLEN( ACNAME ) ) // '".'
                    CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )

                END IF      !  if read3() failed for control matrix

C...............   If this inventory species name matches the
C...............   leading portion of the speciation-matrix-vble name,
C...............   multiply this diagonal species-conversion matrix
C...............   by this inventory-species control matrix
C...............   (recall that INVNAME(V) = VNAME3D(V)(1:5) --
C...............   e.g., the "42101" in "42101toCO"

                DO  42  V = 1, NMPOL
                    IF ( EINAM( S ) .EQ. INVNAME( V ) ) THEN
                        DO  41  I = 1, NASRC
                            SA( I,V ) = CA( I,S ) * SA( I,V )
41                      CONTINUE
                    END IF
42              CONTINUE

43          CONTINUE    ! end loop on species S, for controlarea

        END IF          !  if acflag:  control area


        IF( APFLAG ) THEN

            WRITE( *,92000 ) 
     &      ' ', 'Reading  AREA PROJECTION matrix.', ' '

            DO  46  S = 1, NIPOL

                SCRBUF = 'PROJECT_' // EINAM( S )( 1:8 )
                IF ( .NOT. READ3( APNAME, SCRBUF, 1, 0, 0, 
     &                            CA( 1,S ) ) ) THEN
                       
                    MESG = 'Could not read "' // SCRBUF( 1:12 ) //
     &                     '" from PROJECTION FILE "' //
     &                     APNAME( 1 : TRIMLEN( APNAME ) ) // '".'
                    CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )

                END IF      !  if read3() failed for control matrix

C...............   If this inventory species name matches the
C...............   leading portion of the speciation-matrix-vble name,
C...............   multiply this diagonal species-conversion matrix
C...............   by this inventory-species control matrix
C...............   (recall that INVNAME(V) = VNAME3D(V)(1:5) --
C...............   e.g., the "42101" in "42101toCO"

                DO  45  V = 1, NMPOL
                    IF ( EINAM( S ) .EQ. INVNAME( V ) ) THEN
                        DO  44  I = 1, NASRC
                            SA( I,V ) = CA( I,S ) * SA( I,V )
44                      CONTINUE
                    END IF
45              CONTINUE

46          CONTINUE    ! end loop on species S, for controlarea

        END IF          !  if acflag:  control area

C.......   Mobile source controls and projections:

        IF( MCFLAG ) THEN

            IF ( .NOT. DESC3( MCNAME ) ) THEN
                MESG = 'Could not get description of file ' // MCNAME
                CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )
            END IF

            DO  50  I = 1, NIPOL + NEMODE - 1
                MCNAMES( I ) = VNAME3D( I )( 5:16 )
50          CONTINUE

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
                    CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )       
                END IF

                DO  52  L = 1, NVTYPE
                DO  51  S = 1, NMSRC
                    SM( S,L,I ) = SM( S,L,I ) * CM( S,L,M )  
51              CONTINUE
52              CONTINUE

53          CONTINUE

        END IF          !  if mcflag:  control mobile


        IF( MPFLAG ) THEN

            IF ( .NOT. DESC3( MPNAME ) ) THEN
                MESG = 'Could not get description of file ' // MPNAME
                CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )
            END IF

            DO  54  I = 1, NIPOL + NEMODE - 1
                MCNAMES( I ) = VNAME3D( I )( 5:16 )
54          CONTINUE

C...............   Parse the variable-names from the speciation matrix and 
C...............   map to inventory species and mobile emissions modes:

            DO  57  I = 1, MXMPOL*MXEPRC

                L = MBLIPOL( I )
                M = MBLMODE( I )

                TMPBUF = 'PRJ_' // MEMODE( M )( 1:3 ) // 
     &                   '_' // EINAM( L )
                M = INDEX1( TMPBUF, NIPOL + NEMODE - 1, MCNAMES )
                IF ( M .LE. 0 ) THEN                             
                    MESG = 'Missing projection variable "' //
     &                     TMPBUF( 1:TRIMLEN( TMPBUF ) ) //
     &                     '" FOR spc matrix vble name "' //
     &                     VNAME3D( I )( 1:TRIMLEN( VNAME3D( I ) ) ) //
     &                     '"'
                    CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )       
                END IF

                DO  56  L = 1, NVTYPE
                DO  55  S = 1, NMSRC
                    SM( S,L,I ) = SM( S,L,I ) * CM( S,L,M )  
55              CONTINUE
56              CONTINUE

57          CONTINUE

        END IF          !  if mpflag:  project mobile



C.......   Point source controls and projections:

        IF( PCFLAG ) THEN

            WRITE( *,92000 ) 
     &      ' ', 'Reading POINT CONTROL matrix.', ' '

            DO  63  S = 1, NIPOL

                SCRBUF = 'CONTROL_' // EINAM( S )( 1:8 )
                IF ( .NOT. READ3( PCNAME, SCRBUF, 1, 0, 0, 
     &                            CP( 1,S ) ) ) THEN

                    MESG = 'Could not read "' // SCRBUF( 1:12 ) //
     &                     '" from CONTROL FILE "' //
     &                     PCNAME( 1 : TRIMLEN( PCNAME ) ) // '".'
                    CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )

                END IF      !  if read3() failed for control matrix

C...............   If this inventory species name matches the
C...............   leading portion of the speciation-matrix-vble name,
C...............   multiply this diagonal species-conversion matrix
C...............   by this inventory-species control matrix
C...............   (recall that INVNAME(V) = VNAME3D(V)(1:5) --
C...............   e.g., the "42101" in "42101toCO"

                DO  62  V = 1, NMPOL
                    IF ( EINAM( S ) .EQ. INVNAME( V ) ) THEN
                        DO  61  I = 1, NPSRC
                            SP( I,V ) = CP( I,S ) * SP( I,V )
61                      CONTINUE
                    END IF
62              CONTINUE

63          CONTINUE    !  end loop on species S, for controlpt

        END IF          !  if pcflag:  control point sources

        IF( PPFLAG ) THEN

            WRITE( *,92000 ) 
     &      ' ', 'Reading POINT PROJECTION matrix.', ' '

            DO  66 S = 1, NIPOL

                SCRBUF = 'CONTROL_' // EINAM( S )( 1:8 )
                IF ( .NOT. READ3( PPNAME, SCRBUF, 1, 0, 0, CP( 1,S ) ) 
     &             ) THEN

                    MESG = 'Could not read "' // SCRBUF( 1:12 ) //
     &                     '" from PROJECTION FILE "' //
     &                     PPNAME( 1 : TRIMLEN( PPNAME ) ) // '".'
                    CALL M3EXIT( 'CSGLMERGE', 0, 0, MESG, 2 )

                END IF      !  if read3() failed for control matrix

C...............   If this inventory species name matches the
C...............   leading portion of the speciation-matrix-vble name,
C...............   multiply this diagonal species-conversion matrix
C...............   by this inventory-species control matrix
C...............   (recall that INVNAME(V) = VNAME3D(V)(1:5) --
C...............   e.g., the "42101" in "42101toCO"

                DO  65 V = 1, NMPOL
                    IF ( EINAM( S ) .EQ. INVNAME( V ) ) THEN
                        DO  64  I = 1, NPSRC
                            SP( I,V ) = CP( I,S ) * SP( I,V )
64                      CONTINUE
                    END IF
65              CONTINUE

66          CONTINUE    !  end loop on species S, for controlpt

        END IF          !  if ppflag:  project point sources


C.......   Transform and write out area source emissions values:

        IF ( OUTFLAG ) THEN
            WRITE( *,92000 ) 
     &      ' ', 
     &      'Calculating and writing speciated gridded emissions...',
     &      ' '
        ELSE 
            WRITE( *,92000 ) 
     &      ' ', 
     &      'Calculating speciated gridded emissions...',
     &      ' '
        END IF

        DO  199  T = 1, NSTEPS

C...........   Read and transpose point source layering fraction coeffs:
        
            IF( .NOT. READ3( PLNAME, 'LFRAC', ALLAYS3, 
     &                 JDATE, JTIME, LFRAC ) ) THEN
           
                CALL M3EXIT( 'CSGLMERGE', JDATE, JTIME, 
     &                       'Could not read LFRAC from '// PLNAME, 2 )
           
            END IF          !  if read3() failed

            DO  112  K = 1, NPSRC
            DO  112  L = 1, EMLAYS
                FRACL( L,K ) = LFRAC( K,L )
111         CONTINUE            ! end loop on input variables
112         CONTINUE            ! end loop on input variables

C...........   Read mobile source emissions:

            IF ( MFLAG ) THEN

                IF( .NOT. READ3( MNAME, 'ALL', ALLAYS3, 
     &                     JDATE, JTIME, MOBILV ) ) THEN

                    MESG = 'Could not read timestep from ' // MNAME
                    CALL M3EXIT( 'CSGLMERGE', JDATE, JTIME, MESG, 2 )
                    GO TO  199

                END IF          !  if read3() failed
            END IF      !  if mflag

C...........   For each inventory species, read area and mobile emissons
C...........   as appropriate, and for each model species dependent upon
C...........   that species, apply gridding matrices as appropriate and
C...........   write to output file:

            DO  198  V = 1, NIPOL

                IF( AFLAG ) THEN
                    IF( .NOT. READ3( ANAME, EINAM( V ), ALLAYS3, 
     &                         JDATE, JTIME, AREAV ) ) THEN

                        MESG = 'Could not read "' //
     &                  EINAM( V )( 1:TRIMLEN( EINAM( V ))) //
     &                  '" from file ' // ANAME
                        CALL M3EXIT( 'CSGLMERGE', JDATE,JTIME, MESG, 2 )
                        GO TO  198

                    END IF          !  if read3() failed
                END IF          !  if aflag:  merging area emissions

                IF ( PFLAG ) THEN
                    IF( .NOT. READ3( PNAME, EINAM( V ), ALLAYS3, 
     &                         JDATE, JTIME, POINTV ) ) THEN

                        MESG = 'Could not read "' //
     &                  EINAM( V )( 1:TRIMLEN( EINAM( V ))) //
     &                  '" from ' // PNAME
                        CALL M3EXIT( 'CSGLMERGE', JDATE,JTIME, MESG, 2 )
                        GO TO  198

                    END IF          !  if read3() failed
                END IF          !  if pflag

                DO  197  J = 1,  NSPCOUT( V )

                    K = SPCREF( V,J )

C...................   Initialization of emisg[]:  either read from bio file or 
C...................   set to zero:

                    IF ( BSPCFLAG( K ) ) THEN
                        IF( .NOT. READ3( BNAME, EMNAM( K ), ALLAYS3,
     &                                   JDATE, JTIME, EMISG ) ) THEN
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
                        M = MBLSCNT( K ) * NVTYPE  
                        N = MBLSDEX( K )
                        CALL MSMATVEC( NMSRC, M, SM( 1,1,N ), 
     &                                 MOBILV, MOBILC )
                        CALL AMRGVEC ( NMSRC, MOBILC, NMATX, 
     &                                 NM, IM, MM,
     &                                 NGRID, EMISG )

                    END IF      !  if mflag


                    IF ( PFLAG ) THEN

                        CALL DMATVEC( NPSRC, SP( 1,K ), POINTV, POINTC )
                        CALL PLMVEC ( FRACL, POINTC, NP, IP, EMISG )

                    END IF      !  if pflag

                    IF ( OUTFLAG ) THEN

                        IF ( .NOT. WRITE3( ENAME, EMNAM( K ),
     &                                     JDATE, JTIME, 
     &                                     EMISG ) ) THEN

                            CALL M3WARN( 'CSGLMERGE', JDATE, JTIME,
     &                      'Could not write "' // 
     &                      EMNAM( K )( 1:TRIMLEN( EMNAM( K ) ) ) //
     &                      '" to ' // ENAME )

                        END IF          !  if write3() failed

                    END IF          !  if outputflag

197             CONTINUE    !  end loop on output vbles for this input vble

198         CONTINUE            ! end loop on input variables

            CALL NEXTIME( JDATE, JTIME, TSTEP )

199     CONTINUE          !  end loop on time steps

999     CONTINUE          !  exit program
      
        CALL M3EXIT( 'CSGLMERGE', 0, 0, 
     &               'Normal completion  of PROGRAM CSGLMERGE', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Error and warning message formats..... 91xxx


C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )
92010   FORMAT( 10X, A )

C...........   Formatted file I/O formats............ 93xxx

93010   FORMAT( A16 )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I10, :, 2X ) )


C...........   Miscellaneous formats................. 95xxx

95000   FORMAT ( /5X , A , $ )          !  generic prompt format.


        END

