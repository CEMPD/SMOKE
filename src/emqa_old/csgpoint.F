C copied by: vukovich
C origin: csgpoint.F 3.1

        PROGRAM CSGPOINT

C***********************************************************************
C  program body starts at line 167
C
C  DESCRIPTION:
C       Use gridding-transform, model-speciation, and control matrices 
C       from programs GRDPMAT, SPCPMAT, and CTLPMAT to grid, speciate,
C       and control source-level point source emissions files.
C
C  PRECONDITIONS REQUIRED:  
C       M3IO source-level point source emissions in canonical order used
C       by RAWPOINT
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; GETYN, PROMPTMFILE, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  7/95 by CJC.
C       Report added 4/97 by JMV.
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
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
C************************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'TMDIMS3.EXT'   !  grid parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         DSCGRID
        INTEGER         ENVINT
        LOGICAL         ENVYN
        LOGICAL         GETYN
        CHARACTER*10    HHMMSS
        INTEGER         INDEX1
        CHARACTER*14    MMDDYY
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN
        INTEGER         WKDAY
        REAL            YR2DAY

        EXTERNAL        DSCGRID, ENVINT, ENVYN, GETYN, HHMMSS, INDEX1,
     &                  MMDDYY, PROMPTFFILE, PROMPTMFILE, TRIMLEN, 
     &                  WKDAY, YR2DAY

C...........   PARAMETER

        CHARACTER*16    NONE
        PARAMETER     ( NONE = 'NONE' )
        REAL GM2TON
        PARAMETER ( GM2TON = 1.102311E-06 )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Point Sources input and output arrays (one variable at a time)

        REAL            EMISV( NPSRC )  !  input emissions.
        REAL            EMISC( NPSRC )  !  controlled speciated emissions
        REAL            EMISG( NGRID )  !  gridded output emissions

C...........   Gridding Matrix (all active coeffs are 1.0, since 
C...........   any point source "lives in" just one cell)

        INTEGER         NS( NGRID )
        INTEGER         IS( NMATX )

        COMMON  / GRIDMAT / NS, IS

C...........   Control Matrix

        REAL            CC( NPSRC, NIPOL )

C...........   Speciation Matrix, and description of how it maps input 
C...........   variables to output variables:  SS( S, SPCREF( V,W ) )
C...........   is used to produce species # OUTREF in the output file
C...........   from species V in the input file
C...........   (NOTE that names in the speciation-matrix file are of the 
C...........   form "<INV-SPC-NAME>to<MODEL-SPC-NAME>")

        REAL            SS( NPSRC, NMPOL )
        INTEGER         NSPCOUT( NIPOL )        !  # of active rows in SPCREF
        INTEGER         SPCREF ( NIPOL, NMPOL )

        CHARACTER*16    INVNAME( NMPOL )	!  inv-spc part of name

C...........   Inventory arrays
        INTEGER         IFIP  ( NPSRC )   ! FIPs  code for each source
        INTEGER         INVYR ( NPSRC )   ! Inventory year code for each source
        INTEGER         ISPTR ( NPSRC )   ! Pointer for each state
        INTEGER         STATES( NPSID )   ! Unique list of state IDs

C...........   Report arrays
        REAL            EMSUM( NPSID, NMPOL ) ! Summed emissions by state

C...........   Other local variables

        REAL            FAC     !  temporary conversion factor

        INTEGER         BTIME   !  beginning time of day in run
        INTEGER         IOS
        INTEGER         NSTA               ! actual number of states

        INTEGER         INO, INO2          ! position in list of pollutants
        INTEGER         STA                ! temporay state code

        INTEGER         LDATE, LTIME, JDATE, JTIME, TSTEP
        INTEGER         RDEV               ! output report unit number
        INTEGER         RTIME              ! time for reporting state totals
        INTEGER         LSTA               ! previous state code

        INTEGER         NSTEPS
        INTEGER         LDEV
        INTEGER         I, J, K, F, S, T, V

        LOGICAL         CONTROLFLAG, OUTPUTFLAG, ST_REPORT

        CHARACTER*16    PNAME   !  logical name for point-source      input file
        CHARACTER*16    GNAME   !  logical name for grid       matrix input file
        CHARACTER*16    CNAME   !  logical name for control    matrix input file
        CHARACTER*16    SNAME   !  logical name for speciation matrix input file
        CHARACTER*16    ENAME   !  logical name for emission output file
        CHARACTER*16    INAME   !  logical name for inventory input file

        CHARACTER*8     REPNAM  !  contains default logical report output name
        CHARACTER*5     OUTNAM  !  contains default logical output name
        CHARACTER*16    SCRBUF  !  scratch buffer
        CHARACTER*16    UNITS   !  temporary units buffer
        CHARACTER*256   MESG    !  message buffer


C***********************************************************************
C   begin body of program CSGPOINT

        LDEV = INIT3()

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program CSGPOINT to take a sorted source level point source ',
     &  'emissions file and a gridding-transform matrix computed by',
     &  'program GRDPMAT, a speciation matrix computed by SPCPMAT, ',
     &  'optionally a control matrix produced by CTLPMAT, and produce ',
     &  'gridded emissions.  The output time step structure of the ',
     &  'gridded file will be the same as the time step structure of ',
     &  'the input point source emissions file.',
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
            CALL M3EXIT( 'CSGPOINT', 0, 0, 'Ending Program', 2 )
        END IF

C.......   Get environment variable settings
        ST_REPORT = ENVYN( 'ST_REPORT', 'Per-state report writing flag',
     &                      .TRUE., IOS )

        RTIME = ENVINT( 'ST_REPTIME', 
     &                  'Per-State summation reporting time (HHMMSS)',
     &                   230000, IOS )

C.......   Get file names; open input gridding, speciation, 
C.......   optional control matrix files, point source file, 
C.......   optional output file

        GNAME = PROMPTMFILE( 
     &          'Enter logical name for GRID MATRIX file',
     &          FSREAD3, 'PGMAT', 'CSGPOINT' )

        SNAME = PROMPTMFILE( 
     &          'Enter logical name for SPECIATION MATRIX file',
     &          FSREAD3, 'PSMAT', 'CSGPOINT' )

        CNAME = PROMPTMFILE( 
     &          'Enter logical name for CONTROL MATRIX file ' //
     &          'or "NONE"',
     &          FSREAD3, 'PCMAT', 'CSGPOINT' )

        CONTROLFLAG = ( CNAME .NE. NONE )

        PNAME = PROMPTMFILE( 
     &          'Enter logical name for POINT SOURCE file',
     &          FSREAD3, 'PTMP', 'CSGPOINT' )

C.......   Build description of output file, and optionally open it:

        IF ( .NOT. DESC3( PNAME ) ) THEN
            CALL M3EXIT( 'CSGPOINT', 0, 0, 
     &          'Could not get description of file ' // PNAME, 2 )

        ELSE IF ( NROWS3D .NE. NPSRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  EMISSIONS file:', NROWS3D,
     &      'program:', NPSRC
            CALL M3EXIT( 'CSGPOINT', 0, 0, MESG, 2 )

        END IF

        JDATE = SDATE3D
        JTIME = STIME3D
        TSTEP = TSTEP3D

C.........  Annual inventory file as input with controls
        IF ( TSTEP .EQ. 0 .AND. CONTROLFLAG ) THEN
            NSTEPS  = 1
            FDESC3D( 1 ) = 'Point gridded, controlled model species'
            OUTNAM = 'PGSC'
            REPNAM = 'PGSCREP'
            UNITS  = 'mole/year'

C.........  Annual inventory file as input with NO controls
        ELSEIF( TSTEP .EQ. 0 ) THEN
            NSTEPS  = 1
            FDESC3D( 1 ) = 'Point gridded model species'
            OUTNAM = 'PGS'
            REPNAM = 'PGSREP'
            UNITS  = 'mole/year'

C.........  Hourly inventory file as input with controls
        ELSEIF ( TSTEP .GT. 0 .AND. CONTROLFLAG ) THEN
            NSTEPS  = MXREC3D
            FDESC3D( 1 ) = 
     &        'Point gridded, hourly, controlled model species'
            OUTNAM = 'PGTSC'
            REPNAM = 'PGTSCREP'
            UNITS  = 'mole/hour'

C.........  Hourly inventory file as input with NO controls
        ELSE
            NSTEPS  = MXREC3D
            FDESC3D( 1 ) = 'Point gridded, hourly model species'
            OUTNAM = 'PGTS'
            REPNAM = 'PGTSREP'
            UNITS  = 'mole/hour'

        END IF

C.......   Read FIPS code for reporting purpose
        IF( ST_REPORT ) THEN

C............  If we've read in temporal file, must get inventory file name
            IF ( TSTEP .NE. 0 ) THEN

                INAME = PROMPTMFILE(
     &                  'Enter logical name for POINT INVENTORY file',
     &                  FSREAD3, 'PNTS', 'CSGPOINT' )
            ELSE
                INAME = PNAME

            ENDIF

            IF( .NOT. READ3( INAME, 'IFIP', ALLAYS3, 0, 0, IFIP ) ) THEN

                CALL M3EXIT( 'CSGPOINT', 0, 0,
     &                       'Error reading FIP from file' // PNAME, 2 )

            END IF

            IF( .NOT.
     &          READ3( INAME, 'INVYR', ALLAYS3, 0, 0, INVYR ) ) THEN

                CALL M3EXIT( 'CSGPOINT', 0, 0,
     &               'Error reading INVYR from file' // PNAME, 2 )

            END IF

            RDEV  = PROMPTFFILE(
     &             'Enter logical name for TOTALS REPORT or "NONE"',
     &             .FALSE., .TRUE., REPNAM, 'CSGPOINT' )
C............  Create state codes arrays
            LSTA = -9
            K    = 0
            DO 22 I = 1, NPSRC

                STA = IFIP( I ) / 1000

                IF( STA .GT. LSTA ) THEN

                    K = K + 1
                    STATES( K ) = STA
                    LSTA        = STA

                ENDIF

                ISPTR ( I ) = K

22          CONTINUE

            NSTA = K

            IF( NSTA .GT. NPSID ) THEN

                WRITE( MESG,94010 )
     &                'Number of states is ', NSTA,
     &                'but maximum (NPSID) is', NPSID
                CALL M3EXIT( 'CSGPOINT', JDATE, JTIME, MESG, 2 )
            ENDIF

C............  Initialize emission summary array
            DO 28 I = 1, NPSID
                DO 27 S = 1, NMPOL
                    EMSUM( I,S ) = 0.
27              CONTINUE
28          CONTINUE

C............  Set position of "special" species in name list
            INO  = INDEX1( 'NO' , NMPOL, EMNAM )
            INO2 = INDEX1( 'NO2', NMPOL, EMNAM )

        ENDIF

        IF ( .NOT. DSCGRID( GRDNM, SCRBUF,
     &                      GDTYP3D, P_ALP3D, P_BET3D, P_GAM3D,
     &                      XCENT3D, YCENT3D, XORIG3D, YORIG3D,
     &                      XCELL3D, YCELL3D,
     &                      NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            CALL M3EXIT( 'CSGPOINT', 0, 0,
     &              'Could not get description of grid ' // GRDNM, 2 )

        END IF

        NVARS3D = NMPOL
        FTYPE3D = GRDDED3   !  shares most of file-description with input file.
        NLAYS3D = 1
        GDNAM3D = GRDNM
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3

        DO  33 V = 1, NMPOL
            VNAME3D( V ) = EMNAM( V )
            UNITS3D( V ) = UNITS
            VDESC3D( V ) = 'gridded speciated emissions'
            VTYPE3D( V ) = m3REAL
33      CONTINUE


        ENAME = PROMPTMFILE( 
     &          'Enter logical name for OUTPUT EMIS file ' //
     &          'or "NONE"',
     &          FSUNKN3, OUTNAM, 'CSGPOINT' )

        OUTPUTFLAG = ( ENAME .NE. NONE )


C.......   Read gridding, speciation, and (optional) control matrices:

        CALL M3MSG2( 'Reading GRIDDING matrix...' )

        IF ( .NOT. READ3( GNAME, 'ALL', 1, 0, 0, NS ) ) THEN

            CALL M3ERR( 'CSGPOINT', 0, 0,
     &           'Could not read gridding matrix from file "' //
     &           GNAME( 1 : TRIMLEN( GNAME ) ) // '".',
     &           .TRUE. )

        END IF      !  if read3() failed for gridding matrix

        CALL M3MSG2( 'Reading SPECIATION matrix...' )

        IF ( .NOT. READ3( SNAME, 'ALL', 1, 0, 0, SS ) ) THEN

            CALL M3ERR( 'CSGPOINT', 0, 0,
     &           'Could not read SPECIATION matrix from file "'//
     &           SNAME( 1 : TRIMLEN( SNAME ) ) // '".',
     &           .TRUE. )

        END IF      !  if read3() failed for speciation matrix


C.......   Build description of how to use speciation matrix:

        IF ( .NOT. DESC3( SNAME ) ) THEN
            CALL M3EXIT( 'CSGPOINT', 0, 0, 
     &                   'Could not get description of file ' // SNAME,
     &                  .TRUE. )
        ELSE IF ( NROWS3D .NE. NPSRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  S-MATRIX file:', NROWS3D,
     &      'program:', NPSRC
            CALL M3EXIT( 'CSGPOINT', 0, 0, MESG, 2 )
        END IF

        DO  54  V = 1, NIPOL
            NSPCOUT( V ) = 0
54      CONTINUE

        DO  55  V = 1, NMPOL

            K = INDEX( VNAME3D( V ), '_to_' )
            IF( K .LT. 1 ) THEN
                CALL M3EXIT( 'CSGPOINT', 0, 0,
     &              'Bad conversion factor name "' // 
     &              VNAME3D( V )( 1:TRIMLEN( VNAME3D( V ) ) ) //
     &              '" in ' // SNAME, 2 )
            END IF
            INVNAME( V ) = VNAME3D( V )( 1:K-1 )

            I = INDEX1( INVNAME( V ), NIPOL, EINAM )
            IF ( I .LT. 0 ) THEN
                CALL M3ERR( 'CSGPOINT', 0, 0, 
     &                  'Bad input species "' // 
     &                  INVNAME( V )( 1:TRIMLEN( INVNAME( V ) ) ) // 
     &                  '" in ' // SNAME, 
     &                  .TRUE. )
            END IF

            SCRBUF = VNAME3D( V )( K+4:16 )
            J = INDEX1( SCRBUF, NMPOL, EMNAM )
            IF ( J .LT. 0 ) THEN
                CALL M3ERR( 'CSGPOINT', 0, 0, 
     &                      'Bad output species "' // 
     &                      SCRBUF( 1:TRIMLEN( SCRBUF ) ) // 
     &                      '" in ' // SNAME, 2 )
            END IF

            F = NSPCOUT( I ) + 1
            SPCREF ( I,F ) = J
            NSPCOUT( I )   = F
            
55      CONTINUE


C.......   If CONTROLFLAG, read and apply control matrix:

        IF( CONTROLFLAG ) THEN

            CALL M3MSG2( 'Reading CONTROL matrix...' )

            IF ( .NOT. DESC3( CNAME ) ) THEN
                WRITE( MESG, '(2A)' )
     &          'Could not get description of file ', CNAME
                CALL M3EXIT( 'CSGPOINT', 0, 0, MESG, 2 )
            ELSE IF ( NROWS3D .NE. NPSRC ) THEN
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  C-MATRIX file:', NROWS3D,
     &          'program:', NPSRC
                CALL M3EXIT( 'CSGPOINT', 0, 0, MESG, 2 )
            END IF

            DO  66 S = 1, NIPOL

                SCRBUF = 'CONTROL_' // EINAM( S )
                IF ( .NOT. READ3( CNAME, SCRBUF, 1, 0, 0, CC( 1,S ) ) 
     &             ) THEN

                    CALL M3ERR( 'CSGPOINT', 0, 0,
     &                   'Could not read "' // SCRBUF( 1:12 ) // 
     &                   '" from CONTROL FILE "' //
     &                   CNAME( 1 : TRIMLEN( CNAME ) ) // '".',
     &                   .TRUE. )

                END IF      !  if read3() failed for control matrix

                DO  65 V = 1, NMPOL
                    IF ( EINAM( S ) .EQ. INVNAME( V ) ) THEN
                        DO  64  I = 1, NPSRC
                            SS( I,V ) = CC( I,S ) * SS( I,V )
64                      CONTINUE
                    END IF
65              CONTINUE

66          CONTINUE

        END IF          !  if controlflag


C.......   Transform and write out point source emissions values:

        IF ( CONTROLFLAG .AND. OUTPUTFLAG ) THEN
            CALL M3MSG2( 'Calculating and writing ' //
     &                   'controlled speciated gridded emissions...' )
        ELSE IF ( CONTROLFLAG ) THEN
            CALL M3MSG2( 'Calculating controlled ' //
     &                   'speciated gridded emissions...' )
        ELSE IF ( OUTPUTFLAG ) THEN
            CALL M3MSG2( 'Calculating and writing ' //
     &                   'speciated gridded emissions...' )
        ELSE 
            CALL M3MSG2( 'Calculating speciated gridded emissions...' )
        END IF

        BTIME = JTIME
        LDATE = 0

        DO  199  T = 1, NSTEPS

C...............   If this is a new month, or new day, write message
            IF ( JDATE .GT. 0 .AND. LDATE .NE. JDATE ) THEN
                LDATE=JDATE 
                MESG = 'Processing ' //
     &                 DAYS( WKDAY( JDATE ) ) // MMDDYY( JDATE )
                CALL M3MSG2( MESG( 1:TRIMLEN( MESG ) ) )
 
            END IF

C.............  Write to screen because WRITE3 only writes to LDEV
            IF( JDATE .GT. 0 ) WRITE( *, 93020 ) HHMMSS( JTIME )

            DO  198  V = 1, NIPOL

                IF( READ3( PNAME, EINAM( V ), ALLAYS3, 
     &                     JDATE, JTIME, EMISV ) ) THEN

                    DO  197  J = 1,  NSPCOUT( V )

                        K = SPCREF( V,J )

                        CALL DMATVEC( NPSRC, SS( 1,K ), EMISV, EMISC )

                        CALL PMATVEC( NPSRC, NGRID, NMATX, NS, IS,
     &                                EMISC, EMISG )

                        IF ( OUTPUTFLAG ) THEN

                            IF ( .NOT. WRITE3( ENAME, EMNAM( K ),
     &                                         JDATE, JTIME, 
     &                                         EMISG ) ) THEN

                                CALL M3EXIT( 'CSGPOINT', JDATE, JTIME,
     &                          'Could not write "' // 
     &                          EMNAM( K )( 1: TRIMLEN( EMNAM( K ))) // 
     &                          '" to ' // ENAME, 2 )
      
                            END IF          !  if write3() failed

                        END IF          !  if outputflag

C.........................  When reporting is switched on, calculate emissions
C.........................  summaries by state and report to output QA file
C.........................  EMSUM must be initialized
                        IF( ST_REPORT ) THEN

C.............................  For annual values, convert to day-specific
                            IF( TSTEP .EQ. 0 ) THEN

                                DO 111 S = 1, NPSRC
                                    EMISC( S ) = YR2DAY( INVYR( S ) ) *
     &                                           EMISC( S )
111                             CONTINUE

                            ENDIF

                            CALL PGRDSUM( NPSRC, NGRID, NMATX, NS, IS,
     &                                    ISPTR, EMISC,
     &                                    EMSUM( 1,K ) )
                        ENDIF

197                 CONTINUE    !  end loop on output vbles for this input vble

                ELSE

                    CALL M3EXIT( 'CSGPOINT', JDATE, JTIME, 
     &                           'Could not read "' //  
     &                           EINAM( V )( 1:TRIMLEN( EINAM( K ))) //
     &                           '" from ' // PNAME, 2 )

                END IF          !  if read3() succeeds or not

198         CONTINUE            ! end loop on input variables

            LDATE = JDATE
            LTIME = JTIME

            CALL NEXTIME( JDATE, JTIME, TSTEP )

C.............  Write report if report flag is switched on and previous
C.............  time equals the report time

            IF( ST_REPORT .AND.
     &          ( LTIME .EQ. RTIME .OR. T .EQ. NSTEPS ) ) THEN

C.................  Convert summary emissions from moles/day to tons/day

                DO 155 S = 1, NMPOL

C.....................  Conversion factor for AREA as input file
                    IF( S .NE. INO ) THEN
                        FAC = EMWT( S ) * GM2TON
                    ELSE
                        FAC = EMWT( INO2 ) * GM2TON
                    ENDIF

                    DO 144 I = 1, NSTA

                       EMSUM( I,S ) = EMSUM( I,S ) * FAC

144                 CONTINUE

155             CONTINUE

                IF( TSTEP .NE. 0 .AND. RDEV .GT. 0)
     &          WRITE( RDEV, 94010 ) DAYS( WKDAY( LDATE ) ) //
     &                               MMDDYY( LDATE ) // ' from',
     &                               BTIME, 'to', LTIME

                IF( TSTEP .EQ. 0 ) THEN
                    MESG = 'Point source emissions ' //
     &                     '(tons/ave inventory day) within grid ' //
     &                     GRDNM
                ELSE
                    MESG = 'Point source emissions ' //
     &                     '(tons/specific day) within grid ' // GRDNM
                ENDIF
                IF( RDEV .GT. 0) 
     &          CALL STATELST( RDEV  , MESG , NPSID , NSTA , NMPOL,
     &                         'F7.1', EMNAM, STATES, EMSUM         )

C.................  Reinitialize emissions summary
                DO 298 I = 1, NPSID
                    DO 297 S = 1, NMPOL
                        EMSUM( I,S ) = 0.
297                 CONTINUE
298             CONTINUE

C.................  Set begining time for next day
                BTIME = JTIME

            ENDIF
 

199     CONTINUE          !  end loop on time steps

999     CONTINUE          !  exit program

        CALL M3EXIT( 'CSGPOINT', 0, 0, 
     &               'Normal completion  of PROGRAM CSGPOINT', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )


C...........   Formatted file I/O formats............ 93xxx

93010   FORMAT( A16 )

93020   FORMAT( 8X, 'at time ', A8 )
 

C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT( 10( A, :, I7, :, 1X ) )

 
C...........   Miscellaneous formats................. 95xxx

95000   FORMAT ( /5X , A , $ )          !  generic prompt format.


        END

