
        PROGRAM GRDAREA

C***********************************************************************
C  program body starts at line  113
C
C  DESCRIPTION:
C       Use gridding-transform matrices from program GRDAMAT to grid
C       source-level area source emissions files.
C
C  PRECONDITIONS REQUIRED:  
C       M3IO source-level area source emissions in canonical order used
C       by GRDAMAT
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; GETYN, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  2/95 by CJC.
C
C       Revised 9/96 by CJC:  match with char-string pollutant names,
C	rather than numeric names, in loop 33
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE)
C                Modeling System
C File Version @(#)$Id$
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
C Pathname:    $Source$
C Last updated: $Date$ 
C
C***************************************************************************/

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'ARDIMS3.EXT'   !  area-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'TMDIMS3.EXT'   !  temporal parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         DSCGRID
        LOGICAL         GETYN
        CHARACTER*10    HHMMSS
        INTEGER         INDEX1
        CHARACTER*14    MMDDYY
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN
        INTEGER         WKDAY

        EXTERNAL        DSCGRID, GETYN, HHMMSS, INDEX1, MMDDYY,
     &                  PROMPTMFILE, TRIMLEN, WKDAY


C...........   LOCAL VARIABLES and their descriptions:
C...........   Area Sources input and output arrays

        REAL            EMISV( NASRC )  !  one variable of input emissions.
        REAL            EMISG( NGRID )  !  gridded output variable

C...........   Gridding Matrix

        INTEGER     NS( NGRID )
        INTEGER     IS( NMATX )
        REAL        CS( NMATX )

        COMMON  / GRIDMAT / NS, IS, CS

C...........   Other local variables

        INTEGER         L, T, V
        INTEGER         LDATE, JDATE, JTIME, TSTEP
        INTEGER         LDEV
        INTEGER         NSTEPS
        INTEGER         VTYPE( NIPOL )

        LOGICAL         EFLAG   !  Error flag
 
        CHARACTER*4     OUTNAM  !  contains default  logical output name
        CHARACTER*16    ANAME   !  logical name for area-source      input file
        CHARACTER*16    GNAME   !  logical name for grid matrix      input file
        CHARACTER*16    ENAME   !  logical name for emission source output file

        CHARACTER*16    SCRBUF  !  scratch name buffer
        CHARACTER*16    VNAME( NIPOL )
        CHARACTER*16    UNITS( NIPOL )
        CHARACTER*16    VDESC( NIPOL )
        CHARACTER*256   MESG


C***********************************************************************
C   begin body of program GRDAREA

        LDEV = INIT3()

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program GRDAREA to take a sorted source level area source ',
     &  'emissions file and a gridding-transform matrix computed by',
     &  'program GRDAMAT, and produce gridded emissions.  You will be',
     &  'able to choose which variables are gridded; the output time',
     &  'step structure of the gridded file will be the same as the',
     &  'time step structure of the input area source emissions file.',
     &  ' '
        WRITE( *,92000 ) 
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").  Input files must',
     &  'have been sorted as indicated, prior to program execution.',
     &  'You will be prompted for which species to grid (and do not',
     &  'need to grid them all; in particular, you SHOULD NOT grid',
     &  'the FIP-related or ASC-related variables.)',
     &  ' '
        WRITE( *,92000 ) 
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry.  Default responses are indicated',
     &  'in brackets [LIKE THIS].',
     &  ' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'GRDAREA', 0, 0 , 'Ending program.', 2 )
        END IF

        GNAME = PROMPTMFILE( 
     &          'Enter logical name for input GRID MATRIX file',
     &          FSREAD3, 'AGMAT', 'GRDAREA' )

        ANAME = PROMPTMFILE( 
     &          'Enter logical name for ANNUAL or TIME-STEPPED ' //
     &          'AREA SOURCE file',
     &          FSREAD3, 'ATMP', 'GRDAREA' )


C.......   Borrow much of file description from input area file.
C.......   Get list of output variables (shrink list in file description
C.......   data structure to just those variables desired):

        IF ( .NOT. DESC3( ANAME ) ) THEN
            CALL M3ERR( 'GRDAREA', 0, 0, 
     &                  'Could not get description of file ' // ANAME,
     &                  .TRUE. )

        ELSEIF( NROWS3D .NE. NASRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  AREA SOURCES file:', NROWS3D,
     &      'program:', NASRC

        END IF


        IF ( .NOT.DSCGRID( GRDNM, SCRBUF, GDTYP3D, 
     &                     P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D, 
     &                     XORIG3D, YORIG3D, XCELL3D, YCELL3D, 
     &                     NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            CALL M3EXIT( 'CSGAREA', 0, 0,
     &              'Could not get description of grid ' // GRDNM, 2 )

        END IF

C......... Store variable properties based on DESC3 call
        EFLAG = .FALSE.
        DO  22 V = 1, NIPOL
            L = INDEX1( EINAM( V ), NVARS3D, VNAME3D )
 
            IF( L .LE. 0 ) THEN
 
                EFLAG = .TRUE.
                MESG = 'ERROR: Inventory pollutant ' // EINAM( V ) //
     &                 ' not found in input file!'
                CALL M3MSG2( MESG )
 
            ELSE
                VNAME( V ) = VNAME3D( L )
                UNITS( V ) = UNITS3D( L )
                VTYPE( V ) = VTYPE3D( L )
                VDESC( V ) = VDESC3D( L )
            ENDIF
 
22      CONTINUE

        IF( EFLAG ) THEN
            CALL M3EXIT( 'GRDAREA', 0, 0, 'Bad input file', 2 )
 
        ENDIF

        NVARS3D = NIPOL
        FTYPE3D = GRDDED3   !  shares most of file-description with input file.
        NLAYS3D = 1
        GDNAM3D = GRDNM
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3
        JDATE   = SDATE3D
        JTIME   = STIME3D
        TSTEP   = TSTEP3D

C.........  Annual inventory file as input
        IF ( TSTEP .EQ. 0 ) THEN
            NSTEPS  = 1
            FDESC3D( 1 ) = 'Area gridded inventory pollutants'
            OUTNAM = 'AGRD'
 
C.........  Hourly inventory file as input
        ELSE
            NSTEPS  = MXREC3D
            FDESC3D( 1 ) = 'Area gridded hourly inventory pollutants'
            OUTNAM = 'AGRT'
 
        END IF

C.......... Set up output variable characteristics
        DO  33 V = 1, NIPOL
            VNAME3D( V ) = VNAME( V )
            UNITS3D( V ) = UNITS( V )
            VTYPE3D( V ) = VTYPE( V )
            VDESC3D( V ) = VDESC( V )
33      CONTINUE

C.......   Get file name; open output gridded area sources file

        ENAME = PROMPTMFILE( 
     &          'Enter logical name for OUTPUT EMIS file',
     &          FSUNKN3, OUTNAM, 'GRDAREA' )


C.......   Read gridding-transform matrix:

        CALL M3MSG2( 'Reading GRIDDING matrix...' )

        IF ( .NOT. READ3( GNAME, 'ALL', 1, 0, 0, NS ) ) THEN

            CALL M3ERR( 'GRDAREA', 0, 0,
     &                  'Could not read gridding matrix from file "' //
     &                  GNAME( 1 : TRIMLEN( GNAME ) ) // '".',
     &                  .TRUE. )

        END IF

C.......   Write out the area source emissions values:

        CALL M3MSG2( 'Calculating gridded area source emissions...' )

        LDATE = 0
        DO  199  T = 1, NSTEPS

C...............   If this is a new month, or new day, write message
            IF ( JDATE .GT. 0 .AND. LDATE .NE. JDATE  ) THEN
 
                MESG = 'Processing ' //
     &                 DAYS( WKDAY( JDATE ) ) // MMDDYY( JDATE )
                CALL M3MSG2( MESG( 1:TRIMLEN( MESG ) ) )
 
                LDATE = JDATE

            END IF

C.............  Write to screen because WRITE3 only writes to LDEV
            IF( JDATE .GT. 0 ) WRITE( *, 93020 ) HHMMSS( JTIME )

            DO  198  V = 1, NIPOL

                IF( READ3( ANAME, EINAM( V ), ALLAYS3, 
     &                     JDATE, JTIME, EMISV ) ) THEN

                    CALL SMATVEC( NASRC, NGRID, NMATX, NS, IS, CS,
     &                            EMISV, EMISG )

                    IF ( .NOT. WRITE3( ENAME, EINAM( V ), 
     &                                 JDATE, JTIME, EMISG ) ) THEN

                        CALL M3EXIT( 'GRDAREA', JDATE, JTIME,
     &                       'Could not write "' // 
     &                       EINAM( V )( 1: TRIMLEN( EINAM( V ))) //
     &                      '" to ' // ENAME, 2 )

                    END IF

                ELSE

                    CALL M3EXIT( 'GRDAREA', JDATE, JTIME, 
     &                           'Could not read "' // 
     &                           EINAM( V )( 1:TRIMLEN( EINAM( V ))) //
     &                           '" from ' // ANAME, 2 )

                END IF  !  if read3() succeeds or not

198         CONTINUE

            CALL NEXTIME( JDATE, JTIME, TSTEP )

199     CONTINUE

999   CONTINUE          !  exit program

      CALL M3EXIT( 'GRDAREA', 0, 0, 'Normal completion', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Error and warning message formats..... 91xxx


C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )


C...........   Formatted file I/O formats............ 93xxx

93010   FORMAT( A16 )

93020   FORMAT( 8X, 'at time ', A8 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT ( 10 ( A, :, I10, :, 2X ) )

C...........   Miscellaneous formats................. 95xxx

95000   FORMAT ( /5X , A , $ )          !  generic prompt format.


        END

