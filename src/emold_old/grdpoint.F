
        PROGRAM GRDPOINT

C***********************************************************************
C  program body starts at line 105
C
C  DESCRIPTION:
C       Use gridding-transform matrices from program GRDPMAT to grid
C       source-level point source emissions files.
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
C****************************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'TMDIMS3.EXT'   !  temporal parameters (for DAYS variable)
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


C...........   PARAMETER

        CHARACTER*16    NONE
        PARAMETER     ( NONE = 'NONE' )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Point Sources input and output arrays (one variable at a time)

        REAL            EMISV( NPSRC )  !  input emissions.
        REAL            EMISG( NGRID )  !  gridded output emissions

C...........   Gridding Matrix (all active coeffs are 1.0, since 
C...........   any point source "lives in" just one cell)

        INTEGER         NS( NGRID )
        INTEGER         IS( NMATX )

        COMMON  / GRIDMAT / NS, IS

C...........   Other local variables

        INTEGER         I, J, L, K, F, S, T, V
        INTEGER         JDATE, JTIME, TSTEP
        INTEGER         LDEV
        INTEGER         NSTEPS
        INTEGER         VTYPE( NIPOL )

        LOGICAL         EFLAG   !  Error flag

        CHARACTER*4     OUTNAM  !  contains default  logical output name
        CHARACTER*16    PNAME   !  logical name for point-source      input file
        CHARACTER*16    GNAME   !  logical name for grid       matrix input file
        CHARACTER*16    ENAME   !  logical name for emission output file

        CHARACTER*16    SCRBUF  !  scratch name buffer
        CHARACTER*16    VNAME( NIPOL )
        CHARACTER*16    UNITS( NIPOL )
        CHARACTER*16    VDESC( NIPOL )
        CHARACTER*256   MESG


C***********************************************************************
C   begin body of program GRDPOINT

        LDEV = INIT3()

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program GRDPOINT to take a sorted source level point source ',
     &  'emissions file and a gridding-transform matrix computed by',
     &  'program GRDPMAT and produce gridded emissions.  The output ',
     &  'time step structure of the gridded file will be the same as ',
     &  'the time step structure of the input point source emissions',
     &  'file.',
     &  ' '
        WRITE( *,92000 ) 
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").',
     &  ' ',
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry.  Default responses are indicated',
     &  'in brackets [LIKE THIS].',
     &  ' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            WRITE( *,92000 ) ' ', 'Ending program GRDPOINT', ' '
            CALL EXIT( 0 )
        END IF


C.......   Get file names; open input gridding, speciation, 
C.......   optional control matrix files, point source file, 
C.......   optional output file

        GNAME = PROMPTMFILE( 
     &          'Enter logical name for GRID MATRIX file',
     &          FSREAD3, 'PGMAT', 'GRDPOINT' )

        PNAME = PROMPTMFILE( 
     &          'Enter logical name for ANNUAL or TIME-STEPPED ' //
     &          'POINT SOURCE file',
     &          FSREAD3, 'PTMP', 'GRDPOINT' )

C.......   Build description of output file, and optionally open it:

        IF ( .NOT. DESC3( PNAME ) ) THEN
            CALL M3EXIT( 'GRDPOINT', 0, 0, 
     &          'Could not get description of file ' // PNAME, 2 )

        ELSEIF( NROWS3D .NE. NPSRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  POINT SOURCES file:', NROWS3D,
     &      'program:', NPSRC
            CALL M3EXIT( 'GRDPOINT', 0, 0, MESG, 2 )

        ENDIF

        IF ( .NOT. DSCGRID( GRDNM, SCRBUF,
     &                      GDTYP3D, P_ALP3D, P_BET3D, P_GAM3D, 
     &                      XCENT3D, YCENT3D, XORIG3D, YORIG3D, 
     &                      XCELL3D, YCELL3D, 
     &                      NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            CALL M3EXIT( 'GRDPOINT', 0, 0,
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
            CALL M3EXIT( 'GRDPOINT', 0, 0, 'Bad input file', 2 )

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
            FDESC3D( 1 ) = 'Point gridded inventory pollutants'
            OUTNAM = 'PGRD'

C.........  Hourly inventory file as input
        ELSE
            NSTEPS  = MXREC3D
            FDESC3D( 1 ) = 'Point gridded hourly inventory pollutants'
            OUTNAM = 'PGRT'

        END IF

C.......... Set up output variable characteristics
        DO  33 V = 1, NIPOL
            VNAME3D( V ) = VNAME( V )
            UNITS3D( V ) = UNITS( V )
            VTYPE3D( V ) = VTYPE( V )
            VDESC3D( V ) = VDESC( V )
33      CONTINUE

        ENAME = PROMPTMFILE( 
     &          'Enter logical name for OUTPUT EMIS file',
     &          FSUNKN3, OUTNAM, 'GRDPOINT' )


C.......   Read gridding matrix:

        CALL M3MSG2( 'Reading GRIDDING matrix...' )

        IF ( .NOT. READ3( GNAME, 'ALL', 1, 0, 0, NS ) ) THEN

            CALL M3EXIT( 'GRDPOINT', 0, 0,
     &           'Could not read gridding matrix from file "' //
     &           GNAME( 1 : TRIMLEN( GNAME ) ) // '".',
     &           2 )

        END IF      !  if read3() failed for gridding matrix


C.......   Transform and write out point source emissions values:

        CALL M3MSG2( 'Calculating gridded point source emissions...' )

        DO  199  T = 1, NSTEPS

C...............   If this is a new month, or new day, write message
            IF ( JDATE .GT. 0 .AND. MOD( T,24 ) .EQ. 1 ) THEN
 
                MESG = 'Processing ' //
     &                 DAYS( WKDAY( JDATE ) ) // MMDDYY( JDATE )
                CALL M3MSG2( MESG( 1:TRIMLEN( MESG ) ) )
 
            END IF
 
C.............  Write to screen because WRITE3 only writes to LDEV
            IF( JDATE .GT. 0 ) WRITE( *, 93020 ) HHMMSS( JTIME )

            DO  198  V = 1, NIPOL

                IF( READ3( PNAME, EINAM( V ), ALLAYS3, 
     &                     JDATE, JTIME, EMISV ) ) THEN

                    CALL PMATVEC( NPSRC, NGRID, NMATX, NS, IS,
     &                            EMISV, EMISG )

                    IF ( .NOT. WRITE3( ENAME, EINAM( V ),
     &                                 JDATE, JTIME, 
     &                                 EMISG ) ) THEN
       
                        CALL M3EXIT( 'GRDPOINT', JDATE, JTIME,
     &                      'Could not write "' // 
     &                      EINAM( V )( 1: TRIMLEN( EINAM( V ))) // 
     &                      '" to ' // ENAME, 2 )

                    END IF          !  if write3() failed
      
                ELSE

                    CALL M3EXIT( 'GRDPOINT', JDATE, JTIME, 
     &                           'Could not read "' //  
     &                           EINAM( V )( 1:TRIMLEN( EINAM( V ))) //
     &                           '" from ' // PNAME, 2 )

                END IF          !  if read3() succeeds or not

198         CONTINUE            ! end loop on input variables V

            CALL NEXTIME( JDATE, JTIME, TSTEP )

199     CONTINUE          !  end loop on time steps

999     CONTINUE          !  exit program

        CALL M3EXIT( 'GRDPOINT', 0, 0, 
     &               'Normal completion of PROGRAM GRDPOINT', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

C...........   Formatted file I/O formats............ 93xxx

93020   FORMAT( 8X, 'at time ', A8 )

C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT ( 10 ( A, :, I10, :, 2X ) )

        END

