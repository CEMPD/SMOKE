
        PROGRAM GRDAREA

C***********************************************************************
C  program body starts at line  96
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
C Pathname:    $Source$
C Last updated: $Date$ 
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
C***************************************************************************/

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'ARDIMS3.EXT'   !  area-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         GETYN
        INTEGER         INDEX1
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN

        EXTERNAL    GETYN, INDEX1, PROMPTMFILE, TRIMLEN


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

        INTEGER          JDATE, JTIME, TSTEP
        INTEGER          NSTEPS
        INTEGER          LDEV
        INTEGER          I, J, T, V

        CHARACTER*16    ANAME   !  logical name for area-source      input file
        CHARACTER*16    GNAME   !  logical name for grid matrix      input file
        CHARACTER*16    ENAME   !  logical name for emission source output file

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
            WRITE( *,92000 ) ' ', 'Ending program GRDAREA', ' '
            CALL EXIT( 0 )
        END IF

        GNAME = PROMPTMFILE( 
     &          'Enter logical name for input GRID MATRIX file',
     &          FSREAD3, 'AGMAT', 'GRDAREA' )

        ANAME = PROMPTMFILE( 
     &          'Enter logical name for input AREA SOURCE file',
     &          FSREAD3, 'AREA', 'GRDAREA' )


C.......   Get list of output variables (shrink list in file description
C.......   data structure to just those variables desired):

        IF ( .NOT. DESC3( ANAME ) ) THEN
            CALL M3ERR( 'GRDAREA', 0, 0, 
     &                  'Could not get description of file ' // ANAME,
     &                  .TRUE. )
        END IF


        I = 0
        J = 0

        DO  33  I = 1, NVARS3D

            IF( VTYPE3D( I ) .EQ. M3REAL ) THEN !  only grid REAL vbles 

                IF ( GETYN( 'Grid variable "' // 
     &                      VNAME3D( I )( 1:TRIMLEN( VNAME3D( I ) ) )
     &                      // '" for output?',
     &                      ( 0 .NE. INDEX1( VNAME3D( I ), 
     &                                       NIPOL, EINAM ) ) ) ) THEN

                    J = J + 1
                    VNAME3D( J ) = VNAME3D( I )
                    UNITS3D( J ) = UNITS3D( I )
                    VDESC3D( J ) = VDESC3D( I )
                    VTYPE3D( J ) = M3REAL

                END IF		!  if getyn-prompt accepted

                GO TO  33

            END IF		!  if vtype real

33      CONTINUE

        NVARS3D = J
        FTYPE3D = GRDDED3   !  shares most of file-description with input file.
        SDATE3D = 0
        STIME3D = 0
        TSTEP3D = 0
        NCOLS3D = NCOLS
        NROWS3D = NROWS
        NVARS3D = 1
        NLAYS3D = 1
        NTHIK3D = NASRC
        GDNAM3D = 'NC_UAM'
        GDTYP3D = UTMGRD3
        P_ALP3D = 17.0D0                 !  DBLE( ZONE )
        P_BET3D = DBLE( BADVAL3 )
        P_GAM3D = DBLE( BADVAL3 )
        XCENT3D =       0.0D0
        YCENT3D =       0.0D0
        XCELL3D =    5000.0D0
        YCELL3D =    5000.0D0
        XORIG3D =  350000.0D0
        YORIG3D = 3785000.0D0
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3
        JDATE   = SDATE3D
        JTIME   = STIME3D
        TSTEP   = TSTEP3D
        IF ( TSTEP .EQ. 0 ) THEN
            NSTEPS  = 1
        ELSE
            NSTEPS  = MXREC3D
        END IF


C.......   Get file name; open output gridded area sources file

        ENAME = PROMPTMFILE( 
     &          'Enter logical name for output GRIDDED EMIS file',
     &          FSUNKN3, 'AGRD', 'GRDAREA' )


C.......   Read gridding-transform matrix:

        WRITE( *,92000 ) ' ', 'Reading GRIDDING-TRANSFORM matrix.', ' '
        IF ( .NOT. READ3( GNAME, 'ALL', 1, 0, 0, NS ) ) THEN

            CALL M3ERR( 'GRDAREA', 0, 0,
     &                  'Could not read gridding matrix from file "' //
     &                  GNAME( 1 : TRIMLEN( GNAME ) ) // '".',
     &                  .TRUE. )

        END IF


C.......   Write out the area source emissions values:

        WRITE( *,92000 ) ' ', 'Writing out emissions values...', ' '

        DO  199  T = 1, NSTEPS

            DO  198  V = 1, NVARS3D

                IF( READ3( ANAME, VNAME3D( V ), ALLAYS3, 
     &                     JDATE, JTIME, EMISV ) ) THEN

                    CALL SMATVEC( NASRC, NGRID, NMATX, NS, IS, CS,
     &                            EMISV, EMISG )

                    IF ( .NOT. WRITE3( ENAME, VNAME3D( V ), 
     &                                 JDATE, JTIME, EMISG ) ) THEN

                        CALL M3ERR( 'GRDAREA', JDATE, JTIME,
     &                              'Could not write "' // 
     &                              VNAME3D( V ) // '" to ' // ENAME,
     &                              .FALSE. )

                    END IF

                ELSE

                    CALL M3ERR( 'GRDAREA', JDATE, JTIME, 
     &                          'Could not read "' // 
     &                          VNAME3D( V ) // '" from ' // ANAME,
     &                          .FALSE. )

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


C...........   Internal buffering formats............ 94xxx


C...........   Miscellaneous formats................. 95xxx

95000   FORMAT ( /5X , A , $ )          !  generic prompt format.


        END

