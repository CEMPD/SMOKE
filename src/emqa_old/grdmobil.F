
        PROGRAM GRDMOBIL

C***********************************************************************
C  program body starts at line 103
C
C  DESCRIPTION:
C       Use gridding-transform from program GRDMMAT to grid time-stepped,
C       source level mobile source emission files.
C
C  PRECONDITIONS REQUIRED:  
C       M3IO source-level area source emissions in canonical order used
C       by GRDMMAT
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; PROMPTMFILE
C
C  REVISION  HISTORY:
C       Prototype  2/96 by CJC.
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
        
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'MBDIMS3.EXT'   !  area-source dimensioning parameters
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
C...........   Area Sources input and output arrays (one variable at a time)
        
        REAL            EMISV( NMSRC, NVTYPE )  !  input emissions.
        REAL            EMISG( NGRID, NVTYPE )  !  gridded output emissions
        
C...........   Gridding Matrix, Control Matrix, Speciation Matrix
      
        INTEGER         NX( NGRID )
        INTEGER         IX( NMATX )
        REAL            CX( NMATX )
        
        COMMON  / GRIDMAT / NX, IX, CX
           
C...........   Other local variables
        
        INTEGER         JDATE, JTIME, TSTEP
        INTEGER         NSTEPS
        INTEGER         LDEV
        INTEGER         L, T, V

        CHARACTER*16    ANAME   !  logical name for area-source       input file
        CHARACTER*16    GNAME   !  logical name for grid       matrix input file
        CHARACTER*16    ENAME   !  logical name for emission output file

        CHARACTER*16    SCRBUF  !  scratch name area
        CHARACTER*256   MESG    !  scratch message area


C***********************************************************************
C   begin body of program GRDMOBIL
        
        LDEV = INIT3()
        
        WRITE( *,92000 )  )
     &  ' ',
     &  'Program GRDMOBIL to take a sorted source level mobile source ',
     &  'emissions file and a gridding-transform matrix computed by',
     &  'program GRDMMAT, and produce gridded vehicle-typed, mode',
     &  'specific emissions.  The output time step structure of the ',
     &  'gridded file will be the same as the time step structure of',
     &  'the input area source emissions file.',
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
            CALL M3EXIT( 'GRDMOBIL', 0, 0, 'Ending Program', 2 )
        END IF


C.......   Get file names; open input gridding matrix, mobile source, and
C.......   output files.

        GNAME = PROMPTMFILE( 
     &          'Enter logical name for GRID MATRIX file',
     &          FSREAD3, 'MGMAT', 'GRDMOBIL' )

        ANAME = PROMPTMFILE( 
     &          'Enter logical name for TIME-STEPPED MOBILE EMIS file',
     &          FSREAD3, 'MTMP', 'GRDMOBIL' )

C.......   Build description of output file, and optionally open it:

        IF ( .NOT. DESC3( ANAME ) ) THEN
            CALL M3EXIT( 'GRDMOBIL', 0, 0, 
     &              'Could not get description of file ' // ANAME, 2 )
        END IF

        IF ( .NOT.DSCGRID( GRDNM, SCRBUF, GDTYP3D, 
     &                     P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D, 
     &                     XORIG3D, YORIG3D, XCELL3D, YCELL3D, 
     &                     NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            CALL M3EXIT( 'GRDMOBIL', 0, 0,
     &              'Could not get description of grid ' // GRDNM, 2 )

        END IF

        FTYPE3D = GRDDED3   !  shares most of file-description with input file.
        GDNAM3D = GRDNM
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

        FDESC3D( 1 ) = 'mobile source gridded VMT'
                 
        
        ENAME = PROMPTMFILE( 
     &          'Enter logical name for OUTPUT EMIS file ',
     &          FSUNKN3, 'MGTS', 'GRDMOBIL' )


C.......   Read gridding, speciation, and (optional) control matrices:

        CALL M3MSG2( 'Reading GRIDDING matrix...' )
        IF ( .NOT. READ3( GNAME, 'ALL', 1, 0, 0, NX ) ) THEN

            CALL M3EXIT( 'GRDMOBIL', 0, 0,
     &           'Could not read gridding matrix from file "' //
     &           GNAME( 1 : TRIMLEN( GNAME ) ) // '".', 2 )
            
        END IF      !  if read3() failed for gridding matrix

C.......   Transform and write out area source emissions values:

        CALL M3MSG2( 'Calculating gridded categorized VMT...' )

        DO  199  T = 1, NSTEPS

            DO  198  V = 1, NVARS3D

C...................   Read in entire emissions field

                IF( .NOT. READ3( ANAME, VNAME3D( V ), ALLAYS3, 
     &                     JDATE, JTIME, EMISV ) ) THEN

                    MESG = 'Could not read "' //
     &                     VNAME3D( V )( 1:TRIMLEN( VNAME3D( V ) ) )//
     &                     '" from ' // ANAME
                    CALL M3EXIT( 'GRDMOBIL', JDATE, JTIME, MESG, 2 )

                END IF          !  if read3() succeeds or not

                DO  197  L = 1, NLAYS3D

                    CALL SMATVEC( NMSRC, NGRID, NMATX, NX, IX, CX,
     &                            EMISV( 1,L ), EMISG( 1,L ) )

197             CONTINUE

                IF ( .NOT. WRITE3( ENAME, VNAME3D( V ),
     &                             JDATE, JTIME, 
     &                             EMISG ) ) THEN

                    MESG = 'Could not write "' //
     &                     VNAME3D( V )( 1:TRIMLEN( VNAME3D( V ) ) ) //
     &                     '" to ' // ENAME
                    CALL M3EXIT( 'GRDMOBIL', JDATE, JTIME, MESG, 2 )

                END IF          !  if write3() failed

198         CONTINUE            ! end loop on input variables

            CALL NEXTIME( JDATE, JTIME, TSTEP )

199     CONTINUE          !  end loop on time steps

999     CONTINUE          !  exit program
      
        CALL M3EXIT( 'GRDMOBIL', 0, 0, 
     &               'Normal completion  of PROGRAM GRDMOBIL', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )


        END

