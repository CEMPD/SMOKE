
        PROGRAM MAPPOINT

C***********************************************************************
C  program body starts at line 115
C
C  DESCRIPTION:
C       Construct point source gridding matrix from data contained in
C       EPS-style surrogates file.
C
C  PRECONDITIONS REQUIRED:
C       Sorted input data for point sources; grid name for use with
C       DSCGRID().  Logical names already set.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       PROMPTFFILE, PROMPTMFILE
C
C  REVISION  HISTORY:
C       Prototype  1/97 by CJC.
C
C****************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File Version  @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C 
C See file COPYRIGHT for conditions of use.
C 
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C 
C  env_progs@mcnc.org
C 
C Pathname: $Source$
C Last updated: $Date$ 
C
C****************************************************************************/

      IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        LOGICAL         DSCGRID
        INTEGER         FIND1
        LOGICAL         GETYN
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN
        
        EXTERNAL  DSCGRID, FIND1, GETYN, PROMPTFFILE, PROMPTMFILE,
     &            TRIMLEN


C...........   LOCAL VARIABLES and their descriptions:
C...........   NOTE that ASC (area-source-category) ID's are 10-digit
C...........   unsigned integers which may be treated as a leading 7-digit
C...........   field, and a trailing 3-digit field.  *7 and *3 arrays below
C...........   follow this scheme with parallel arrays

        INTEGER     COL
        INTEGER     ROW
        REAL        XX, YY, XDUM, YDUM
        REAL        XX0, YY0
        REAL        XX1, YY1
        REAL        DDX, DDY
                                 
C...........   Point Sources Table
        
        INTEGER     INDEXA( NPSRC )
        INTEGER     IFIP  ( NPSRC )
        INTEGER     ISCC  ( NPSRC )
        INTEGER     IPLANT( NPSRC )
        REAL        XLOCA ( NPSRC )
        REAL        YLOCA ( NPSRC )

C...........   SCC list        

        INTEGER     NSCC
        INTEGER     SCCS( NPSCC )  !  scc list to be counted

C...........   Output counts

        REAL		COUNT( NCOLS, NROWS )

C...........   Other local variables
        
        INTEGER         C, R, I, J, K, L  !  loop counters.
        CHARACTER*16    ENAME   !  logical name for emission source  input file
        CHARACTER*16    MNAME   !  logical name for grid matrix     output file
        INTEGER         LDEV    !  log-device unit number
        INTEGER         SDEV    !  SCC-file unit number
        CHARACTER*16    CORDN   !  GRIDDESC coord-name -- PNTS  file
        CHARACTER*16    COORD   !  GRIDDESC coord-name -- GRIDDESC file
        CHARACTER*256   MESG    !  message buffer

C***********************************************************************
C   begin body of program MAPPOINT
        
        LDEV = INIT3()

        CALL INITEM( LDEV )
        
        WRITE( *,92000 ) 
     &  ' ',
     &  'Program MAPPOINT to take the POINT SOURCE file produced by',
     &  'RAWPOINT and ANY sorted list of SCC codes, and produces a',
     &  'gridded map of plant-counts for PLANTS (not stacks)',
     &  'found with SCC codes matching those in the list.',
     &  ' ',
     &  'You will need to enter the logical names for the input and',
     &  'output files and for the grid descriptions file GRIDDESC;',  
     &  'you need to have set them prior to program launch using ',
     &  ' ',
     &  '    setenv <logicalname> <pathname>',
     &  ' ',
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry.  Default responses are given in',
     &  'brackets [LIKE THIS] and can be accepted by hitting the',
     &  '<RETURN> key.',
     &  ' '
        
        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'MAPPOINT', 0, 0, 'Ending program', 2 )
        END IF


C.......   Get file name; open points file
        
        ENAME = PROMPTMFILE( 
     &          'Enter logical name for POINT SOURCE input file',
     &          FSREAD3, 'PNTS', 'MAPPOINT' )
        
        IF( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'MAPPOINT', 0,0, 
     &                   'Could not read description for "' //
     &                   ENAME( 1:TRIMLEN( ENAME ) ) // '"', 2 )
        END IF
        CORDN = GDNAM3D

        IF( NROWS3D .NE. NPSRC ) THEN
            CALL M3EXIT( 'MAPPOINT', 0,0,
     &                   'Bad source count in input PNTS file', 2 )
            
        ENDIF

        SDEV = PROMPTFFILE( 'Enter file with a list of SCC codes',
     &          .TRUE., .TRUE., 'PSCC', 'MAPPOINT' )

C.......   Initializations:
        
        IF ( .NOT. DSCGRID( GRDNM, COORD, GDTYP3D, 
     &                 P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &                 XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &                 NCOLS3D, NROWS3D, NTHIK3D ) ) THEN
            CALL M3EXIT( 'MAPPOINT', 0,0,
     &                   'DSCGRID() could not get grid description', 2 )
        END IF

        IF( NROWS .NE. NROWS3D .OR. 
     &      NCOLS .NE. NCOLS3D .OR.
     &      CORDN .NE. COORD        ) THEN
            CALL M3EXIT( 'MAPPOINT', 0,0,
     &                   'Bad grid parameters for ' // GRDNM, 2 )

        ENDIF

        DO 22  R = 1, NROWS
        DO 11  C = 1, NCOLS
            COUNT( C,R ) = 0.0
11      CONTINUE
22      CONTINUE


C.......   Store grid parameters for later processing
        
        XX0   = SNGL( XORIG3D )
        YY0   = SNGL( YORIG3D )
        XX1   = XX0 + FLOAT( NCOLS3D ) * SNGL( XCELL3D )
        YY1   = YY0 + FLOAT( NROWS3D ) * SNGL( YCELL3D )
        DDX   = 1.0 / SNGL( XCELL3D )
        DDY   = 1.0 / SNGL( YCELL3D )


C.......   Get file name; open output site count file

        FTYPE3D = GRDDED3
        STIME3D = 0
        TSTEP3D = 0
        NVARS3D = 1
        NLAYS3D = 1
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3
        VGTYP3D = IMISS3
        GDNAM3D = GRDNM         !  grid characteristics from DSCGRID() above

        VNAME3D( 1 ) = 'PLANT_COUNT'
        VTYPE3D( 1 ) = M3REAL       ! Must be real for PAVE's sake
        UNITS3D( 1 ) = 'NONE'
        VDESC3D( 1 ) = 
     &  'Number of plants per grid cell matching specified SCC list'

        FDESC3D( 1 ) = 'Count of plants matching SCC list'
        DO  33  K = 2, MXDESC3
            FDESC3D( K ) = ' '
33      CONTINUE

        MNAME = PROMPTMFILE( 
     &  'Enter logical name for SITE COUNT output file',
     &          FSUNKN3, 'PGSCC', 'MAPPOINT' )


C.......   Read in the area source emissions FIP and ASC tables:
        
        CALL M3MSG2( 'Reading in POINT SOURCES file...' )

        IF ( .NOT. READ3( ENAME, 'IFIP', ALLAYS3, 0, 0,  IFIP ) ) THEN
            CALL M3EXIT( 'MAPPOINT', 0, 0, 
     &      'Error reading variable IFIP from POINT SOURCES file', 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'ISCC', ALLAYS3, 0, 0,  ISCC ) ) THEN
            CALL M3EXIT( 'MAPPOINT', 0, 0, 
     &      'Error reading variable ISCC from POINT SOURCES file', 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'IPLANT', ALLAYS3,0,0,  IPLANT ) ) THEN
            CALL M3EXIT( 'MAPPOINT', 0, 0, 
     &      'Error reading variable IPLANT from POINT SOURCES file', 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'XLOCA', ALLAYS3, 0, 0,  XLOCA ) ) THEN
            CALL M3EXIT( 'MAPPOINT', 0, 0, 
     &      'Error reading variable XLOCA from POINT SOURCES file', 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'YLOCA', ALLAYS3, 0, 0,  YLOCA ) ) THEN
            CALL M3EXIT( 'MAPPOINT', 0, 0, 
     &      'Error reading variable YLOCA from POINT SOURCES file', 2 )
        END IF

C...........   Initialize sorting index
        DO 39 I = 1, NPSRC
            INDEXA( I ) = I
39      CONTINUE
       
C...........   Sort by plant before SCC, still need FIPS code
        CALL SORTI3( NPSRC, INDEXA, IFIP, IPLANT, ISCC ) 

C...........   Read SCC list

        I = 0

44	CONTINUE		!  head of read-scc-list loop

            READ( SDEV,*,END=45 ) K
                
            I = I + 1

            IF( I .LE. NPSCC ) THEN
                SCCS( I ) = K
            ENDIF

        GO TO 44

45	CONTINUE		!  end read-scc-list loop

        NSCC = I

        IF( NSCC .GT. NPSCC ) THEN
            WRITE( MESG,94010 )
     &             'SCC list overflow; actual', NSCC,
     &             'dimensioned NPSCC=', NPSCC
            CALL M3EXIT( 'MAPPOINT', 0, 0, MESG, 2 ) 

        ENDIF 
        
        L = -1
        DO  111  I = 1, NPSRC

            J = INDEXA( I )

            IF ( FIND1( ISCC( J ), NSCC, SCCS ) .LE. 0 ) GO TO 111

            XX = XLOCA( J )
            YY = YLOCA( J ) 
            
            XDUM = DDX * ( XX - XX0 )
            IF ( XDUM .LT. 0.0  )    GO TO  111
               
            COL = 1 + INT( XDUM )
            IF ( COL .GT. NCOLS )    GO TO  111
               
            YDUM = DDY * ( YY - YY0 )
            IF ( YDUM .LT. 0.0  )    GO TO  111
               
            ROW = 1 + INT( YDUM )
            IF ( ROW .GT. NROWS )    GO TO  111
            
            K = IPLANT( J )
            IF ( K .EQ. L )          GO TO  111

            L = K
            COUNT( COL,ROW ) = COUNT( COL,ROW ) + 1.0
               
111     CONTINUE        !  end loop on sources I, computing gridding matrix.
        
                   
        CALL M3MSG2( 'Writing out SOURCE COUNT file...' )
           
        IF ( .NOT. WRITE3( MNAME, 'ALL', 0, 0, COUNT ) ) THEN
            CALL M3EXIT( 'MAPPOINT', 0, 0, 
     &          'Error writing GRIDDING MATRIX file.', 2 )
        END IF
     

C...............   End of program

999     CONTINUE
      
        CALL M3EXIT( 'MAPPOINT', 0, 0, 
     &           'SUCCESSFUL COMPLETION of program MAPPOINT', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT ( 10 ( A, :, I10, :, 2X ) )

        END

