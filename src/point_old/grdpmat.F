
        PROGRAM GRDPMAT

C***********************************************************************
C  program body starts at line 124
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
C       DSCGRID, LL2UTM, UTM2LL, PROMPTMFILE, GETYN, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  2/95 by CJC.
C       Revised   11/95 by CJC -- nonUTM coordinate systems
C       Revised    5/96 by M Houyoux -- more compact temporary storage
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
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        LOGICAL         DSCGRID
        LOGICAL         GETYN
        CHARACTER*16    PROMPTMFILE
        INTEGER         GETNUM
        INTEGER         STR2INT
        INTEGER         TRIMLEN
        
        EXTERNAL  DSCGRID, GETNUM, GETYN, PROMPTMFILE, STR2INT, TRIMLEN


C...........   LOCAL VARIABLES and their descriptions:
C...........   NOTE that ASC (area-source-category) ID's are 10-digit
C...........   unsigned integers which may be treated as a leading 7-digit
C...........   field, and a trailing 3-digit field.  *7 and *3 arrays below
C...........   follow this scheme with parallel arrays

        INTEGER     COL
        INTEGER     ROW
        INTEGER     CMAX, CMIN, CSUM
        REAL        XX, YY, XDUM, YDUM
        REAL        XX0, YY0
        REAL        XX1, YY1
        REAL        DDX, DDY
                                 
C...........   point Sources Table
        
        INTEGER     IFIPS( NPSRC )  !  source FIPS (county) ID
        REAL        XLOCA( NPSRC )
        REAL        YLOCA( NPSRC )
        
C...........   Scratch Gridding Matrix

        INTEGER     INDX( NPSRC )
        INTEGER     GN  ( NPSRC )
        INTEGER     SN  ( NPSRC )

C...........   Gridding Matrix

        INTEGER     NX( NGRID )
        INTEGER     IX( NPSRC )

        COMMON  / GRIDMAT / NX, IX

C...........   Other local variables
        
        INTEGER         C, R, I, J, K  !  loop counters.
        INTEGER         NBLANK         !  number of sources not in grid.  Will
                                       !  appear as blanks in SN( ) and GN( )

        CHARACTER*16    ENAME   !  logical name for emission source  input file
        CHARACTER*16    GNAME   !  logical name for grid matrix     output file
        INTEGER         LDEV    !  log-device
        LOGICAL         GFLAG   !  generate output gridding matrix if true
        LOGICAL         OFLAG   !  Did array overflow occur

        CHARACTER*16    CORDN   !  GRIDDESC coord-name -- PNTS  file
        CHARACTER*16    COORD   !  GRIDDESC coord-name -- PGMAT file
        CHARACTER*16    SCRBUF  !  Scratch buffer
        CHARACTER*256   MESG    !  message buffer for m3err(), m3exit()


C***********************************************************************
C   begin body of program GRDPMAT
        
        LDEV = INIT3()
        
        WRITE( *,92000 ) 
     &  ' ',
     &  'Program GRDPMAT to take the POINT SOURCE file produced by',
     &  'RAWPOINT and produce the POINT SOURCE GRIDDING-COEFFICIENT',
     &  'MATRIX file', 
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
            CALL M3EXIT( 'GRDPMAT', 0, 0, 'Ending Program.', 2 )
        END IF


C.......   Get file name; open actual-scc file
        
        ENAME = PROMPTMFILE( 
     &          'Enter logical name for POINT SOURCE input file',
     &          FSREAD3, 'PNTS', 'GRDPMAT' )
        
        IF( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'GRDPMAT', 0,0, 
     &                   'Could not read description for "' //
     &                   ENAME( 1:TRIMLEN( ENAME ) ) // '"', 2 )
        END IF
        CORDN = GDNAM3D

C.......   Initializations:
        
        DO  11  C = 1, NGRID
            NX( C ) = 0         !  for Gridding Matrix
11      CONTINUE

C.......   Read the grid parameters from GRIDDESC

        IF ( .NOT. DSCGRID( GRDNM, COORD, GDTYP3D, 
     &                 P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &                 XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &                 NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            SCRBUF = GRDNM
            MESG   = '"' // SCRBUF( 1:TRIMLEN( SCRBUF ) ) //
     &               '" not found in GRIDDESC file'
            CALL M3EXIT( 'GRDPMAT', 0, 0, MESG, 2 )


        ELSE IF ( CORDN .NE. COORD ) THEN

            MESG = 'Coord system mismatch--input "' // 
     &              CORDN( 1:TRIMLEN( CORDN ) ) //
     &              '", output "' // COORD( 1:TRIMLEN( COORD ) )
            CALL M3EXIT( 'GRDPMAT', 0, 0, MESG, 2 )

        ELSE IF ( NCOLS .NE. NCOLS3D ) THEN

            WRITE( MESG, 94050 ) 
     &              'Grid-cols mismatch:  dim=', NCOLS,
     &              'GRIDDESC=', NCOLS3D
            CALL M3EXIT( 'GRDPMAT', 0, 0, MESG, 2 )

        ELSE IF ( NROWS .NE. NROWS3D ) THEN

            WRITE( MESG, 94050 ) 
     &              'Grid-rows mismatch:  dim=', NROWS,
     &              'GRIDDESC=', NROWS3D
            CALL M3EXIT( 'GRDPMAT', 0, 0, MESG, 2 )

        END IF      !  if dscgrid() failed

C.......   Store grid parameters for later processing
        
        XX0   = SNGL( XORIG3D )
        YY0   = SNGL( YORIG3D )
        XX1   = XX0 + FLOAT( NCOLS3D ) * SNGL( XCELL3D )
        YY1   = YY0 + FLOAT( NROWS3D ) * SNGL( YCELL3D )
        DDX   = 1.0 / SNGL( XCELL3D )
        DDY   = 1.0 / SNGL( YCELL3D )


C.......   Get file name; open output gridding matrix file

        FTYPE3D = SMATRX3
        STIME3D = 0
        TSTEP3D = 0
        NVARS3D = 0
        NCOLS3D = NPSRC
        NROWS3D = NGRID
        NLAYS3D = 1
        NTHIK3D = NPSRC
        NLAYS3D = 1
        NTHIK3D = NPSRC
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3
        GDNAM3D = GRDNM         !  grid characteristics from DSCGRID() above

        FDESC3D( 1 ) = 'NC Point Source gridding-coefficient matrix'
        DO  33  K = 2, MXDESC3
            FDESC3D( K ) = ' '
33      CONTINUE

        GNAME = PROMPTMFILE( 
     &  'Enter logical name for GRIDDING MATRIX output file',
     &          FSUNKN3, 'PGMAT', 'GRDPMAT' )

        GFLAG = ( GNAME( 1:5 ) .NE. 'NONE ' )


C.......   Read in the area source emissions FIP and ASC tables:
        
        CALL M3MSG2( 'Reading in POINT SOURCES file...' )

        IF ( .NOT. READ3( ENAME, 'IFIP', ALLAYS3, 0, 0,  IFIPS ) ) THEN
            CALL M3EXIT( 'GRDPMAT', 0, 0, 
     &      'Error reading variable IFIP from POINT SOURCES file', 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'XLOCA', ALLAYS3, 0, 0,  XLOCA ) ) THEN
            CALL M3EXIT( 'GRDPMAT', 0, 0, 
     &      'Error reading variable XLOCA from POINT SOURCES file', 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'YLOCA', ALLAYS3, 0, 0,  YLOCA ) ) THEN
            CALL M3EXIT( 'GRDPMAT', 0, 0, 
     &      'Error reading variable YLOCA from POINT SOURCES file', 2 )
        END IF


C...........   Compute gridding matrix
        
        CALL M3MSG2( 'Computing gridding matrix...' )

        OFLAG = .FALSE.
        NBLANK = 0

        DO  111  I = 1, NPSRC

C.............  Initialize temporary storage table
            GN  ( I ) = 0
            SN  ( I ) = 0
            INDX( I ) = I

            XX = XLOCA( I )
            YY = YLOCA( I ) 
            
            XDUM = DDX * ( XX - XX0 )
            IF ( XDUM .LT. 0.0  )  THEN
                NBLANK = NBLANK + 1
                GO TO  111
            ENDIF

            COL = 1 + INT( XDUM )
            IF ( COL .GT. NCOLS )  THEN
                NBLANK = NBLANK + 1
                GO TO  111
            ENDIF

            YDUM = DDY * ( YY - YY0 )
            IF ( YDUM .LT. 0.0  )  THEN
                NBLANK = NBLANK + 1
                GO TO  111
            ENDIF

            ROW = 1 + INT( YDUM )
            IF ( ROW .GT. NROWS )  THEN
                NBLANK = NBLANK + 1
                GO TO  111
            ENDIF
            
            C = COL  +  NCOLS * ( ROW - 1 )

            NX( C ) = NX( C ) + 1

            GN  ( I ) = C
            SN  ( I ) = I

111     CONTINUE        !  end loop on sources I, computing gridding matrix.
        
        CMAX = NX( 1 )
        CMIN = CMAX
        CSUM = CMAX
        
        IF ( GFLAG ) THEN          ! write gridding matrix to file

C.......   Sort temporary gridding matrix
           CALL SORTI2( NPSRC, INDX, GN, SN )
            
C.......   Compress matrix into I/O representation from scratch representation.
C.......   Compute statistics.
            
            K = 0
                
            DO  144  R = 1, NGRID
            
                J = NX( R )
                   
                IF (      J .GT. CMAX ) THEN
                    CMAX = J
                ELSE IF ( J .LT. CMIN ) THEN
                    CMIN = J
                END IF

                IF( J .NE. 0 ) THEN
                   
                    DO 143  C = 1, J
                        K       = K + 1
                        IX( K ) = SN( INDX( NBLANK + K ) )
143                 CONTINUE

                ENDIF
                   
144         CONTINUE    !  end of loop on cells K for this FIP
            
            IF ( K .GT. NPSRC ) THEN
                WRITE( MESG,94050 ) 
     &                  'Number of coeffs', K, 
     &                  'exceeds NPSRC=', NPSRC
                CALL M3EXIT( 'GRDPMAT', 0, 0, MESG, 3 ) 
            END IF
                   
            CALL M3MSG2( 'Writing out GRIDDING MATRIX file...' )

            IF ( .NOT. WRITE3( GNAME, 'ALL', 0, 0, NX ) ) THEN
                CALL M3EXIT( 'GRDPMAT', 0, 0, 
     &              'Error writing GRIDDING MATRIX file.', 2 )
            END IF

        ELSE            !....  not GFLAG:  report matrix stats only
            
            K = 0
                
            DO  155  R = 2, NGRID
            
                J = NX( R )
                K = K + J
                   
                IF      ( J .GT. CMAX ) THEN
                    CMAX = J
                ELSE IF ( J .LT. CMIN ) THEN
                    CMIN = J
                END IF

155         CONTINUE    !  end of loop on cells K for this FIP

        END IF  !  if gflag, or not

                        
C.......   Report statistics:

        XDUM = FLOAT( K ) / FLOAT( NGRID )
        WRITE( LDEV,92010 ) 
     &      'Total number of coefficients   ', K   ,
     &      'Max  number of sources per cell', CMAX,
     &      'Min  number of sources per cell', CMIN
        WRITE( LDEV,92020 ) 
     &      'Mean number of sources per cell', XDUM,
     &      ' '
        WRITE( LDEV,92010 ) 
     &      'NOTE: Number of sources excluded from grid was', NBLANK,
     &      ' '
           

C...............   End of program

999     CONTINUE
      
        IF ( OFLAG .AND. GFLAG ) THEN
            CALL M3EXIT( 'GRDPMAT', 0, 0, 
     &               'ARRAY OVERFLOW in program GRDPMAT', 2 )
        ELSE 
            CALL M3EXIT( 'GRDPMAT', 0, 0, 
     &               'SUCCESSFUL COMPLETION of program GRDPMAT', 0 )
        END IF

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, :, I12 )

92020   FORMAT( 5X, A, :, F17.4 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A16 )


C...........   Internal buffering formats............ 94xxx

94050   FORMAT( 10( A, :, I10, :, 2X ) )


        END

