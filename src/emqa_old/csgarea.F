
        PROGRAM CSGAREA

C***********************************************************************
C  program body starts at line 128
C
C  DESCRIPTION:
C       Use gridding-transform, model-speciation, and control matrices 
C       from programs GRDAMAT, CHMAMAT, and CTLAMAT to grid, speciate,
C       and control source-level area source emissions files.
C
C  PRECONDITIONS REQUIRED:  
C       M3IO source-level area source emissions in canonical order used
C       by GRDAMAT
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; GETYN, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  3/95 by CJC.
C       Extra dimension checks 12/96 by CJC
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

C...........   PARAMETER

        CHARACTER*16    NONE
        PARAMETER     ( NONE = 'NONE' )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Area Sources input and output arrays (one variable at a time)
        
        REAL            EMISV( NASRC )  !  input emissions.
        REAL            EMISC( NASRC )  !  controlled speciated emissions
        REAL            EMISG( NGRID )  !  gridded output emissions
        
C...........   Gridding Matrix
      
        INTEGER         NS( NGRID )
        INTEGER         IS( NMATX )
        REAL            CS( NMATX )
        
        COMMON  / GRIDMAT / NS, IS, CS
           
C...........   Control Matrix

        REAL            CC( NASRC, NIPOL )

C...........   Speciation Matrix, and description of how it maps input 
C...........   variables to output variables:  SS( S, SPCREF( V,W ) )
C...........   is used to produce species # OUTREF in the output file
C...........   from species V in the input file
C...........   (NOTE that names in the speciation-matrix file are of the 
C...........   form "<INV-SPC-NAME>to<MODEL-SPC-NAME>")

        REAL            SS( NASRC, NMPOL )
        INTEGER         NSPCOUT( NIPOL )        !  # of active rows in SPCREF
        INTEGER         SPCREF ( NIPOL, NMPOL )
        
        CHARACTER*16    INVNAME( NMPOL )	!  inv-spc part of name

C...........   Other local variables
        
        INTEGER         JDATE, JTIME, TSTEP
        INTEGER         NSTEPS
        INTEGER         LDEV
        INTEGER         I, J, K, F, S, T, V
        
        LOGICAL         CONTROLFLAG, OUTPUTFLAG

        CHARACTER*16    ANAME   !  logical name for area-source       input file
        CHARACTER*16    GNAME   !  logical name for grid       matrix input file
        CHARACTER*16    CNAME   !  logical name for control    matrix input file
        CHARACTER*16    SNAME   !  logical name for speciation matrix input file
        CHARACTER*16    ENAME   !  logical name for emission output file
            
        CHARACTER*16    SCRBUF  !  scratch area
        CHARACTER*256   MESG    !  message buffer

C***********************************************************************
C   begin body of program CSGAREA
        
        LDEV = INIT3()
        
        WRITE( *,92000 ) 
     &  ' ',
     &  'Program CSGAREA to take a sorted source level area source ',
     &  'emissions file and a gridding-transform matrix computed by',
     &  'program GRDAMAT, a speciation matrix computed by SPCAMAT, ',
     &  'optionally a control matrix produced by CTLAMAT, and produce ',
     &  'gridded emissions.  The output time step structure of the ',
     &  'gridded file will be the same as the time step structure of ',
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
            CALL M3EXIT( 'CSGAREA', 0, 0, 'Ending program.', 2 )
        END IF

C.......   Get file names; open input gridding, speciation, 
C.......   optional control matrix files, area source file, 
C.......   optional output file

        GNAME = PROMPTMFILE( 
     &          'Enter logical name for GRID MATRIX file',
     &          FSREAD3, 'AGMAT', 'CSGAREA' )

        SNAME = PROMPTMFILE( 
     &          'Enter logical name for SPECIATION MATRIX file',
     &          FSREAD3, 'ASMAT', 'CSGAREA' )

        CNAME = PROMPTMFILE( 
     &          'Enter logical name for CONTROL MATRIX file ' //
     &          'or "NONE"',
     &          FSREAD3, 'ACMAT', 'CSGAREA' )

        CONTROLFLAG = ( CNAME .NE. NONE )

        ANAME = PROMPTMFILE( 
     &          'Enter logical name for AREA SOURCE file',
     &          FSREAD3, 'ATMP', 'CSGAREA' )

C.......   Build description of output file, and optionally open it:

        IF ( .NOT. DESC3( ANAME ) ) THEN
            CALL M3EXIT( 'CSGAREA', 0, 0, 
     &              'Could not get description of file ' // ANAME, 2 )
        ELSE IF ( NROWS3D .NE. NASRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  EMISSIONS file:', NROWS3D,
     &      'program:', NASRC
            CALL M3EXIT( 'CSGAREA', 0, 0, MESG, 2 )
        END IF

        IF ( .NOT.DSCGRID( GRDNM, SCRBUF, GDTYP3D, 
     &                     P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D, 
     &                     XORIG3D, YORIG3D, XCELL3D, YCELL3D, 
     &                     NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            CALL M3EXIT( 'CSGAREA', 0, 0,
     &              'Could not get description of grid ' // GRDNM, 2 )

        END IF

        NVARS3D = NMPOL
        FTYPE3D = GRDDED3   !  shares most of file-description with input file.
        NLAYS3D = 1
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

        DO  33 V = 1, NMPOL
            VNAME3D( V ) = EMNAM( V )
            UNITS3D( V ) = 'moles/hour'
            VDESC3D( V ) = 'gridded speciated emissions'
            VTYPE3D( V ) = m3REAL
33      CONTINUE
        FDESC3D( 1 ) = 'area source speciated gridded emissions.'
        
        ENAME = PROMPTMFILE( 
     &          'Enter logical name for OUTPUT EMIS file ' //
     &          'or "NONE"',
     &          FSUNKN3, 'AGTS', 'CSGAREA' )

        OUTPUTFLAG = ( ENAME .NE. NONE )


C.......   Read gridding, speciation, and (optional) control matrices:

        CALL M3MSG2( 'Reading GRIDDING matrix...' )

        IF ( .NOT. READ3( GNAME, 'ALL', 1, 0, 0, NS ) ) THEN

            CALL M3EXIT( 'CSGAREA', 0, 0,
     &           'Could not read gridding matrix from file "' //
     &           GNAME( 1 : TRIMLEN( GNAME ) ) // '".', 2 )
            
        END IF      !  if read3() failed for gridding matrix

        CALL M3MSG2( 'Reading SPECIATION matrix...' )

        IF ( .NOT. READ3( SNAME, 'ALL', 1, 0, 0, SS ) ) THEN

            CALL M3EXIT( 'CSGAREA', 0, 0,
     &           'Could not read SPECIATION matrix from file "'//
     &           SNAME( 1 : TRIMLEN( SNAME ) ) // '".', 2 )
            
        END IF      !  if read3() failed for speciation matrix


C.......   Build description of how to use speciation matrix:

        IF ( .NOT. DESC3( SNAME ) ) THEN
            WRITE( MESG, '(2A)' ) 
     &      'Could not get description of file ', SNAME
            CALL M3EXIT( 'CSGAREA', 0, 0, MESG, 2 )
        ELSE IF ( NROWS3D .NE. NASRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  S-MATRIX file:', NROWS3D,
     &      'program:', NASRC
            CALL M3EXIT( 'CSGAREA', 0, 0, MESG, 2 )
        END IF

        DO  54  V = 1, NIPOL
            NSPCOUT( V ) = 0
54      CONTINUE

        DO  55  V = 1, NMPOL
                
            K = INDEX( VNAME3D( V ), '_to_' )
            IF( K .LT. 1 ) THEN
                CALL M3EXIT( 'CSGAREA', 0, 0,
     &              'Bad conversion factor name "' // 
     &              VNAME3D( V )( 1:TRIMLEN( VNAME3D( V ) ) ) //
     &              '" in ' // SNAME, 2 )
            END IF
            INVNAME( V ) = VNAME3D( V )( 1 : K-1 )

            I = INDEX1( INVNAME( V ), NIPOL, EINAM )
            IF ( I .LT. 0 ) THEN
                CALL M3EXIT( 'CSGAREA', 0, 0, 
     &              'Bad input species "' // 
     &              INVNAME( V )( 1:TRIMLEN( INVNAME( V ) ) ) // 
     &              '" in ' // SNAME, 2 )
            END IF

            SCRBUF = VNAME3D( V )( K+4:16 )
            J = INDEX1( SCRBUF, NMPOL, EMNAM )
            IF ( J .LT. 0 ) THEN
                CALL M3EXIT( 'CSGAREA', 0, 0, 
     &              'Bad output species "' // 
     &              SCRBUF( 1:TRIMLEN( SCRBUF ) ) //
     &              '" in ' // SNAME,  2 )
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
                CALL M3EXIT( 'CSGAREA', 0, 0, MESG, 2 )
            ELSE IF ( NROWS3D .NE. NASRC ) THEN
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  C-MATRIX file:', NROWS3D,
     &          'program:', NASRC
                CALL M3EXIT( 'CSGAREA', 0, 0, MESG, 2 )
            END IF

            DO  66 S = 1, NIPOL

                SCRBUF = 'CONTROL_' // EINAM( S )( 1:8 )
                IF ( .NOT. READ3( CNAME, SCRBUF, 1, 0, 0, CC( 1,S ) ) 
     &             ) THEN
                       
                    CALL M3EXIT( 'CSGAREA', 0, 0,
     &                   'Could not read "' // SCRBUF( 1:12 ) // 
     &                   '" from CONTROL FILE "' //
     &                   CNAME( 1 : TRIMLEN( CNAME ) ) // '".',
     &                   2 )
                       
                END IF      !  if read3() failed for control matrix
         
C...............   If this inventory species name matches the
C...............   leading portion of the speciation-matrix-vble name,
C...............   multiply this diagonal species-conversion matrix
C...............   by this inventory-species control matrix
C...............   (recall that INVNAME(V) = VNAME3D(V)(1:5) --
C...............   e.g., the "42101" in "42101toCO"

                DO  65 V = 1, NMPOL
                    IF ( EINAM( S ) .EQ. INVNAME( V ) ) THEN
                        DO  64  I = 1, NASRC
                            SS( I,V ) = CC( I,S ) * SS( I,V )
64                      CONTINUE
                    END IF
65              CONTINUE

66          CONTINUE

        END IF          !  if controlflag
        
C.......   Transform and write out area source emissions values:

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

        DO  199  T = 1, NSTEPS

C...............   If this is a new month, or new day, write message
            IF ( MOD( T,24 ) .EQ. 1 ) THEN
 
                MESG = 'Processing ' //
     &                 DAYS( WKDAY( JDATE ) ) // MMDDYY( JDATE )
                CALL M3MSG2( MESG( 1:TRIMLEN( MESG ) ) )
 
            END IF

C.............  Write to screen because WRITE3 only writes to LDEV
            WRITE( *, 93020 ) HHMMSS( JTIME )

            DO  198  V = 1, NIPOL

                IF( READ3( ANAME, EINAM( V ), ALLAYS3, 
     &                     JDATE, JTIME, EMISV ) ) THEN

                    DO  197  J = 1,  NSPCOUT( V )

                        K = SPCREF( V,J )

                        CALL DMATVEC( NASRC, SS( 1,K ), EMISV, EMISC )

                        CALL SMATVEC( NASRC, NGRID, NMATX, NS, IS, CS,
     &                                EMISC, EMISG )

                        IF ( OUTPUTFLAG ) THEN

                            IF ( .NOT. WRITE3( ENAME, EMNAM( K ),
     &                                         JDATE, JTIME, 
     &                                         EMISG ) ) THEN

                                CALL M3EXIT( 'CSGAREA', JDATE, JTIME,
     &                          'Could not write "' // 
     &                          EMNAM( K )( 1:TRIMLEN( EMNAM( K ))) //
     &                          '" to ' // ENAME, 2 )
      
                            END IF          !  if write3() failed

                        END IF          !  if outputflag

197                 CONTINUE    !  end loop on output vbles for this input vble

                ELSE

                    CALL M3EXIT( 'CSGAREA', JDATE, JTIME, 
     &                          'Could not read "' //  EINAM( V ) // 
     &                          '" from ' // ANAME, 2 )

                END IF          !  if read3() succeeds or not

198         CONTINUE            ! end loop on input variables

            CALL NEXTIME( JDATE, JTIME, TSTEP )

199     CONTINUE          !  end loop on time steps

999     CONTINUE          !  exit program
      
        CALL M3EXIT( 'CSGAREA', 0, 0, 
     &               'Normal completion  of PROGRAM CSGAREA', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Error and warning message formats..... 91xxx


C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )


C...........   Formatted file I/O formats............ 93xxx

93010   FORMAT( A16 )

93020   FORMAT( 8X, 'at time ', A8 )
 

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 2X ) )


C...........   Miscellaneous formats................. 95xxx

95000   FORMAT ( /5X , A , $ )          !  generic prompt format.


        END

