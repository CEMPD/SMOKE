C ****************************************************************************
C *
C * Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C *                System
C * File: @(#)$Id$
C *
C * COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
C * All Rights Reserved
C *
C * See file COPYRIGHT for conditions of use.
C *
C * Environmental Programs Group
C * MCNC--North Carolina Supercomputing Center
C * P.O. Box 12889
C * Research Triangle Park, NC  27709-2889
C *
C * env_progs@mcnc.org
C *
C * Pathname: $Source$
C * Last updated: $Date$ 
C *
C ****************************************************************************/
 
        PROGRAM GETRECS

C***********************************************************************
C  program body starts at line 164
C
C  DESCRIPTION:
C       Searches for a specific source, for all sources in a 
C       specific cell, or for combinations of source keys.  It
C       creates a ascii file which lists all details about the
C       source including source number, grid cell, if found in
C       gridding matrix, temporalization factors, control factors,
C       inventory pollutant emissions, and model species emissions.
C
C  PRECONDITIONS REQUIRED:
C       Completed files for all point source processing stages
C       including PNTS, PTMP, PGMAT, PSMAT, PCMAT 
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       FIND1, GETNUM, INDEX1, PROMPTFFILE, PROMPTMFILE, GETYN, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  10/95 by M Houyoux
C***********************************************************************
 
        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER       FIND1
        INTEGER       GETNUM
        LOGICAL       GETYN
        INTEGER       INDEX1
        INTEGER       PROMPTFFILE
        CHARACTER*16  PROMPTMFILE
        INTEGER       TRIMLEN

        EXTERNAL      FIND1, GETNUM, GETYN, INDEX1, PROMPTFFILE, 
     &                PROMPTMFILE, TRIMLEN

C.........  Define PNTS variables
        INTEGER     IFIP  ( NPSRC )  !  source FIPS (county) ID
        INTEGER     ISCC  ( NPSRC )  !  source SCC
        INTEGER     ISIC  ( NPSRC )  !  source SIC
        INTEGER     IPLANT( NPSRC )  !  Plant ID
        INTEGER     ISTACK( NPSRC )  !  Stack ID
        INTEGER     TZONES( NPSRC )  !  time zones
        INTEGER     TPFLAG( NPSRC )  !  applicability of temporal profile types
        INTEGER     INVYR ( NPSRC )  !  Inventory year
        REAL        XLOCA ( NPSRC )  !  UTM X-location (m)
        REAL        YLOCA ( NPSRC )  !  UTM Y-location (m)
        REAL        STKHT ( NPSRC )  !  stack height   (m)
        REAL        STKDM ( NPSRC )  !  stack diameter (m)
        REAL        STKTK ( NPSRC )  !  exhaust temperature (deg K)
        REAL        STKVE ( NPSRC )  !  exhaust velocity    (m/s)
        REAL        RULPEN( NPSRC, NIPOL )  !  rule penetration   fraction
        REAL        RULEFF( NPSRC, NIPOL )  !  rule effectiveness fraction
        REAL        CTLEFF( NPSRC, NIPOL )  !  control efficiency fraction
        REAL        EMISV ( NPSRC, NIPOL )  !  emissions values (tons/yr) 

        COMMON / EMISREC / IFIP, ISIC, ISCC, IPLANT, ISTACK, TZONES,
     &                    TPFLAG, INVYR, XLOCA, YLOCA, STKHT, STKDM, 
     &                    STKTK, STKVE, CTLEFF, RULEFF, RULPEN, EMISV

C.........  Define PTMP variables
        REAL        EMIST ( NPSRC )  !  ptmp emissions values (tons/yr) 
        REAL        EV

C.........  Define SMAT variables

        REAL            SS( NPSRC, NMPOL )
        INTEGER         NSPCOUT( NIPOL )        !  # of active rows in SPCREF
        INTEGER         SPCREF ( NIPOL, NMPOL )
        CHARACTER*16    INVNAME( NMPOL )        !  inv-spc part of name


C.........  Define GMAT variables

        INTEGER         NS( NGRID )
        INTEGER         IS( NMATX )
 
        COMMON  / GRIDMAT / NS, IS

C.........  Define CMAT variables

        REAL            CC( NPSRC, NIPOL )

C.........  Local parameters
        INTEGER         NEOUT
        PARAMETER     ( NEOUT = 500 )

C.........  Local variables
        INTEGER         EINX ( NPSRC )   ! index for output sources
        INTEGER         ISLOC( NEOUT )   ! Cell specific IS values
        INTEGER         NOUT             ! actual number of output emis
        INTEGER         LDEV             ! IO/API initialization unit #
        INTEGER         ODEV             ! number for output file
        INTEGER         SPNT ( NMPOL )   ! pointer to speciation factors
        INTEGER         XNUM ( NEOUT )   ! X-cell for output sources
        INTEGER         YNUM ( NEOUT )   ! Y-cell for output sources

        REAL            EOUT( NEOUT, NMPOL, 24) ! output emissions
        REAL            ETMP( NEOUT, NIPOL, 24) ! temporal factors

        CHARACTER*7     GSTATE( NEOUT ) ! state of gridding matrix for output

C.........  Other variables
        INTEGER         I, ICNT, IJ, IWCH 
        INTEGER         J, JDATE, JTIME, K, NINCELL, NSTEPS, NUMBER
        INTEGER         CELLNUM, F, FIP, PLANT, RN, S, SCC, SCNT, SDATE,
     &                  SIC, STACK, STIME, T, TSTEP, V, XN, YN

        REAL            CFAC, DDX, DDY, DX, DY, SFAC
        REAL            X0, XMIN, XMAX, Y0, YMIN, YMAX

        LOGICAL         CELLTRUE, CONTROLFLAG, RECTRUE


        CHARACTER*16    PNAME   !  logical name for point-source      input file
        CHARACTER*16    GNAME   !  logical name for grid       matrix input file
        CHARACTER*16    CNAME   !  logical name for control    matrix input file
        CHARACTER*16    SNAME   !  logical name for speciation matrix input file
        CHARACTER*16    TNAME   !  logical name for temporal point source file

        CHARACTER*16    SCRBUF  !  scratch point
        CHARACTER*256   MESG    !  scratch message buffer 
 

C***********************************************************************
C   begin body of program GETRECS

        LDEV = INIT3()

        WRITE( *,92000 )
     &  ' ',
     &  'Program GETRECS searches for a specific source, for all',
     &  'sources in a specific cell, or for combinations of source',
     &  'keys.  It creates a ascii file which lists all details about',
     &  'the source including source number, grid cell, if found in',
     &  'gridding matrix, temporalization factors, control factors,',
     &  'inventory pollutant emissions, and model species emissions.',
     &  ' ',
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").',
     &  ' ',
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry. Default responses are given in',
     &  'brackets [LIKE THIS] and can be accepted by hitting the',
     &  '<RETURN> key.',
     &  ' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            WRITE( *,92000 ) ' ', 'Ending program GETRECS', ' '
            CALL EXIT( 0 )
        END IF


C.........  Prompt for PNTS

        PNAME = PROMPTMFILE(
     &          'Enter logical name for POINT SOURCE input file >> ',
     &          FSREAD3, 'PNTS', 'GETRECS' )

C.........  Prompt for PTMP

        TNAME = PROMPTMFILE(
     &          'Enter logical name for POINT SOURCE file',
     &          FSREAD3, 'PTMP', 'GETRECS' )

C.........  Get header information from PNTS

        IF ( .NOT. DESC3( TNAME ) ) THEN
            CALL M3EXIT( 'GETRECS', 0, 0,
     &                  'Could not get description of file "'
     &                  // TNAME( 1:TRIMLEN( TNAME ) ) // '"', 2 )
        ENDIF
        JDATE = SDATE3D
        JTIME = STIME3D
        TSTEP = TSTEP3D

C.........  Prompt for starting date

        JDATE  = GETNUM( 0, 9999999, JDATE,
     &                  'Enter starting date (YYYYDDD)' )
        SDATE  = JDATE

C.........  Prompt for starting time

        JTIME  = GETNUM( 0, 235959, JTIME,
     &                  'Enter starting time  (HHMMSS)' )
        STIME  = JTIME

C.........  Prompt for time number of time steps

        NSTEPS = 24
        WRITE( MESG, 94000 ) 'Enter number of time steps of ',
     &                       TSTEP,' in HHMMSS.'
        NSTEPS = GETNUM( 1, 24, NSTEPS, MESG )
        
C.........  Prompt for SMAT

        SNAME = PROMPTMFILE(
     &          'Enter logical name for SPECIATION MATRIX file',
     &          FSREAD3, 'PSMAT', 'GETRECS' )

C.........  Prompt for GMAT

        GNAME = PROMPTMFILE(
     &          'Enter logical name for GRID MATRIX file',
     &          FSREAD3, 'PGMAT', 'GETRECS' )

C.........  Prompt for CMAT

        CNAME = PROMPTMFILE(
     &          'Enter logical name for CONTROL MATRIX file ' //
     &          'or "NONE"',
     &          FSREAD3, 'PCMAT', 'GETRECS' )

        CONTROLFLAG = ( CNAME .NE. 'NONE' )

C.........  Prompt for Cell or Source or Record number
        CELLTRUE = .FALSE.
        RECTRUE  = .FALSE.
        IWCH = 1
        IWCH = GETNUM( 1, 3, IWCH,
     &     'Enter 1 for cell, 2 for source info, 3 for record no.' )
        IF( IWCH .EQ. 1 ) CELLTRUE = .TRUE. 
        IF( IWCH .EQ. 3 ) RECTRUE  = .TRUE.

C.........  If Cell, Prompt for x-dir and y-dir cell numbers

        IF( CELLTRUE ) THEN

            XN  = 1
            XN  = GETNUM( 1, NCOLS, XN,
     &                    'Enter X-dir cell number' )
            YN  = 1
            YN  = GETNUM( 1, NROWS, YN,
     &                      'Enter Y-dir cell number' )

C.........  If Source, Prompt for FIP, SIC, SCC, PLANT, STACK

        ELSEIF( RECTRUE ) THEN

            RN  = 1
            RN  = GETNUM( 1, NPSRC, RN,
     &                    'Enter record number' )

        ELSE

            FIP   = 0
            FIP   = GETNUM( 0, 99999, FIP, 'Enter FIP' )

            SIC   = 0
            SIC   = GETNUM( 0, 9999, SIC, 'Enter SIC' )
           
            SCC   = 0
            SCC   = GETNUM( 0, 99999999, SCC, 'Enter SCC' )

            PLANT = 0
            PLANT = GETNUM( 0, 9999999, PLANT, 'Enter PLANT' )
           
            STACK = 0
            STACK = GETNUM( 0, 9999999, STACK, 'Enter STACK' )

        ENDIF

C.........  Get header information from PNTS

        IF ( .NOT. DESC3( GNAME ) ) THEN
            CALL M3EXIT( 'GETRECS', 0, 0,
     &                  'Could not get description of file "'
     &                  // GNAME( 1:TRIMLEN( GNAME ) ) // '"', 2 )
        END IF

C.........  Read ALL variables in PNTS

        IF ( .NOT. READ3( PNAME, 'ALL', ALLAYS3, 0, 0, IFIP ) ) THEN
            CALL M3EXIT( 'GETRECS', 0, 0,
     &                   'Error reading from file' // PNAME, 2)
        END IF

C.........  Store parameters
        DX  = SNGL( XCELL3D )
        DY  = SNGL( YCELL3D )
        DDX = 1.0 / DX
        DDY = 1.0 / DY
        X0  = XORIG3D
        Y0  = YORIG3D

C.........  Calculate cell location if it cell method selected
        IF( CELLTRUE ) THEN

            XMIN = X0 + ( XN-1 ) * DX
            XMAX = XMIN + DX
            YMIN = Y0 + ( YN-1 ) * DY
            YMAX = YMIN + DY

        ENDIF
C.........  Create pointer based on cell or on source parameters

        ICNT = 0
        DO 11 I = 1, NPSRC

            IF( CELLTRUE ) THEN

                IF( XLOCA( I ) .GE. XMIN .AND.
     &              XLOCA( I ) .LE. XMAX .AND.
     &              YLOCA( I ) .GE. YMIN .AND.
     &              YLOCA( I ) .LE. YMAX       ) THEN

                    ICNT = ICNT + 1
                    EINX( ICNT ) = I

                ENDIF
            
            ELSEIF( RECTRUE ) THEN

                ICNT = 1
                EINX( ICNT ) = RN    

            ELSE
           
                IF( ( FIP   .EQ. IFIP  ( I ) .OR. FIP   .EQ. 0 ) .AND.
     &              ( SIC   .EQ. ISIC  ( I ) .OR. SIC   .EQ. 0 ) .AND.
     &              ( SCC   .EQ. ISCC  ( I ) .OR. SCC   .EQ. 0 ) .AND.
     &              ( PLANT .EQ. IPLANT( I ) .OR. PLANT .EQ. 0 ) .AND.
     &              ( STACK .EQ. ISTACK( I ) .OR. STACK .EQ. 0 ) ) THEN

                    ICNT = ICNT + 1
                    IF( ICNT .LE. NEOUT ) THEN
                        EINX( ICNT ) = I
                    ENDIF
 
                ENDIF
                
            ENDIF

   11   CONTINUE

        NOUT = ICNT
        IF( NOUT .GT. NEOUT ) THEN
            WRITE( MESG, 94000 ) 'Sources found = ',NOUT,
     &                           ' but maximum = ',NEOUT,
     &                           '. Reset NEOUT and try again.'
            CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )

        ENDIF

C.........  Read SMAT

        IF ( .NOT. READ3( SNAME, 'ALL', 1, 0, 0, SS ) ) THEN
 
            CALL M3ERR( 'GETRECS', 0, 0,
     &           'Could not read SPECIATION matrix from file "'//
     &           SNAME( 1 : TRIMLEN( SNAME ) ) // '".',
     &           .TRUE. )
 
        END IF      !  if read3() failed for speciation matrix

C.......   Build description of how to use speciation matrix:
 
        IF ( .NOT. DESC3( SNAME ) ) THEN
            CALL M3EXIT( 'GETRECS', 0, 0,
     &                   'Could not get description of file ' // SNAME,
     &                  .TRUE. )
        END IF
 
        DO  13  V = 1, NIPOL
            NSPCOUT( V ) = 0
13      CONTINUE
 
        DO  15  V = 1, NMPOL
 
            K = INDEX( VNAME3D( V ), '_to_' )
            IF( K .LT. 1 ) THEN
                CALL M3EXIT( 'GETRECS', 0, 0,
     &              'Bad conversion factor name "' //
     &              VNAME3D( V )( 1:TRIMLEN( VNAME3D( V ) ) ) //
     &              '" in ' // SNAME, 2 )
            END IF
            INVNAME( V ) = VNAME3D( V )( 1:K-1 )
 
            I = INDEX1( INVNAME( V ), NIPOL, EINAM )
            IF ( I .LT. 0 ) THEN
                CALL M3ERR( 'GETRECS', 0, 0,
     &                  'Bad input species "' //
     &                  INVNAME( V )( 1:TRIMLEN( INVNAME( V ) ) ) //
     &                  '" in ' // SNAME,
     &                  .TRUE. )
            END IF
 
            SCRBUF = VNAME3D( V )( K+4:16 )
            J = INDEX1( SCRBUF, NMPOL, EMNAM )
            IF ( J .LT. 0 ) THEN
                CALL M3ERR( 'GETRECS', 0, 0,
     &                      'Bad output species "' //
     &                      SCRBUF( 1:TRIMLEN( SCRBUF ) ) //
     &                      '" in ' // SNAME, 2 )
            END IF
 
            F = NSPCOUT( I ) + 1
            SPCREF ( I,F ) = J
            NSPCOUT( I )   = F
 
15      CONTINUE

C.........  Read GMAT

        IF ( .NOT. READ3( GNAME, 'ALL', 1, 0, 0, NS ) ) THEN
 
            CALL M3ERR( 'GETRECS', 0, 0,
     &           'Could not read gridding matrix from file "' //
     &           GNAME( 1 : TRIMLEN( GNAME ) ) // '".',
     &           .TRUE. )
 
        END IF      !  if read3() failed for gridding matrix

C......... Check to make sure that the number records in cell matches GMAT info

        IF( CELLTRUE ) THEN

            CELLNUM = ( YN-1 )*NCOLS + XN
            NINCELL= NS( CELLNUM )
            IF( NOUT .NE. NINCELL ) THEN

                WRITE( * , 92100 ) 
     &          'Warning: the number of counted records for',
     &          'cell (', XN, ',', YN,') is ', NOUT,
     &          ' but the number in the gridding matrix is ', NUMBER 
            ENDIF

C.............  Count up sources up to cell of interest
            SCNT = 1
            DO 17 S = 1, CELLNUM-1
                SCNT = SCNT + NS( S )
   17       CONTINUE
C.............  Create sorted source list for cell
            DO 18 S = 1, NS( CELLNUM )
                ISLOC( S ) = S
   18       CONTINUE

            CALL SORTI1( 
     &           NS( CELLNUM ), ISLOC, IS( SCNT ) ) 

            DO 19 S = 1, NS( CELLNUM )
                ISLOC( S )= IS( SCNT + ISLOC( S )-1 )
   19       CONTINUE

        ENDIF

C.........  Read CMAT

        IF( CONTROLFLAG ) THEN

            DO  22 S = 1, NIPOL
 
                SCRBUF = 'CONTROL_' // EINAM( S )( 1:8 )
                IF ( .NOT. READ3( CNAME, SCRBUF, 1, 0, 0, CC( 1,S ) )
     &             ) THEN
 
                    CALL M3ERR( 'GETRECS', 0, 0,
     &                   'Could not read "' // SCRBUF( 1:12 ) //
     &                   '" from CONTROL FILE "' //
     &                   CNAME( 1 : TRIMLEN( CNAME ) ) // '".',
     &                   .TRUE. )
 
                END IF      !  if read3() failed for control matrix

22          CONTINUE
 
        END IF          !  if controlflag

C.........  Loop over time

        DO 66 T = 1, NSTEPS

C.............  Loop over inventory species
            DO 55 V = 1, NIPOL    

C.................  Read PTMP for given hour
                IF( READ3( TNAME, EINAM( V ), ALLAYS3,
     &                     JDATE, JTIME, EMIST ) ) THEN

C.....................  Loop over list of sources

                    DO 44 I = 1, NOUT

                        IF( .NOT. CELLTRUE ) THEN

                            XNUM( I ) = INT(DDX * 
     &                                  ( XLOCA( EINX( I ) ) -X0 ) ) + 1
                            YNUM( I ) = INT(DDY * 
     &                                  ( YLOCA( EINX( I ) ) -Y0 ) ) + 1
                        ELSE

                            XNUM( I ) = XN
                            YNUM( I ) = YN

                        ENDIF

                        CELLNUM = ( YNUM( I )-1 )*NCOLS + XNUM( I )

                        IF( CELLNUM .LE. 0     .OR. 
     &                      CELLNUM .GT. NGRID      ) THEN

                            GSTATE( I ) = 'OUTSIDE'
                            GO TO 44
                        ENDIF

                        IF( .NOT. CELLTRUE ) THEN

C.............................  Count up sources up to cell of interest
                            SCNT = 1
                            DO 29 S = 1, CELLNUM - 1
                                SCNT = SCNT + NS( S )
   29                       CONTINUE

C.............................  Create sorted source list for cell
                            DO 30 S = 1, NS( CELLNUM )
                                ISLOC( S ) = S
   30                       CONTINUE

                            CALL SORTI1( 
     &                           NS( CELLNUM ), ISLOC, IS( SCNT ) ) 

                            DO 31 S = 1, NS( CELLNUM )
                                ISLOC( S )= IS( SCNT + ISLOC( S )-1 )
   31                       CONTINUE

                        ENDIF

                        F = FIND1( EINX( I ),NS( CELLNUM ), ISLOC )
                                 
                        IF( F .GT. 0 ) THEN
                            GSTATE( I ) = 'FOUND'     
                        ELSE
                            GSTATE( I ) = 'UNFOUND'
                        ENDIF

C.........................  Loop over model species per inventory species

                        DO 33 J = 1,  NSPCOUT( V )

                            K = SPCREF( V,J )

C.............................  Multiply PTMP, SMAT, and CMAT for source

                            SFAC = SS( EINX( I ), K ) 
                            IF( CONTROLFLAG ) THEN
                                CFAC = CC( EINX( I ), V )
                            ELSE
                                CFAC = 1.0
                            ENDIF
                             
                            EOUT( I, K, T ) = 
     &                            EMIST( EINX( I ) ) * SFAC * CFAC

                            EV =  EMISV( EINX( I ), V )

                            IF( EV .GT. 0. ) THEN
                                ETMP( I, V, T ) = EMIST( EINX( I ) )/ EV
                            ELSE
                                ETMP( I, V, T ) = 0
                            ENDIF

   33                   CONTINUE
   44               CONTINUE

                ELSE ! Read failed

                    CALL M3EXIT( 'GETRECS', JDATE, JTIME,
     &                           'Could not read "' //
     &                           EINAM( V )( 1:TRIMLEN( EINAM( K ))) //
     &                           '" from ' // TNAME, 2 )
 
                ENDIF
   55       CONTINUE

            CALL NEXTIME( JDATE, JTIME, TSTEP )

   66   CONTINUE

C.........  Prompt for output temporal x-ref file name

        ODEV = PROMPTFFILE(
     &  'Enter name for ASCII OUTPUT file',
     &           .FALSE., .TRUE., 'LISTOUT', 'GETRECS' )
 
        IF ( ODEV .LE. 0 ) THEN
            CALL M3EXIT( 'GETRECS', 0, 0,
     &                   'Could not open ASCII OUTPUT file', 2 )
        END IF

C.........  Loop over list of sources

        DO 88 I = 1, NOUT
            IJ = EINX( I )

C.............  Write out all time INdependant data
            WRITE( ODEV, 93000 ) 
     &      I, IJ, XNUM( I ), YNUM( I ), GSTATE( I ), XLOCA( IJ ), 
     &      YLOCA( IJ ), IFIP( IJ ), 
     &      ISIC( IJ ), ISCC( IJ ), IPLANT( IJ ), ISTACK( IJ ), 
     &      TZONES( IJ ), TPFLAG( IJ ), INVYR( IJ ), STKVE( IJ ), 
     &      STKDM( IJ ), STKTK( IJ ), STKHT( IJ )

C.............  Write out time INdependant inventory emissions
            WRITE( ODEV, 93005 ) ( EINAM( V ), V = 1, NIPOL )
            WRITE( ODEV, 93010 ) 'INV', 
     &                          ( EMISV( EINX( I ), V ), V = 1, NIPOL )
   
C.............  Write out control coefficients

            IF( CONTROLFLAG ) 
     &      WRITE( ODEV, 93010 ) 'CTL', 
     &                          ( CC( IJ, V ), V = 1, NIPOL )

C.............  Create pointer for model species
            ICNT = 0
            DO 70 V = 1, NIPOL
                 DO 69 J = 1,  NSPCOUT( V )

                     K = SPCREF( V,J )
                     ICNT = ICNT + 1
                     SPNT( ICNT ) = K

   69            CONTINUE
   70       CONTINUE

C.............  Write out speciation coefficients

            WRITE( ODEV, 93005 ) ( EMNAM( SPNT( J ) ), J = 1, NMPOL )
            WRITE( ODEV, 93010 ) 'SPC',
     &                       ( SS( IJ, SPNT( J ) ), J=1,NMPOL )

C.............  Write out all time dependant data
            JDATE = SDATE
            JTIME = STIME
            DO 77 T = 1, NSTEPS

                WRITE( ODEV, 93020 ) 
     &                 JDATE, JTIME, ( ETMP( I, V, T ), V=1,NIPOL )

                WRITE( ODEV, 93030 )
     &                 ( EOUT( I, SPNT( J ), T ), J=1,NMPOL )

                CALL NEXTIME( JDATE, JTIME, TSTEP )

   77       CONTINUE
   88   CONTINUE 

        CALL M3EXIT( 'GETRECS', 0, 0,
     &               'Normal completion  of PROGRAM GETRECS', 0 )

C******************  FORMAT  STATEMENTS   ******************************
 
C...........   Informational (LOG) message formats... 92xxx
 
92000   FORMAT( 5X, A )
92100   FORMAT( 10x A, / 10X, A, I3, A, I3, A, I5 / 10X, A, I5)
 
C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( I3, 1X, I4, ' (',I3,',',I3,') ', A7, 2( 1X, E14.7), 1X,
     &          I5.5, 1X, I5, 1X, I8.8, 
     &          1X, I8, 1X, I8, 3( 1X, I4 ), 1X, 'V:', F6.2, 1X,
     &          'D:', F6.2, 1X, 'T:', F5.1, 1X, 'H:', F6.2 )

93005   FORMAT( 9X, 30( 4X, A5, 2X ) )
93010   FORMAT( 5X, A3 ,': ', 30( 1X, E9.3 ) )

93020   FORMAT( 5X, I7.7, 1X, I6.6, 30( 1X, E9.3 ) )

93030   FORMAT( 19X, 30( 1X, E9.3 ) )

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( A, I5, A, I5, A, I5, A, I5)


        END
