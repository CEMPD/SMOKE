
C Version "@(#)$Id$ $Source$ $Date$ 

        PROGRAM CSGMOBIL

C***********************************************************************
C  program body starts at line 104
C
C  DESCRIPTION:
C       Use gridding-transform, model-speciation, and control matrices 
C       from programs GRDMMAT, SPCMMAT, and CTLMMAT to grid, speciate,
C       and control source-level area source emissions files.
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

      IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
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
        
        REAL            EMISV( NMSRC, NVTYPE, NETYPE )  !  input emissions.
        REAL            EMISC( NMSRC )  !  controlled speciated emissions
        REAL            EMISG( NGRID )  !  gridded output emissions
        
C...........   Gridding Matrix, Control Matrix, Speciation Matrix
      
        INTEGER         NX( NGRID )
        INTEGER         IX( MXGCOF )
        REAL            CX( MXGCOF )
        
        COMMON  / GRIDMAT / NX, IX, CX
           
        REAL            CC( NMSRC, NVTYPE, NETYPE )
        REAL            SS( NMSRC, NVTYPE, NSMATX )

C...........   Other local variables
        
        INTEGER         JDATE, JTIME, TSTEP
        INTEGER         NSTEPS
        INTEGER         LDEV
        INTEGER         E, I, J, S, T, V
        REAL            CFAC
        
        LOGICAL         CONTROLFLAG, OUTPUTFLAG

        CHARACTER*16    ANAME   !  logical name for area-source       input file
        CHARACTER*16    GNAME   !  logical name for grid       matrix input file
        CHARACTER*16    CNAME   !  logical name for control    matrix input file
        CHARACTER*16    SNAME   !  logical name for speciation matrix input file
        CHARACTER*16    ENAME   !  logical name for emission output file

        CHARACTER*16    SCRBUF  !  scratch name area
        CHARACTER*256   MESG    !  scratch message area


C***********************************************************************
C   begin body of program CSGMOBIL
        
        LDEV = INIT3()
        
        WRITE( *,92000 ) 
     &  ' ',
     &  'Program CSGMOBIL to take a sorted source level area source ',
     &  'emissions file and a gridding-transform matrix computed by',
     &  'program GRDMMAT, a speciation matrix computed by SPCMMAT, ',
     &  'optionally a control matrix produced by CTLMMAT, and produce',
     &  'gridded emissions.  The output time step structure of the ',
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
            WRITE( *,92000 ) ' ', 'Ending program CSGMOBIL', ' '
            CALL EXIT( 0 )
        END IF


C.......   Get file names; open input gridding, speciation, 
C.......   optional control matrix files, area source file, 
C.......   optional output file

        GNAME = PROMPTMFILE( 
     &          'Enter logical name for GRID MATRIX file',
     &          FSREAD3, 'MGMAT', 'CSGMOBIL' )

        SNAME = PROMPTMFILE( 
     &          'Enter logical name for SPECIATION MATRIX file',
     &          FSREAD3, 'MSMAT', 'CSGMOBIL' )

        CNAME = PROMPTMFILE( 
     &          'Enter logical name for CONTROL MATRIX file ' //
     &          'or "NONE"',
     &          FSREAD3, 'MCMAT', 'CSGMOBIL' )

        CONTROLFLAG = ( CNAME .NE. NONE )

        ANAME = PROMPTMFILE( 
     &          'Enter logical name for MOBIL SOURCE file',
     &          FSREAD3, 'MTMP', 'CSGMOBIL' )

C.......   Build description of output file, and optionally open it:

        IF ( .NOT. DESC3( ANAME ) ) THEN
            CALL M3EXIT( 'CSGMOBIL', 0, 0, 
     &              'Could not get description of file ' // ANAME, 2 )
        END IF

        IF ( .NOT.DSCGRID( GRDNM, SCRBUF, GDTYP3D, 
     &                     P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D, 
     &                     XORIG3D, YORIG3D, XCELL3D, YCELL3D, 
     &                     NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            CALL M3EXIT( 'CSGMOBIL', 0, 0,
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
     &          FSUNKN3, 'MGTS', 'CSGMOBIL' )

        OUTPUTFLAG = ( ENAME .NE. NONE )


C.......   Read gridding, speciation, and (optional) control matrices:

        WRITE( *,92000 ) ' ', 'Reading GRIDDING matrix.', ' '
        IF ( .NOT. READ3( GNAME, 'ALL', 1, 0, 0, NX ) ) THEN

            CALL M3EXIT( 'CSGMOBIL', 0, 0,
     &           'Could not read gridding matrix from file "' //
     &           GNAME( 1 : TRIMLEN( GNAME ) ) // '".', 2 )
            
        END IF      !  if read3() failed for gridding matrix

        WRITE( *,92000 ) ' ','Reading SPECIATION matrix.',' '
        IF ( .NOT. READ3( SNAME, 'ALL', 1, 0, 0, SS ) ) THEN

            CALL M3EXIT( 'CSGMOBIL', 0, 0,
     &           'Could not read SPECIATION matrix from file "'//
     &           SNAME( 1 : TRIMLEN( SNAME ) ) // '".', 2 )
            
        END IF      !  if read3() failed for speciation matrix


C.......   If CONTROLFLAG, read and apply control matrix:

        IF( CONTROLFLAG ) THEN

            WRITE( *,92000 ) 
     &      ' ', 'Reading CONTROL matrix.', ' '

            IF ( .NOT. READ3( CNAME, 'ALL', 1, 0, 0, CC ) ) THEN
                MESG = 'Could not read control matrix from '//CNAME
                CALL M3EXIT( 'CSGMOBIL', 0,0, MESG, 2 )
            END IF      !  if read3() failed for control matrix

C...........   CO and NOX controls, by source and vehicle type:

            DO  62  V = 1, NVTYPE
            DO  61  S = 1, NMSRC

                SS( S,V,1 ) = SS( S,V,1 ) * CC( S,V,1 )		!  CO
                SS( S,V,2 ) = SS( S,V,2 ) * CC( S,V,2 )		!  NO
                SS( S,V,3 ) = SS( S,V,3 ) * CC( S,V,2 )		!  NO2

61          CONTINUE
62          CONTINUE


C...........   VOC controls, by source, vehicle type, and emissions mode:

            DO  66  E = 4, NETYPE
            DO  65  V = 1, NVTYPE
            DO  64  S = 1, NMSRC

                CFAC = CC( S,V,E )

                DO  63  I = 1, NMPOL - 3

                    J = ( E - 4 ) * ( NMPOL - 3 ) + I
                    SS( S,V,J ) = SS( S,V,J ) * CFAC

63              CONTINUE
64          CONTINUE
65          CONTINUE
66          CONTINUE

        END IF          !  if controlflag
        

C.......   Transform and write out area source emissions values:

        IF ( CONTROLFLAG .AND. OUTPUTFLAG ) THEN
            WRITE( *,92000 ) 
     &      ' ', 
     &      'Calculating and writing ' //
     &      'controlled speciated gridded emissions...'
        ELSE IF ( CONTROLFLAG ) THEN
            WRITE( *,92000 ) 
     &      ' ', 
     &      'Calculating controlled speciated gridded emissions...'
        ELSE IF ( OUTPUTFLAG ) THEN
            WRITE( *,92000 ) 
     &      ' ', 
     &      'Calculating and writing speciated gridded emissions...'
        ELSE 
            WRITE( *,92000 ) 
     &      ' ', 
     &      'Calculating speciated gridded emissions...'
        END IF

        DO  199  T = 1, NSTEPS

C...................   Read in entire emissions field

            IF( .NOT. READ3( ANAME, 'ALL', ALLAYS3, 
     &                 JDATE, JTIME, EMISV ) ) THEN

                MESG = 'Could not read emissions from ' // ANAME
                CALL M3EXIT( 'CSGMOBIL', JDATE, JTIME, MESG, 2 )

            END IF          !  if read3() succeeds or not
            
C...................   Compute output CO

            CALL MSMATVEC( NMSRC, NVTYPE, 
     &                     SS   ( 1,1,1 ),
     &                     EMISV( 1,1,1 ), EMISC )

            CALL SMATVEC( NMSRC, NGRID, MXGCOF, NX, IX, CX,
     &                    EMISC, EMISG )

            IF ( OUTPUTFLAG ) THEN

                IF ( .NOT. WRITE3( ENAME, 'CO',
     &                             JDATE, JTIME, 
     &                             EMISG ) ) THEN

                    MESG = 'Could not write CO to ' // ENAME
                    CALL M3EXIT( 'CSGMOBIL', JDATE, JTIME, MESG, 2 )

                END IF          !  if write3() failed

            END IF          !  if outputflag


C...................   Compute output NO

            CALL MSMATVEC( NMSRC, NVTYPE, 
     &                     SS   ( 1,1,2 ),
     &                     EMISV( 1,1,2 ), EMISC )

            CALL SMATVEC( NMSRC, NGRID, MXGCOF, NX, IX, CX,
     &                    EMISC, EMISG )

            IF ( OUTPUTFLAG ) THEN

                IF ( .NOT. WRITE3( ENAME, 'NO',
     &                             JDATE, JTIME, 
     &                             EMISG ) ) THEN

                    MESG = 'Could not write NO to ' // ENAME
                    CALL M3EXIT( 'CSGMOBIL', JDATE, JTIME, MESG, 2 )

                END IF          !  if write3() failed

            END IF          !  if outputflag


C...................   Compute output NO2

            CALL MSMATVEC( NMSRC, NVTYPE, 
     &                     SS   ( 1,1,3 ),
     &                     EMISV( 1,1,3 ), EMISC )

            CALL SMATVEC( NMSRC, NGRID, MXGCOF, NX, IX, CX,
     &                    EMISC, EMISG )

            IF ( OUTPUTFLAG ) THEN

                IF ( .NOT. WRITE3( ENAME, 'NO2',
     &                             JDATE, JTIME, 
     &                             EMISG ) ) THEN

                    MESG = 'Could not write CO to ' // ENAME
                    CALL M3EXIT( 'CSGMOBIL', JDATE, JTIME, MESG, 2 )

                END IF          !  if write3() failed

            END IF          !  if outputflag


C...................   Compute output VOC-derived species:

            DO  198  V = 4, NMPOL

                CALL MSMATVEC( NMSRC, NVTYPE * NEMODE, 
     &                         SS   ( 1,1,V ),
     &                         EMISV( 1,1,4 ), EMISC )

                CALL SMATVEC( NMSRC, NGRID, MXGCOF, NX, IX, CX,
     &                        EMISC, EMISG )

                IF ( OUTPUTFLAG ) THEN

                    IF ( .NOT. WRITE3( ENAME, EMNAM( V ),
     &                                 JDATE, JTIME, 
     &                                 EMISG ) ) THEN

                        MESG = 'Could not write CO to ' // ENAME
                        CALL M3EXIT( 'CSGMOBIL', JDATE, JTIME, MESG, 2 )

                    END IF          !  if write3() failed

                END IF          !  if outputflag

198         CONTINUE            ! end loop on input variables

            CALL NEXTIME( JDATE, JTIME, TSTEP )

199     CONTINUE          !  end loop on time steps

999     CONTINUE          !  exit program
      
        CALL M3EXIT( 'CSGMOBIL', 0, 0, 
     &               'Normal completion  of PROGRAM CSGMOBIL', 0 )


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

