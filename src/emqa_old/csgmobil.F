
        PROGRAM CSGMOBIL

C***********************************************************************
C  program body starts at line 134
C
C  DESCRIPTION:
C       Use gridding-transform, model-speciation, and control matrices 
C       from programs GRDMMAT, SPCMMAT, and CTLMMAT to grid, speciate,
C       and control source-level mobile source emissions files.
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
C       Modified   6/96 by M Houyoux
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
        
        INCLUDE 'CHDGEN3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'MBDIMS3.EXT'   !  area-source dimensioning parameters
        INCLUDE 'TMDIMS3.EXT'   !  temporal dimensioning parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        INTEGER         CSVNAME
        LOGICAL         DSCGRID
        LOGICAL         GETYN
        CHARACTER*10    HHMMSS
        INTEGER         INDEX1
        CHARACTER*14    MMDDYY
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN
        INTEGER         WKDAY           !  day of week (1...7)
        
        EXTERNAL        CSVNAME, DSCGRID, GETYN, HHMMSS, INDEX1, MMDDYY,
     &                  PROMPTFFILE, PROMPTMFILE, TRIMLEN, WKDAY

C...........   PARAMETER

        CHARACTER*16    NONE
        PARAMETER     ( NONE = 'NONE' )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Area Sources input and output arrays (one variable at a time)
        
        REAL            EMISV( NMSRC, NVTYPE, MXEMIS )  !  input emissions.
        REAL            EMISC( NMSRC )  !  controlled speciated emissions
        REAL            EMISG( NGRID )  !  gridded output emissions
        
C...........   Gridding Matrix, Control Matrix, Speciation Matrix
      
        INTEGER         NX( NGRID )
        INTEGER         IX( NMATX )
        REAL            CX( NMATX )
        
        COMMON  / GRIDMAT / NX, IX, CX
           
c       REAL            CC( NMSRC, NVTYPE, MXEMIS )
        REAL            SS( NMSRC, NVTYPE, MXMPOL*MXEPRC )

C..........    Naming and matrix/vector matching arrays
        INTEGER         NPPM  ( MXMPOL )  ! Number of spcs facs per model spcs

        CHARACTER*16    EMTNAM( MXEMIS )  ! Names of emission types
        CHARACTER*16    LNAME             ! Var name from previous loop iteration
        CHARACTER*16    NMATCH            ! Name for matching
        CHARACTER*16    SPCFNM( MXMPOL*MXEPRC )  ! Names of speciation factors

C...........   Other local variables
        
        INTEGER         ICNT               ! Counter
        INTEGER         LDATE, JDATE, JTIME, TSTEP
        INTEGER         NETYPE             ! No. of scenario-spec MV emis types
        INTEGER         NMPOL              ! No. of scenario-spec model species
        INTEGER         NSNAME             ! No. of speciation factor names
        INTEGER         NSTEPS
        INTEGER         LDEV
        INTEGER         SIDX               ! Speciation matrix index no.
        INTEGER         TDEV               ! speciation list file unit
        INTEGER         I, J, K, L, N, M, S, T, V
        
        LOGICAL         CONTROLFLAG, OUTPUTFLAG

        CHARACTER*16    ANAME   !  logical name for area-source       input file
        CHARACTER*16    GNAME   !  logical name for grid       matrix input file
c       CHARACTER*16    CNAME   !  logical name for control    matrix input file
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
            CALL M3EXIT( 'CSGMOBIL', 0, 0, 'Ending program', 2 )
        END IF


C.......   Get file names; open input gridding, speciation, 
C.......   optional control matrix files, area source file, 
C.......   optional output file

        GNAME = PROMPTMFILE( 
     &          'Enter logical name for GRID MATRIX file',
     &          FSREAD3, 'MGMAT', 'CSGMOBIL' )

C.......   Get description of GNAME for dimensioning check
        IF ( .NOT. DESC3( GNAME ) ) THEN
            CALL M3EXIT( 'CSGMOBIL', 0, 0, 
     &              'Could not get description of file ' // GNAME, 2 )
        END IF

        IF( NTHIK3D .NE. NMSRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  MGMAT file:', NTHIK3D,
     &      'program:', NMSRC
            CALL M3EXIT( 'CSGMOBIL', 0, 0, MESG, 2 )
        END IF

        SNAME = PROMPTMFILE( 
     &          'Enter logical name for SPECIATION MATRIX file',
     &          FSREAD3, 'MSMAT', 'CSGMOBIL' )

C.......   Get description of SNAME for dimensioning check
        IF ( .NOT. DESC3( SNAME ) ) THEN
            CALL M3EXIT( 'CSGMOBIL', 0, 0, 
     &              'Could not get description of file ' // SNAME, 2 )
        END IF

        IF( NROWS3D .NE. NMSRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  MSMAT file:', NROWS3D,
     &      'program:', NMSRC
            CALL M3EXIT( 'CSGMOBIL', 0, 0, MESG, 2 )
        END IF

c       CNAME = PROMPTMFILE( 
c    &          'Enter logical name for CONTROL MATRIX file ' //
c    &          'or "NONE"',
c    &          FSREAD3, 'MCMAT', 'CSGMOBIL' )

c       CONTROLFLAG = ( CNAME .NE. NONE )
        CONTROLFLAG = .FALSE.

        TDEV = PROMPTFFILE(
     &           'Enter logical name for SPECIATION LIST file',
     &           .TRUE., .TRUE., 'SPCS', 'CSGMOBIL' )

        ANAME = PROMPTMFILE( 
     &          'Enter logical name for TIME-STEPPED MOBILE EMIS file',
     &          FSREAD3, 'MTMP', 'CSGMOBIL' )


C.......   Get description of ANAME for building output file and checking header
        IF ( .NOT. DESC3( ANAME ) ) THEN
            CALL M3EXIT( 'CSGMOBIL', 0, 0, 
     &              'Could not get description of file ' // ANAME, 2 )
        END IF

        IF( NROWS3D .NE. NMSRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  MTMP file:', NROWS3D,
     &      'program:', NMSRC
            CALL M3EXIT( 'CSGMOBIL', 0, 0, MESG, 2 )
        END IF

C.......   Get description of grid for building output file
        IF ( .NOT.DSCGRID( GRDNM, SCRBUF, GDTYP3D, 
     &                     P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D, 
     &                     XORIG3D, YORIG3D, XCELL3D, YCELL3D, 
     &                     NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            CALL M3EXIT( 'CSGMOBIL', 0, 0,
     &              'Could not get description of grid ' // GRDNM, 2 )

        END IF

C.......   Get speciation naming information from list file TDEV
C..........  Get number of scenario-specific MV emission types
C..........     ( Use VNAME3D as dummy so don't have to declare new array)
 
        NETYPE  = CSVNAME( TDEV, .TRUE., .FALSE., .TRUE., EMTNAM )
        NSNAME  = CSVNAME( TDEV, .TRUE., .TRUE. , .TRUE., SPCFNM )

C..........  Get number and names of scenario-specific model species

        NMPOL   = CSVNAME( TDEV, .FALSE., .TRUE. , .FALSE., VNAME3D ) 

C.........  Initialize "last name" holder "LNAME"
        L     = TRIMLEN( SPCFNM( 1 ) )
        J     = INDEX  ( SPCFNM( 1 )(   1:L ), '_' )
        K     = INDEX  ( SPCFNM( 1 )( J+1:L ), '_' ) + J
        LNAME = SPCFNM( 1 )( K+1:L )

        V     = 0
        ICNT  = 0
C.........  Count tail part of speciation factor variable names
        DO 22 I = 1, NSNAME

C.............  Find second "_" in variable names
            L   = TRIMLEN( SPCFNM( I ) )
            J   = INDEX  ( SPCFNM( I )(   1:L ), '_' )
            K   = INDEX  ( SPCFNM( I )( J+1:L ), '_' ) + J

C.............  If current name and last name don't match or on last 
            IF ( SPCFNM( I )( K+1:L ) .NE. LNAME ) THEN

C.................  Store count of names
                V = V + 1
                NPPM( V ) = ICNT

C.................  Reset last name and counter
                LNAME = SPCFNM( I )( K+1:L )
                ICNT  = 1

            ELSE
                ICNT = ICNT + 1

            ENDIF
            
   22   CONTINUE  ! End counting model species in speciation factor names

C.........  Store last counter value
        V = V + 1
        NPPM( V ) = ICNT

C.......   Finish the description of output file name
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
            UNITS3D( V ) = 'moles/hour'
            VDESC3D( V ) = 'Gridded speciated emissions'
            VTYPE3D( V ) = m3REAL
  33    CONTINUE
        FDESC3D( 1 ) = 'Mobile source speciated gridded emissions.'
                 
C.......   Optionally open an output file        

        ENAME = PROMPTMFILE( 
     &          'Enter logical name for OUTPUT EMIS file ' //
     &          'or "NONE"',
     &          FSUNKN3, 'MGTS', 'CSGMOBIL' )

        OUTPUTFLAG = ( ENAME .NE. NONE )

C.......   Read gridding, speciation, and (optional) control matrices:

        CALL M3MSG2( 'Reading GRIDDING matrix...' )

        IF ( .NOT. READ3( GNAME, ALLVAR3, 1, 0, 0, NX ) ) THEN

            CALL M3EXIT( 'CSGMOBIL', 0, 0,
     &           'Could not read gridding matrix from file "' //
     &           GNAME( 1 : TRIMLEN( GNAME ) ) // '".', 2 )
            
        END IF      !  if read3() failed for gridding matrix

        CALL M3MSG2( 'Reading SPECIATION matrix...' )

        IF ( .NOT. READ3( SNAME, ALLVAR3, ALLAYS3, 0, 0, SS ) ) THEN

            CALL M3EXIT( 'CSGMOBIL', 0, 0,
     &           'Could not read SPECIATION matrix from file "'//
     &           SNAME( 1 : TRIMLEN( SNAME ) ) // '".', 2 )
            
        END IF      !  if read3() failed for speciation matrix

C.......   If CONTROLFLAG, read and apply control matrix:
C.......   NOTE: This section removed until we decide in what way
C.......         a control matrix makes sense for MV processing

C.......   Transform and write out mobile source emissions values:

        IF ( CONTROLFLAG .AND. OUTPUTFLAG ) THEN
            CALL M3MSG2( 'Calculating and writing ' //
     &                   'controlled speciated gridded emissions...' )
        ELSE IF ( CONTROLFLAG ) THEN
            CALL M3MSG2( 
     &       'Calculating controlled speciated gridded emissions...' )
        ELSE IF ( OUTPUTFLAG ) THEN
            CALL M3MSG2(
     &       'Calculating and writing speciated gridded emissions...' )
        ELSE 
            CALL M3MSG2( 'Calculating speciated gridded emissions...' )

        END IF

C.........  Calculate emissions for all time steps available from MTMP
        LDATE = 0
        DO  199  T = 1, NSTEPS

C...............   If this is a new month, or new day, write message
            IF ( LDATE .NE. JDATE ) THEN
                    
                MESG = 'Processing ' //
     &                 DAYS( WKDAY( JDATE ) ) // MMDDYY( JDATE )
                CALL M3MSG2( MESG( 1:TRIMLEN( MESG ) ) )

            END IF

C.............  Write to screen because WRITE3 only writes to LDEV
            WRITE( *, 93020 ) HHMMSS( JTIME )

C.............   Read in entire emissions field
            IF( .NOT. READ3( ANAME, ALLVAR3, ALLAYS3, 
     &                 JDATE, JTIME, EMISV ) ) THEN

                MESG = 'Could not read emissions from ' // ANAME
                CALL M3EXIT( 'CSGMOBIL', JDATE, JTIME, MESG, 2 )

            END IF          !  if read3() succeeds or not

C.............  Loop through all output model species
            SIDX = 0
            DO  189 V = 1, NMPOL

C.................  Initialize output emissions array to 0.
                DO 177 S = 1, NMSRC
                    EMISC( S ) = 0.
  177           CONTINUE

C.................  Loop through processes associated w/ model species
                DO 179 J = 1, NPPM( V )

C.....................  Set index of SS for model species and process
                    SIDX = SIDX + 1

C.....................  Get index of EMISV for model species and process
C.........................  First find the location of the second "_"
                    L = TRIMLEN( SPCFNM ( SIDX ) ) 
                    N = INDEX  ( SPCFNM ( SIDX )(   1:L ), '_' )
                    M = INDEX  ( SPCFNM ( SIDX )( N+1:L ), '_' ) + N

C.........................  Then match the first 2 name sections to EMTNAM
                    NMATCH = SPCFNM( SIDX )( 1:M-1 )
                    K      = INDEX1( NMATCH, NETYPE, EMTNAM )

C.....................  Apply SS to EMISV and add result to current EMISC

                    CALL MSMATVEC( NMSRC, NVTYPE, 
     &                             SS   ( 1,1,SIDX ),
     &                             EMISV( 1,1,K    ), EMISC )

  179           CONTINUE

C.................  Apply gridding matrix
                CALL SMATVEC( NMSRC, NGRID, NMATX, NX, IX, CX,
     &                        EMISC, EMISG )

C.................  Write out model species emissions

                IF ( OUTPUTFLAG ) THEN

                    IF ( .NOT. WRITE3( ENAME, VNAME3D( V ),
     &                                 JDATE, JTIME, 
     &                                 EMISG ) ) THEN

                        MESG = 'Could not write ' // VNAME3D( V ) //
     &                         ' to ' // ENAME
                        CALL M3EXIT( 'CSGMOBIL', JDATE, JTIME, MESG, 2 )

                    END IF          !  if write3() failed

                END IF          !  if outputflag

  189       CONTINUE  !  End of loop on output model species

            LDATE = JDATE
            CALL NEXTIME( JDATE, JTIME, TSTEP )

  199   CONTINUE          !  end loop on time steps

  999   CONTINUE          !  exit program
      
        CALL M3EXIT( 'CSGMOBIL', 0, 0, 
     &               'Normal completion  of PROGRAM CSGMOBIL', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Error and warning message formats..... 91xxx


C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )


C...........   Formatted file I/O formats............ 93xxx

93010   FORMAT( A16 )

93020   FORMAT( 8X, 'at time ', A8 )
 

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 2X ) )
 

        END

