
        PROGRAM CSELEV

C***********************************************************************
C  program body starts at line 128
C
C  DESCRIPTION:
C       Use model-speciation and control matrices from programs 
C       SPCPMAT and CTLPMAT to speciateand control source-level
C       elevated point source emissions files.
C
C  PRECONDITIONS REQUIRED:  
C       M3IO elevated source-level point source emissions produced
C       by TMPPOINT
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; GETYN, PROMPTMFILE, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  10/95 by CJC.
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
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
        INCLUDE 'TMDIMS3.EXT'   !  temporal
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         GETYN
        CHARACTER*10    HHMMSS
        INTEGER         INDEX1
        CHARACTER*14    MMDDYY
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN
        INTEGER         WKDAY

        EXTERNAL        GETYN, HHMMSS, INDEX1, MMDDYY, PROMPTMFILE,
     &                  TRIMLEN, WKDAY


C...........   PARAMETER

        CHARACTER*16    NONE
        PARAMETER     ( NONE = 'NONE' )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Point Sources input and output arrays (one variable at a time)

        INTEGER       NELVIN             !  current source-count
        INTEGER       INDXIN( NPELV )    !  subscript into netCDF inventory file
        REAL          EMISIN( NPELV, NIPOL )  !  emissions values (tons/yr)

        COMMON / ELEVIN / NELVIN, INDXIN, EMISIN

        INTEGER       NELPT            !  current source-count
        INTEGER       INDXE( NPELV )    !  subscript into netCDF inventory file
        REAL          EMISE( NPELV, NMPOL )  !  emissions values (tons/yr)

        COMMON / ELEVREC / NELPT, INDXE, EMISE

C...........   Control Matrix

        REAL            CC( NPSRC, NIPOL )

C...........   Speciation Matrix, and description of how it maps input 
C...........   variables to output variables:  SS( S, SPCREF( V,W ) )
C...........   is used to produce species # OUTREF in the output file
C...........   from species V in the input file
C...........   (NOTE that names in the speciation-matrix file are of the 
C...........   form "<INV-SPC-NAME>to<MODEL-SPC-NAME>")

        REAL            SS( NPSRC, NMPOL )
        INTEGER         NSPCOUT( NIPOL )        !  # of active rows in SPCREF
        INTEGER         SPCREF ( NIPOL, NMPOL )

        CHARACTER*16    INVNAME( NMPOL )	!  inv-spc part of name

C...........   Other local variables

        INTEGER         JDATE, JTIME, TSTEP
        INTEGER         NSTEPS
        INTEGER         LDEV
        INTEGER         I, J, K, F, S, T, V

        LOGICAL         CONTROLFLAG, OUTPUTFLAG

        CHARACTER*16    PNAME   !  logical name for point-source      input file
        CHARACTER*16    CNAME   !  logical name for control    matrix input file
        CHARACTER*16    SNAME   !  logical name for speciation matrix input file
        CHARACTER*16    ENAME   !  logical name for emission output file
        CHARACTER*16    SCRBUF  !  scratch point
        CHARACTER*256   MESG    !  message buffer


C***********************************************************************
C   begin body of program CSELEV

        LDEV = INIT3()

        CALL INITEM( LDEV )

        WRITE( *,92000 ) 
     &' ',
     &'Program CSELEV to take an elevated point source emissions file',
     &'a speciation matrix computed by SPCPMAT, optionally a control',
     &'matrix produced by CTLPMAT, and produce speciated controlled',
     &'elevated emissions.  The output time step structure of the ',
     &'output file will be the same as the time step structure of the',
     &'input point source emissions file.',
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
            CALL M3EXIT( 'CSELEV', 0, 0, 'Ending Program.', 2 )
        END IF


C.......   Get file names; open input gridding, speciation, 
C.......   optional control matrix files, point source file, 
C.......   optional output file

        SNAME = PROMPTMFILE( 
     &          'Enter logical name for SPECIATION MATRIX file',
     &          FSREAD3, 'PSMAT', 'CSELEV' )

        CNAME = PROMPTMFILE( 
     &          'Enter logical name for CONTROL MATRIX file ' //
     &          'or "NONE"',
     &          FSREAD3, 'PCMAT', 'CSELEV' )

        CONTROLFLAG = ( CNAME .NE. NONE )

        PNAME = PROMPTMFILE( 
     &          'Enter logical name for ELEVATED SOURCE file',
     &          FSREAD3, 'ETMP', 'CSELEV' )

C.......   Build description of output file, and optionally open it.
C.......   Borrow most of the file description from the input file.

        IF ( .NOT. DESC3( PNAME ) ) THEN
            CALL M3EXIT( 'CSELEV', 0, 0, 
     &          'Could not get description of file ' // PNAME, 2 )
        END IF
                           
        NSTEPS  = MXREC3D
        JDATE   = SDATE3D
        JTIME   = STIME3D
        TSTEP   = TSTEP3D
        NVARS3D = NMPOL

        IF ( CONTROLFLAG ) THEN

            DO  33 V = 1, NMPOL
                VNAME3D( V ) = EMNAM( V )
                UNITS3D( V ) = 'moles/hour'
                VDESC3D( V ) = 'controlled speciated elevated emissions'
                VTYPE3D( V ) = m3REAL
33          CONTINUE
            FDESC3D( 1 ) = 'controlled speciated elevated emissions.'

        ELSE

            DO  34 V = 1, NMPOL
                VNAME3D( V ) = EMNAM( V )
                UNITS3D( V ) = 'moles/hour'
                VDESC3D( V ) = 'speciated elevated emissions'
                VTYPE3D( V ) = m3REAL
34          CONTINUE
            FDESC3D( 1 ) = 'speciated elevated emissions.'

        END IF


        MESG = 'Enter logical name for OUTPUT EMIS file, or "NONE"'

        ENAME = PROMPTMFILE( MESG, FSUNKN3, 'ETSC', 'CSELEV' )

        OUTPUTFLAG = ( ENAME .NE. NONE )


C.......   Read gridding, speciation, and (optional) control matrices:

        CALL M3MSG2( 'Reading SPECIATION matrix...' )

        IF ( .NOT. READ3( SNAME, 'ALL', 1, 0, 0, SS ) ) THEN

            CALL M3EXIT( 'CSELEV', 0, 0,
     &           'Could not read SPECIATION matrix from file "'//
     &           SNAME, 2 )

        END IF      !  if read3() failed for speciation matrix


C.......   Build description of how to use speciation matrix:

        IF ( .NOT. DESC3( SNAME ) ) THEN
            CALL M3EXIT( 'CSELEV', 0, 0, 
     &      'Could not get description of file ' // SNAME, 2 )
        ELSE IF ( NROWS3D .NE. NPSRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  S-MATRIX file:', NROWS3D,
     &      'program:', NPSRC
            CALL M3EXIT( 'CSELEV', 0, 0, MESG, 2 )
        END IF

        DO  54  V = 1, NIPOL
            NSPCOUT( V ) = 0
54      CONTINUE

        DO  55  V = 1, NMPOL

            K = INDEX( VNAME3D( V ), '_to_' )
            IF( K .LT. 1 ) THEN
                CALL M3EXIT( 'CSELEV', 0, 0,
     &              'Bad conversion factor name "' // 
     &              VNAME3D( V )( 1:TRIMLEN( VNAME3D( V ) ) ) //
     &              '" in ' // SNAME, 2 )
            END IF
            INVNAME( V ) = VNAME3D( V )( 1 : K-1 )

            I = INDEX1( INVNAME( V ), NIPOL, EINAM )
            IF ( I .LT. 0 ) THEN
                CALL M3EXIT( 'CSELEV', 0, 0, 
     &          'Bad input species "' //
     &          INVNAME( V )( 1:TRIMLEN( INVNAME( V ) ) ) //
     &          '" in ' // SNAME, 2 )
            END IF

            SCRBUF = VNAME3D( V )( K+4:16 )
            J = INDEX1( SCRBUF, NMPOL, EMNAM )
            IF ( J .LT. 0 ) THEN
                CALL M3EXIT( 'CSELEV', 0, 0, 
     &          'Bad output species "' // 
     &          SCRBUF( 1:TRIMLEN( SCRBUF ) ) // 
     &          '" in ' // SNAME, 2 )
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
                CALL M3EXIT( 'CSELEV', 0, 0, MESG, 2 )
            ELSE IF ( NROWS3D .NE. NPSRC ) THEN
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  C-MATRIX file:', NROWS3D,
     &          'program:', NPSRC
                CALL M3EXIT( 'CSELEV', 0, 0, MESG, 2 )
            END IF

            DO  66 S = 1, NIPOL

                SCRBUF = 'CONTROL_' // EINAM( S )( 1:8 )
                IF ( .NOT. READ3( CNAME, SCRBUF, 1, 0, 0, CC( 1,S ) ) 
     &             ) THEN

                    CALL M3EXIT( 'CSELEV', 0, 0,
     &                   'Could not read "' // SCRBUF( 1:12 ) // 
     &                   '" from CONTROL FILE "' // CNAME, 2 )

                END IF      !  if read3() failed for control matrix

                DO  65 V = 1, NVARS3D
                    IF ( EINAM( S ) .EQ. INVNAME( V ) ) THEN
                        DO  64  I = 1, NPSRC
                            SS( I,V ) = CC( I,S ) * SS( I,V )
64                      CONTINUE
                    END IF
65              CONTINUE

66          CONTINUE

        END IF          !  if controlflag


C.......   Transform and write out point source emissions values:

        IF ( CONTROLFLAG .AND. OUTPUTFLAG ) THEN
            CALL M3MSG2( 'Calculating and writing ' //
     &      'controlled speciated elevated emissions...' )
        ELSE IF ( CONTROLFLAG ) THEN
            CALL M3MSG2( 
     &      'Calculating controlled speciated elevated emissions...' )
        ELSE IF ( OUTPUTFLAG ) THEN
            CALL M3MSG2(
     &      'Calculating and writing speciated elevated emissions...' )
        ELSE 
            CALL M3MSG2( 'Calculating speciated elevated emissions...' )
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

            IF( READ3( PNAME, ALLVAR3, ALLAYS3, 
     &                 JDATE, JTIME, NELVIN ) ) THEN

                IF ( T .EQ. 1 ) THEN	!  setup output counter, indexes

                    NELPT = NELVIN

                    DO  191  S = 1, NELVIN
                        INDXE( S ) = INDXIN( S )
191                 CONTINUE

                END IF		!  end setup step (T = 1)

                DO  198  S = 1, NELVIN

                    I = INDXE( S )

                    DO  197  V = 1, NIPOL
                    DO  196  J = 1, NSPCOUT( V )

                        K = SPCREF( V,J )
                        EMISE( S,K ) = SS( I,K ) * EMISIN( S,V )

196                 CONTINUE    !  end loop on output vbles for this input vble
197                 CONTINUE    !  end loop on input variables

198             CONTINUE        !  end loop on sources

                IF ( OUTPUTFLAG ) THEN

                    IF ( .NOT. WRITE3( ENAME, ALLVAR3,
     &                                 JDATE, JTIME, 
     &                                 NELPT ) ) THEN

                        CALL M3EXIT( 'CSELEV', JDATE, JTIME,
     &                  'Could not write time step to ' //
     &                  ENAME, 2 )

                    END IF          !  if write3() failed

                END IF          !  if outputflag

            ELSE		!  else read3() failed

                CALL M3EXIT( 'CSELEV', JDATE, JTIME, 
     &                       'Could not read time step from ' // PNAME,
     &                       2 )

            END IF          !  if read3() succeeds or not

            CALL NEXTIME( JDATE, JTIME, TSTEP )

199     CONTINUE          !  end loop on time steps

999     CONTINUE          !  exit program

        CALL M3EXIT( 'CSELEV', 0, 0, 
     &               'Normal completion  of PROGRAM CSELEV', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

C...........   Formatted file I/O formats............ 93xxx
 
93020   FORMAT( 8X, 'at time ', A8 )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 3( A, :, I7, :, 1X ) )


        END

