
        PROGRAM PCAREA

C***********************************************************************
C  program body starts at line 
C
C  DESCRIPTION:
C       Updates the SMOKE area source inventory file with a growth or
C       control matrix and outputs a new file.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
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
C****************************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'ARDIMS3.EXT'   !  area-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'TMDIMS3.EXT'   !  emis temporal (need for DAYS variable)
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         ENVINT
        LOGICAL         ENVYN
        LOGICAL         GETYN
        CHARACTER*10    HHMMSS
        INTEGER         INDEX1
        CHARACTER*14    MMDDYY
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         STR2INT
        INTEGER         TRIMLEN
        INTEGER         WKDAY
        REAL            YR2DAY

        EXTERNAL        ENVINT, ENVYN, GETYN, HHMMSS, INDEX1, 
     &                  MMDDYY, PROMPTFFILE, PROMPTMFILE, STR2INT,
     &                  TRIMLEN, WKDAY, YR2DAY

C...........   LOCAL VARIABLES and their descriptions:
C...........   Area Sources input and output arrays 

        INTEGER         INTVAL ( NASRC )  !  inventory integer values
        INTEGER         INVYR  ( NASRC )  !  stored inventory year
        REAL            LEAPFAC( NASRC )  !  adjustment factors for leap year
        REAL            REALVAL( NASRC )  !  inventory float values
        REAL            EMISV  ( NASRC )  !  input/output emissions

C...........   Projection/Control Matrix

        REAL            FACS( NASRC,NIPOL )

        INTEGER         NVNAMS
        CHARACTER*16    CPVNAMS( NIPOL ) ! variable names in matrix

C...........   Index from EINAM to CPVNAMS
       INTEGER          CPVIDX ( NIPOL )

C...........   Logical names and unit numbers
        INTEGER         LDEV    !  i/o api log file

        CHARACTER*16    CNAME   !  logical name for projeciton/control file
        CHARACTER*16    ENAME   !  logical name for emission output file
        CHARACTER*16    ANAME   !  logical name for inventory file

C...........   Other local variables

        INTEGER         I, J, L, S, V
        INTEGER         JDATE, JTIME
        INTEGER         NSTEPS
        INTEGER         OUTYR   ! output year (if any difference)
        INTEGER         VTYP    ! tmp i/o api variable type

        REAL            RTMP                ! tmp real value

        LOGICAL         EFLAG               ! Error flag
        DATA            EFLAG / .FALSE. /

        CHARACTER*5     OUTNAM  !  contains default  logical output name
        CHARACTER*16    VBUF    !  tmp i/o api variable buffer
        CHARACTER*256   MESG

C***********************************************************************
C   begin body of program PCAREA

        LDEV = INIT3()

        CALL INITEM( LDEV )

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program PCAREA to take a SMOKE area source inventory file',
     &  'and a projection or control matrix, apply the factors to ',
     &  'the input inventory file and output the result to a ',
     &  'SMOKE file.',
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
            CALL M3EXIT( 'PCAREA', 0, 0, 'Ending Program.', 2 )
        END IF

C.......   Get file names; open input gridding, speciation, 
C.......   optional control matrix files, area source file, 
C.......   optional output file

        ANAME= PROMPTMFILE( 
     &          'Enter logical name for SMOKE AREA INVENTORY file', 
     &          FSREAD3, 'AREA', 'PCAREA' )

        CNAME = PROMPTMFILE( 
     &          'Enter logical name for PROJECTION or CONTROL MATRIX',
     &          FSREAD3, 'APMAT', 'PCAREA' )


C.......   Check dimensions of projection or control matrix and store variable
C          names
        IF ( .NOT. DESC3( CNAME ) ) THEN
            CALL M3EXIT( 'PCAREA', 0, 0, 
     &          'Could not get description of file ' // CNAME, 2 )

        ELSEIF( NROWS3D .NE. NASRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  PROJECTION/CONTROL file:', NROWS3D,
     &      'program:', NASRC
            CALL M3EXIT( 'PCAREA', 0, 0, MESG, 2 )

        ELSEIF( NVARS3D .GT. NIPOL ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch. PROJECTION/CONTROL file has:', NVARS3D,
     &      'variables, which is more than no. of pollutants ', NIPOL
            CALL M3EXIT( 'PCAREA', 0, 0, MESG, 2 )

        ENDIF

        NVNAMS = NVARS3D
        DO 22 V = 1, NVNAMS
            J = INDEX( VNAME3D( V ), '_' )
            L = TRIMLEN( VNAME3D( V ) )
            CPVNAMS( V ) = VNAME3D( V )( J+1:L )
22      CONTINUE

C.........  Build index from master pollutant list to pollutants in matrix
        DO 33 V = 1, NIPOL
            I = INDEX1( EINAM( V ), NVNAMS, CPVNAMS )
            CPVIDX( V ) = I
33      CONTINUE

C.........  Get output year from projection matrix, if any
        OUTYR = 0
        DO 35 I = 1, MXDESC3
            J = INDEX( FDESC3D( I ), '/YEAR/' )
            IF( J .GT. 0 ) THEN
                L = TRIMLEN( FDESC3D( I ) )
                OUTYR = STR2INT( FDESC3D( I )( J+6:L ) )
            END IF
35      CONTINUE
       
C.........  Build description of output file, and optionally open it:

        IF ( .NOT. DESC3( ANAME) ) THEN
            CALL M3EXIT( 'PCAREA', 0, 0, 
     &          'Could not get description of file ' // ANAME, 2 )

        ELSEIF( NROWS3D .NE. NASRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  AREA SOURCES file:', NROWS3D,
     &      'program:', NASRC
            CALL M3EXIT( 'PCAREA', 0, 0, MESG, 2 )

        ENDIF

        JDATE = SDATE3D
        JTIME = STIME3D

C.........  Based on input files and control status so far, set default
C........   output file names

C.........  Annual inventory file as input
        NSTEPS  = 1
        FDESC3D( 11 ) = 'Emissions adjusted using PCAREA'
        OUTNAM = 'AREAC'


        ENAME = PROMPTMFILE( 
     &          'Enter logical name for OUTPUT EMIS file',
     &          FSUNKN3, OUTNAM, UPNAM3D )

C.......   Read and write any variables in the input file that are not 
C          inventory pollutant data
        DO 66 I = 1, NVARS3D

            VBUF = VNAME3D( I )
            VTYP = VTYPE3D( I )
            L = TRIMLEN( VBUF )

            IF( INDEX1( VBUF, NIPOL, EINAM ) .GT. 0 ) GO TO 66

            IF( VTYP .EQ. M3INT ) THEN

                IF( .NOT. READ3( ANAME, VBUF, ALLAYS3, 
     &                           JDATE, JTIME, INTVAL ) ) THEN

                    MESG = 'Could not read "' // VBUF( 1:L ) //
     &                     '" from file ' // ANAME
                    CALL M3EXIT( 'PCAREA', JDATE, JTIME, MESG, 2 )
  
                END IF          !  if read3() failed

            ELSEIF( VTYP .EQ. M3REAL ) THEN

                IF( .NOT. READ3( ANAME, VBUF, ALLAYS3, 
     &                           JDATE, JTIME, REALVAL ) ) THEN

                    MESG = 'Could not read "' // VBUF( 1:L ) //
     &                     '" from file ' // ANAME
                    CALL M3EXIT( 'PCAREA', JDATE, JTIME, MESG, 2 )
  
                END IF          !  if read3() failed

            ENDIF


            IF( VBUF .EQ. 'INVYR' .AND. OUTYR .GT. 0 ) THEN

                DO 55 S = 1, NASRC
                    INVYR ( S ) = INTVAL( S )
                    INTVAL( S ) = OUTYR
55              CONTINUE

            ENDIF

            IF( VTYP .EQ. M3INT ) THEN

                IF ( .NOT. WRITE3( ENAME, VBUF, 
     &                             0, 0, INTVAL ) ) THEN

                    MESG = 'Could not write "' // VBUF( 1:L ) //
     &                     '" to file ' // ENAME
                    CALL M3EXIT( 'PCAREA', JDATE, JTIME, MESG, 2 )
  
                END IF          !  if write3() failed

            ELSEIF( VTYP .EQ. M3REAL ) THEN

                IF ( .NOT. WRITE3( ENAME, VBUF, 
     &                             0, 0, REALVAL ) ) THEN

                    MESG = 'Could not write "' // VBUF( 1:L ) //
     &                     '" to file ' // ENAME
                    CALL M3EXIT( 'PCAREA', JDATE, JTIME, MESG, 2 )
  
                END IF          !  if write3() failed

            ENDIF
 
66      CONTINUE

C.........  Set year-to-year adjustment to account for one of the years
C           being a leap-year, and one of the years not. 
        RTMP = YR2DAY( OUTYR ) 
        DO 77 S = 1, NASRC
            LEAPFAC( S ) = YR2DAY( INVYR( S ) ) / RTMP
77      CONTINUE

C.........   Read projection/control matrix:

        CALL M3MSG2( 'Reading PROJECTION/CONTROL matrix...' )

        IF ( .NOT. READ3( CNAME, 'ALL', 1, 0, 0, FACS ) ) THEN

            CALL M3EXIT( 'PCAREA', 0, 0,
     &           'Could not read gridding matrix from file "' //
     &           CNAME( 1 : TRIMLEN( CNAME ) ) // '".', 2 )

        END IF      !  if read3() failed for gridding matrix

C.......   Transform and write out area source emissions values:

        CALL M3MSG2( 'Calculating area source emissions...' )

C.........  Loop through pollutants, read, and apply projection/control matrix

        DO  198  V = 1, NIPOL

            VBUF = EINAM( V )
            L    = TRIMLEN( VBUF )
            IF( READ3( ANAME, VBUF, ALLAYS3, 
     &                 JDATE, JTIME, EMISV ) ) THEN

                J = CPVIDX( V )
                IF( J .GT. 0 ) THEN

                    DO 101 S = 1, NASRC

                        EMISV( S ) = EMISV( S ) * LEAPFAC( S ) * 
     &                               FACS( S,J )

101                 CONTINUE

                    IF ( .NOT. WRITE3( ENAME, VBUF, 
     &                                 JDATE, JTIME, EMISV ) ) THEN

                        CALL M3EXIT( 'PCAREA', JDATE, JTIME,
     &                           'Could not write "' // VBUF( 1:L ) //
     &                           '" to ' // ENAME, 2 )

                    END IF          !  if write3() failed

                END IF
            ELSE

                CALL M3EXIT( 'PCAREA', JDATE, JTIME, 
     &                       'Could not read "' //  VBUF( 1:L ) //
     &                       '" to ' // ENAME, 2 )

            END IF          !  if read3() succeeds or not

198     CONTINUE            ! end loop on input variables

C.........  Normal completion
        CALL M3EXIT( 'PCAREA', 0, 0, ' ', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

C...........   Formatted file I/O formats............ 93xxx

93020   FORMAT( 8X, 'at time ', A8 )

C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT ( 10 ( A, :, I10, :, 2X ) )

        END

