
        PROGRAM MULTMAT

C***********************************************************************
C  program body starts at line 101
C
C  DESCRIPTION:
C       Read in control matrices and multiply them
C       together to get a final control matrice.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; GETYN, PROMPTMFILE, TRIMLEN
C
C  REVISION  HISTORY:
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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

        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning paramaters
        INCLUDE 'ARDIMS3.EXT'   !  area-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'TMDIMS3.EXT'   !  grid parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         CLOSE3
        LOGICAL         DSCGRID
        INTEGER         ENVINT
        LOGICAL         ENVYN
        LOGICAL         GETYN
        INTEGER         INDEX1
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN

        EXTERNAL        CLOSE3, DSCGRID, ENVINT, ENVYN, GETYN, INDEX1,
     &                  PROMPTMFILE, TRIMLEN 

C...........   PARAMETER

        CHARACTER*16    NONE
        PARAMETER     ( NONE = 'NONE' )

C...........   Control Matrix

        REAL            CC( NPSRC, NIPOL )
        REAL            CMAT( NPSRC, NIPOL ) 

        REAL            CCA( NASRC, NIPOL )   
        REAL            CMATA( NASRC, NIPOL )

C...........   Other local variables

        INTEGER         MSRCS
        INTEGER         IFILE
        INTEGER         LDEV
        INTEGER         S, V

        CHARACTER*16    CNAME   !  logical name for 1st control matrix input file
        CHARACTER*16    MNAME   !  logical name for control matrix output file
        CHARACTER*8     CFILE   !  logical name for next control matrix input file
        CHARACTER*16    SCRBUF  !  scratch buffer
        CHARACTER*256   MESG    !  message buffer

        LOGICAL         PFLAG   !  true for point src, false for area src
        LOGICAL         ENDFLAG !  true when list of control matrices ended 

C***********************************************************************
C   begin body of program MULTMAT

        LDEV = INIT3()

        CALL INITEM( LDEV )

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program MULTMAT takes sorted source level control ',
     &  'matrices and multiply them together and output a new ',
     &  'control matrix.  This program works for either point ',
     &  'or area source control matrices.',
     &  ' '
        WRITE( *,92000 ) 
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").  Use "NONE" as the',
     &  'name for the control-matrix file if you wish to end the list ',
     &  'of control matrices to use.',
     &  ' '
        WRITE( *,92000 ) 
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry.  Default responses are indicated',
     &  'in brackets [LIKE THIS].',
     &  ' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'MULTMAT', 0, 0, 'Ending Program', 2 )
        END IF

C.......   Get environment variable settings

C.......   Initialize control matrix array

        DO 10 S = 1, NIPOL

           DO 20 V = 1, NPSRC

              CMAT( V,S ) = 1.00 

20         CONTINUE

           DO 30 V = 1, NASRC
 
             CMATA( V,S ) = 1.00

30         CONTINUE

 
10      CONTINUE

C.......   Get file names; open  
C.......   control matrix files 


        CNAME = PROMPTMFILE( 
     &          'Enter logical name for first CONTROL MATRIX file ' //
     &          'or "NONE" to',
     &          FSREAD3, 'PCMAT1', 'MULTMAT' )

        IFILE = 1

C.......   Read description of input file 
C.......   Check to see if this is for point or area sources

        CALL M3MSG2( 'Reading CONTROL matrix...' )

        IF ( .NOT. DESC3( CNAME ) ) THEN
                WRITE( MESG, '(2A)' )
     &          'Could not get description of file ', CNAME
                CALL M3EXIT( 'MULTMAT', 0, 0, MESG, 2 )
        ENDIF

        IF ( NROWS3D .EQ. NPSRC ) THEN
 
          PFLAG = .TRUE.
          MSRCS = NPSRC

        ELSE IF ( NROWS3D .EQ. NASRC ) THEN
 
          PFLAG = .FALSE.
          MSRCS = NASRC
  
        ELSE 
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  C-MATRIX file:', NROWS3D,
     &          'program NPSRC :', NPSRC, ' program NASRC :',
     &          NASRC
      
                CALL M3EXIT( 'MULTMAT', 0, 0, MESG, 2 )
        END IF

        ENDFLAG = .FALSE.

C....  Do while NONE used for name of control matrix filename

        DOWHILE ( .NOT. ENDFLAG  )

C....  Read control factors from input file(s)
C....  if point source control matrices are used...
 
           IF ( PFLAG ) THEN
 
             DO  66 S = 1, NIPOL
               
                SCRBUF = 'CONTROL_' // EINAM( S )
                IF ( .NOT. READ3( CNAME, SCRBUF, 1, 0, 0, CC( 1,S ) ) 
     &             ) THEN

                    CALL M3ERR( 'MULTMAT', 0, 0,
     &                   'Could not read "' // SCRBUF( 1:12 ) // 
     &                   '" from CONTROL FILE "' //
     &                   CNAME( 1 : TRIMLEN( CNAME ) ) // '".',
     &                   .TRUE. )

                 END IF      !  if read3() failed for control matrix

66           CONTINUE

C......   Multiply matrices to get updated factors

             CALL M3MSG2( 'Updating factors ...' )


             DO 80 V = 1, MSRCS

              DO 70 S = 1, NIPOL

                CMAT( V,S ) = CMAT( V,S ) * CC( V,S ) 

70            CONTINUE

80           CONTINUE

           ELSE

C..........  if area source control matrices are input.

             DO  86 S = 1, NIPOL

                SCRBUF = 'CONTROL_' // EINAM( S )
                IF ( .NOT. READ3( CNAME, SCRBUF, 1, 0, 0, CCA( 1,S ) )
     &             ) THEN

                    CALL M3ERR( 'MULTMAT', 0, 0,
     &                   'Could not read "' // SCRBUF( 1:12 ) //
     &                   '" from CONTROL FILE "' //
     &                   CNAME( 1 : TRIMLEN( CNAME ) ) // '".',
     &                   .TRUE. )

                 END IF      !  if read3() failed for control matrix

86           CONTINUE

C......   Multiply matrices to get updated factors

             CALL M3MSG2( 'Updating factors ...' )

             DO 90 V = 1, MSRCS

              DO 95 S = 1, NIPOL

                CMATA( V,S ) = CMATA( V,S ) * CCA( V,S )

95            CONTINUE

90          CONTINUE

          ENDIF

          IF ( .NOT. CLOSE3( CNAME ) ) THEN
                     CALL M3ERR( 'MULTMAT', 0, 0,
     &                   CNAME( 1 : TRIMLEN( CNAME ) ) // '".',
     &                   .TRUE. )

          END IF      !  if close3() failed for control matrix
  

C......   Check to see if there is another matrice to be read in  

          IFILE = IFILE + 1
 
          WRITE ( CFILE, 94020 ) 'PCMAT' , IFILE
          CNAME = PROMPTMFILE(
     &          'Enter logical name for next CONTROL MATRIX file ' //
     &          'or "NONE" to end list',
     &          FSREAD3, CFILE , 'MULTMAT' )

          CALL M3MSG2( 'Reading CONTROL matrix...' )
           
          IF ( CNAME .NE. NONE ) THEN
            IF ( .NOT. DESC3( CNAME ) ) THEN
                WRITE( MESG, '(2A)' )
     &          'Could not get description of file ', CNAME
                CALL M3EXIT( 'MULTMAT', 0, 0, MESG, 2 )
            ELSE IF ( NROWS3D .NE. MSRCS ) THEN
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  C-MATRIX file:', NROWS3D,
     &          'program:', MSRCS
                CALL M3EXIT( 'MULTMAT', 0, 0, MESG, 2 )
            END IF

          ELSE

            ENDFLAG = .TRUE.

          ENDIF

        ENDDO   ! end dowhile no more matrices to be read in

C.... Open output file : new control matrix

        MNAME = PROMPTMFILE(
     &          '...for CONTROL MATRIX output file',
     &          FSUNKN3, 'PCMAT', 'MULTMAT' )

C............................................................................
C.......   Write out the control matrix:

900    CONTINUE          !  Done processing; write output file:

        IF ( PFLAG ) THEN

          WRITE( *,92000 )' ','Writing new point CONTROL MATRIX...',' '

          IF ( .NOT. WRITE3( MNAME, 'ALL', 0, 0,CMAT ) ) THEN
            CALL M3EXIT( 'MULTMAT', 0, 0,
     &                   'Error writing CONTROL MATRIX file.', 2 )
          END IF
 
        ELSE

          WRITE( *,92000 )' ','Writing new area CONTROL MATRIX...',' '

          IF ( .NOT. WRITE3( MNAME, 'ALL', 0, 0,CMATA ) ) THEN
            CALL M3EXIT( 'MULTMAT', 0, 0,
     &                   'Error writing CONTROL MATRIX file.', 2 )
          END IF

        ENDIF

999     CONTINUE          !  exit program

        CALL M3EXIT( 'MULTMAT', 0, 0, 
     &               'Normal completion  of PROGRAM MULTMAT', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )


C...........   Formatted file I/O formats............ 93xxx

93010   FORMAT( A16 )

93020   FORMAT( 8X, 'at time ', A8 )
 

C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT( 10( A, :, I7, :, 1X ) )

94020   FORMAT( A5, I1 ) 

C...........   Miscellaneous formats................. 95xxx

95000   FORMAT ( /5X , A , $ )          !  generic prompt format.


        END

