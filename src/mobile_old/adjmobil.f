
        PROGRAM ADJMOBIL

C***********************************************************************
C  program body starts at line 
C
C  DESCRIPTION:
C       Adjusts variables in the MOBL or MTMP file and outputs another file
C       NOTE: Initial implementation limited to MTMP and road class adjustment
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; PROMPTMFILE
C
C  REVISION  HISTORY:
C       Prototype 1/99 by M Houyoux
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
        
        INCLUDE 'MBDIMS3.EXT'   !  mobil-source dimensioning parameters
        INCLUDE 'TMDIMS3.EXT'   !  temporal dimensioning parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        INTEGER         FIND1
        LOGICAL         GETYN
        CHARACTER*10    HHMMSS
        INTEGER         INDEX1
        CHARACTER*14    MMDDYY
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN
        INTEGER         WKDAY           !  day of week (1...7)

        EXTERNAL        FIND1, GETYN, HHMMSS, INDEX1, MMDDYY, 
     &                  PROMPTFFILE, PROMPTMFILE, TRIMLEN, WKDAY

C...........   LOCAL VARIABLES and their descriptions:
C...........   Mobil Sources input and output arrays (one variable at a time)
        
        REAL            EMISV( NMSRC, NVTYPE )  !  input emissions.
        
C...........   Inventory arrays
        INTEGER         IRCLAS( NMSRC )   ! Inventory year code for each source
 
C...........   Factor arrays
        INTEGER ADJRDT( NRCLAS )

        REAL    ADJFAC ( NRCLAS )
        REAL    OLDEMIS( NVTYPE, NRCLAS )
        REAL    NEWEMIS( NVTYPE, NRCLAS )

C...........   Indicator for which variable to adjust
        LOGICAL LADJUST( MXVARS3 )

C...........   Other local variables
        INTEGER         I, J, K, L, S, T, V

        INTEGER         IOS
        INTEGER         JDATE, JTIME, TSTEP
        INTEGER         LDATE, LTIME
        INTEGER         LDEV
        INTEGER         LRDT               ! previous road class code
        INTEGER         NSTEPS
        INTEGER         MDEV               ! input factors unit number
        INTEGER         RDT                ! temporary road class code

        REAL            FAC  ! temporary adjustment factor

        CHARACTER*7     OUTNAM  !  contains default  logical output name
        CHARACTER*16    ANAME   !  logical name for mobile-source input file
        CHARACTER*16    INAME   !  logical name for input inventory name
        CHARACTER*16    ONAME   !  logical name for output file name
        CHARACTER*16    POLADJ  !  pollutant name to adjust

        CHARACTER*256   MESG    !  scratch message area

C***********************************************************************
C   begin body of program ADJMOBIL
        
        LDEV = INIT3()

        CALL INITEM( LDEV )
        
        WRITE( *,92000 )
     &  ' ',
     &  'Program ADJMOBIL to adjust mobile source emissions by road',
     &  'class by a percentage in the MTMP file.',
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
            CALL M3EXIT( 'ADJMOBIL', 0, 0, 'Ending Program', 2 )
        END IF

C.......   Get environment variable string of pollutant that needs adjusting.
        MESG = 'Name of pollutant to adjust'
        CALL ENVSTR( 'ADJUST_POLLUTANT', MESG, 'NOX', POLADJ, IOS )
 
C.......   Get file names; open input mobile source, mods table, and
C.......   output files.
        INAME = PROMPTMFILE( 
     &          'Enter logical name for MOBILE INVENTORY ' // 
     &          'file', FSREAD3, 'MOBL', 'ADJMOBIL' )

        ANAME = PROMPTMFILE( 
     &          'Enter logical name for TIME-STEPPED ' // 
     &          'MOBILE EMIS file', FSREAD3, 'MTMP', 'ADJMOBIL' )

C.......   Open table for modifying input 
        MDEV = PROMPTFFILE( 
     &           'Enter logical name for adjustments file',
     &           .TRUE., .TRUE., 'MAJFACS', 'ADJMOBIL' )

C.......   Build description of output file, and optionally open it:

        IF ( .NOT. DESC3( INAME ) ) THEN
            CALL M3EXIT( 'ADJMOBIL', 0, 0, 
     &              'Could not get description of file ' // INAME, 2 )

        ELSEIF( NROWS3D .NE. NMSRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  MOBILE INVEN file:', NROWS3D,
     &      'program:', NMSRC
            CALL M3EXIT( 'ADJMOBIL', 0, 0, MESG, 2 )
 
        END IF

        IF ( .NOT. DESC3( ANAME ) ) THEN
            CALL M3EXIT( 'ADJMOBIL', 0, 0, 
     &              'Could not get description of file ' // ANAME, 2 )

        ELSEIF( NROWS3D .NE. NMSRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  MOBILE SOURCES file:', NROWS3D,
     &      'program:', NMSRC
            CALL M3EXIT( 'ADJMOBIL', 0, 0, MESG, 2 )
 
        END IF

        JDATE   = SDATE3D
        JTIME   = STIME3D
        TSTEP   = TSTEP3D
        NSTEPS  = MXREC3D

C.........  Annual inventory VMT file as input
        IF ( TSTEP .EQ. 0 ) THEN
            MESG = 'Cannot run using MOBL file!'
            CALL M3EXIT( 'ADJMOBIL', 0, 0, MESG, 2 )
        END IF

C.......   Read adjustment factors
        LRDT = -9
        DO 11 I = 1, NRCLAS
            READ( MDEV, *, END=12, IOSTAT=IOS ) RDT, FAC

            IF( IOS .GT. 0 ) THEN
                WRITE( MESG,94010 ) 'I/O error', IOS,
     &                 'reading factors at line', I
                CALL M3EXIT( 'ADJMOBIL', 0, 0, MESG, 2 )
            ELSEIF( RDT .LE. LRDT ) THEN

                WRITE( MESG,94010 ) 'Road classes are not sorted ' //
     &                 'at line', I, 'of factors file.'
                CALL M3EXIT( 'ADJMOBIL', 0, 0, MESG, 2 )
            ENDIF

            LRDT = RDT

            ADJRDT( I ) = RDT
            ADJFAC( I ) = FAC
           
11      CONTINUE
12      CONTINUE

C.........  Read road class codes
        IF( .NOT. READ3( INAME, 'IRCLAS', ALLAYS3, 
     &                   0, 0, IRCLAS )         ) THEN

            MESG = 'Could not read "IRCLAS" from ' // INAME
            CALL M3EXIT( 'ADJMOBIL', JDATE, JTIME, MESG, 2 )

        END IF          !  if read3() succeeds or not

C.........  Set up to open output file
        L = TRIMLEN( FDESC3D( 1 ) )
        FDESC3D( 1 ) = FDESC3D( 1 )( 1:L ) // ' adjusted by ADJMOBIL'
        IF ( TSTEP .EQ. 0 ) THEN
            OUTNAM = 'MOBLOUT'

            DO 81 I = 1, NVARS3D
                LADJUST( I ) = .FALSE.
81          CONTINUE
  
            K = INDEX1( 'VMT', NVARS3D, VNAME3D )
            LADJUST( K ) = .TRUE.
            MESG = 'NOTE: Adjusting variable VMT'
            CALL M3MSG2( MESG )

        ELSE  
            OUTNAM = 'MTMPOUT'

C.............  Scan the input variable names and flag the ones that
C               contain the name of the pollutant to adjust.
            L = TRIMLEN( POLADJ )
            DO 91 I = 1, NVARS3D
                K = INDEX( VNAME3D( I ), POLADJ( 1:L ) )
                IF( K .GT. 0 ) THEN
                    LADJUST( I ) = .TRUE.
                    MESG = 'NOTE: Adjusting variable ' // VNAME3D( I )
                    CALL M3MSG2( MESG )
                ELSE
                    LADJUST( I ) = .FALSE.
                END IF
91          CONTINUE

        END IF

        ONAME = PROMPTMFILE( 'Enter logical name for output file',
     &                       FSUNKN3, OUTNAM, 'ADJMOBIL' )

C.......   Initialize old and new total emissions by road clas
       DO 103 J = 1, NRCLAS
           DO 101 I = 1, NVTYPE
               OLDEMIS( I, J ) = 0.
               NEWEMIS( I, J ) = 0.
101        CONTINUE
103    CONTINUE

C.......   Transform and write out mobile source emissions values:
        LDATE = 0
        DO  199  T = 1, NSTEPS

C...............   If this is a new month, or new day, write message
            IF ( JDATE .GT. 0 .AND. LDATE. NE. JDATE ) THEN
 
                MESG = 'Processing ' //
     &                 DAYS( WKDAY( JDATE ) ) // MMDDYY( JDATE )
                CALL M3MSG2( MESG( 1:TRIMLEN( MESG ) ) )
 
            END IF

C.............  Write to screen because WRITE3 only writes to LDEV
            IF( JDATE .GT. 0 ) WRITE( *, 93020 ) HHMMSS( JTIME )

C.............   Read in entire emissions field
            DO 188 V = 1, NVARS3D

                IF( .NOT. READ3( ANAME, VNAME3D( V ), ALLAYS3, 
     &                           JDATE, JTIME, EMISV )         ) THEN

                    MESG = 'Could not read "' //
     &                     VNAME3D( V )( 1:TRIMLEN( VNAME3D( V ) ) ) //
     &                     '" from ' // ANAME
                    CALL M3EXIT( 'ADJMOBIL', JDATE, JTIME, MESG, 2 )

                END IF          !  if read3() succeeds or not

C.................  If this is the variable to adjust, then do it!
                IF( LADJUST( V ) ) THEN

                    DO 111 K = 1, NVTYPE
                        DO 109 S = 1, NMSRC

                            RDT = IRCLAS( S )

                            J = FIND1( RDT, NRCLAS, ADJRDT )

                            OLDEMIS( K,J ) = OLDEMIS( K,J )+EMISV( S,K )

                            EMISV( S,K ) = EMISV( S,K ) * ADJFAC( J )

                            NEWEMIS( K,J ) = NEWEMIS( K,J )+EMISV( S,K )

109                     CONTINUE
111                 CONTINUE

                ENDIF

                IF ( .NOT. WRITE3( ONAME, VNAME3D( V ),
     &                             JDATE, JTIME, EMISV  ) ) THEN

                    MESG = 'Could not write "' //
     &                     VNAME3D( V )( 1:TRIMLEN( VNAME3D( V ) ) )
     &                     // '" to ' // ONAME
                    CALL M3EXIT( 'ADJMOBIL', JDATE, JTIME, MESG, 2 )

                END IF          !  if write3() failed

188         CONTINUE            ! end loop on input variables

            LDATE = JDATE
            LTIME = JTIME

            CALL NEXTIME( JDATE, JTIME, TSTEP )

199     CONTINUE          !  end loop on time steps

C.........  Write old and new emissions by road class
        DO 203 J = 1, NRCLAS
            DO 201 I = 1, NVTYPE
                WRITE( LDEV, * ) ' RDT=', ADJRDT( J ), 
     &                           ' VTYPE=', VTYPE3( I ),
     &                           ' Old emis=', OLDEMIS( I,J ),
     &                           ' New emis=', NEWEMIS( I,J ) 
201         CONTINUE
203     CONTINUE
      
        CALL M3EXIT( 'ADJMOBIL', 0, 0, 
     &               'Normal completion of PROGRAM ADJMOBIL', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

C...........   Formatted file I/O formats............ 93xxx
 
93020   FORMAT( 8X, 'at time ', A8 )
 

C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT ( 10 ( A, :, I10, :, 2X ) )

        END

