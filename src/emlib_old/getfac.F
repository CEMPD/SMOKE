
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE GETFACR( NAME, DEFAULT, NMAX, N, 
     &                      INDEX, NCOL, RTAB1 )

C***********************************************************************
C  program body starts at line
C
C  DESCRIPTION:
C      Program prompts, opens, and reads an ascii file with an array 
C      of integer indices, and real table of indeterminant size
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Prototype 12/95 by MRH
C
C***********************************************************************

        IMPLICIT NONE

        INCLUDE 'PARMS3.EXT'     !  I/O API parameters
        INCLUDE 'IODECL3.EXT'    !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'     !  I/O API file description data structures

C.........  External
        INTEGER  GETEFILE
        INTEGER  TRIMLEN
        EXTERNAL GETEFILE, TRIMLEN

C.........  Program arguments
        INTEGER N, NCOL, NMAX
        INTEGER INDEX( NMAX )       ! Sorted profile ID codes

        REAL    RTAB1( NCOL, NMAX ) ! Sorted profile of real numbers
 
        CHARACTER*(*) NAME, DEFAULT

C.........  Parameters
        INTEGER  CMAX
        PARAMETER  ( CMAX = 200 )   ! Constant maximum thanks to FORTRAN

C.........  Local variables
        INTEGER  INDEXA( CMAX )     ! Unsorted profile ID codes
        INTEGER  MOVETO( CMAX )     ! Index for rearranging profiles
        INTEGER  SORTER( CMAX )     ! Index for sorting profile ID codes

        INTEGER IDEV, IREC, J1, K   ! Counters and pointers
        INTEGER IDUM                ! Dummy integer code
        INTEGER IOS                 ! Read status
        INTEGER READTHRU            ! Loop counter for reading through file
        REAL    RPDUM               ! Dummy real profile
        LOGICAL EFLAG               ! Error flag
        DATA    EFLAG  / .FALSE. /

        CHARACTER*256 BUFFER
        CHARACTER*256 MESG

C***********************************************************************
C   begin body of program GETFACR

C......... Set buffer for writing messages
        BUFFER = NAME

C.........  Open input file
        IDEV = GETEFILE( DEFAULT, .TRUE., .TRUE., 'GETFACR' )

        IF ( IDEV .LT. 0 ) THEN     !  failure to open
 
            MESG = 'Could not open input file "' //
     &             DEFAULT( 1 : TRIMLEN( DEFAULT ) ) // '".'
            CALL M3EXIT( 'GETFACR', 0, 0, MESG, 2 )

        ELSE

            CALL M3MSG2( 'NOTE: ' // BUFFER ( 1:TRIMLEN( BUFFER  ) ) // 
     &                   ' file "' // DEFAULT( 1:TRIMLEN( DEFAULT ) ) //
     &                   '" opened successfully.' )
 
        END IF      !  if getefile() failed

C.........  First time through file, establish sorted order, second time
C.........  through file, store profiles in that sorted order
        DO 44 READTHRU = 1, 2

            IREC = 0
   11       CONTINUE

                IREC = IREC + 1

                IF( IREC .GT. NMAX ) THEN

                    WRITE( MESG, 94010 ) 
     &                   'Maximum dimension =', NMAX,
     &                   'was exceeded while reading '// 
     &                   BUFFER( 1:TRIMLEN( BUFFER ) )  // ' file.'
                    CALL  M3EXIT( 'GETFACR', 0, 0, MESG, 2 )

                ELSEIF( IREC .GT. CMAX ) THEN

                    WRITE( MESG, 94010 ) 
     &                   'Code parameter CMAX (=', CMAX,
     &                   ') was exceeded while reading '// 
     &                   BUFFER( 1:TRIMLEN( BUFFER ) )  // ' file.'
                    CALL  M3MESG( MESG )
                    CALL  M3EXIT( 'GETFACR', 0, 0, 
     &                      'Change CMAX in GETFACR and recompile', 2 )

                ELSE

                    IF( READTHRU .EQ. 1 ) THEN
                      READ( IDEV, *, END=33, IOSTAT=IOS ) 
     &                      INDEXA( IREC ), ( RPDUM, J1=1,NCOL )

                      SORTER( IREC ) = IREC

                    ELSE
                      READ( IDEV, *, END=33, IOSTAT=IOS ) IDUM,
     &                    ( RTAB1( J1, MOVETO( IREC )), J1=1,NCOL )

                    ENDIF
    
                    IF( IOS .NE. 0 ) THEN

                        EFLAG = .TRUE.
                        WRITE( MESG, 94010 )
     &                     'I/O Error:', IOS, ', while reading ' //
     &                     BUFFER( 1:TRIMLEN( BUFFER ) ) // 
     &                     ' file record ', IREC
                        CALL M3MESG( MESG )

                    ENDIF

                ENDIF

            GO TO 11

   33       IF( EFLAG ) 
     &      CALL M3EXIT( 'GETFACR', 0, 0, 'Error in file format', 2 )

            IF( READTHRU .EQ. 1 ) THEN

                N = IREC - 1

                CALL SORTI1( N, SORTER, INDEXA )

                DO 41 K = 1, N
                    INDEX( K ) = INDEXA( SORTER( K ) )
                    MOVETO( SORTER( K ) ) = K
   41           CONTINUE

                REWIND( IDEV )

            ENDIF

   44   CONTINUE

        RETURN

94010   FORMAT( 10 (A, :, 1X, I5, :, 1X ) )

        END

        SUBROUTINE GETFACRR( NAME, DEFAULT, NMAX, N, 
     &                      INDEX, NCOL1, RTAB1, NCOL2, RTAB2 )

C***********************************************************************
C  program body starts at line 47
C
C  DESCRIPTION:
C      Program prompts, opens, and reads an ascii file with an array 
C      of integer indices, and two real tables of indeterminant size
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Prototype 12/95 by MRH
C
C***********************************************************************

        IMPLICIT NONE

        INCLUDE 'PARMS3.EXT'     !  I/O API parameters
        INCLUDE 'IODECL3.EXT'    !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'     !  I/O API file description data structures

C.........  External
        INTEGER  GETEFILE
        INTEGER  TRIMLEN
        EXTERNAL GETEFILE, TRIMLEN

C.........  Program arguments
        INTEGER N, NCOL1, NCOL2, NMAX
        INTEGER INDEX( NMAX )        ! Sorted profile ID codes

        REAL    RTAB1( NCOL1, NMAX ) ! Sorted 1st profile of real numbers
        REAL    RTAB2( NCOL2, NMAX ) ! Sorted 2nd profile of real numbers
 
        CHARACTER*(*) NAME, DEFAULT

C.........  Parameters
        INTEGER  CMAX
        PARAMETER  ( CMAX = 200 )   ! Constant maximum thanks to FORTRAN

C.........  Local variables
        INTEGER  INDEXA( CMAX )     ! Unsorted profile ID codes
        INTEGER  MOVETO( CMAX )     ! Index for rearranging profiles
        INTEGER  SORTER( CMAX )     ! Index for sorting profile ID codes

        INTEGER IDEV, IREC, J1, J2, K   ! Counters and pointers
        INTEGER IDUM                ! Dummy integer code
        INTEGER IOS                 ! Read status
        INTEGER READTHRU            ! Loop counter for reading through file
        REAL    RPDUM               ! Dummy real profile
        LOGICAL EFLAG               ! Error flag
        DATA    EFLAG  / .FALSE. /

        CHARACTER*256 BUFFER
        CHARACTER*256 MESG

C***********************************************************************
C   begin body of program GETFACRR

C......... Set buffer for writing messages
        BUFFER = NAME

C.........  Open input file
        IDEV = GETEFILE( DEFAULT, .TRUE., .TRUE., 'GETFACRR' )
 
        IF ( IDEV .LT. 0 ) THEN     !  failure to open
 
            MESG = 'Could not open input file "' //
     &             DEFAULT( 1 : TRIMLEN( DEFAULT ) ) // '".'
            CALL M3EXIT( 'GETFACRR', 0, 0, MESG, 2 )
 
        ELSE
 
            CALL M3MSG2( 'NOTE: ' // BUFFER ( 1:TRIMLEN( BUFFER  ) ) // 
     &                   ' file "' // DEFAULT( 1:TRIMLEN( DEFAULT ) ) //
     &                   '" opened successfully.' )
 
        END IF      !  if getefile() failed

C.........  First time through file, establish sorted order, second time
C.........  through file, store profiles in that sorted order
        DO 44 READTHRU = 1, 2

            IREC = 0
   11       CONTINUE

                IREC = IREC + 1

                IF( IREC .GT. NMAX ) THEN

                    WRITE( MESG, 94010 ) 
     &                   'Maximum dimension =', NMAX,
     &                   'was exceeded while reading '// 
     &                   BUFFER( 1:TRIMLEN( BUFFER ) ) // ' file.'
                    CALL  M3EXIT( 'GETFACRR', 0, 0, MESG, 2 )

                ELSEIF( IREC .GT. CMAX ) THEN

                    WRITE( MESG, 94010 ) 
     &                   'Code parameter CMAX (=', CMAX,
     &                   ') was exceeded while reading '// 
     &                   BUFFER( 1:TRIMLEN( BUFFER ) )  // ' file.'
                    CALL  M3MESG( MESG )
                    CALL  M3EXIT( 'GETFACRR', 0, 0, 
     &                      'Change CMAX in GETFACRR and recompile', 2 )

                ELSE

                    IF( READTHRU .EQ. 1 ) THEN
                      READ( IDEV, *, END=33, IOSTAT=IOS ) 
     &                    INDEXA( IREC ), ( RPDUM, J1=1,NCOL1 ),
     &                                    ( RPDUM, J2=1,NCOL2 )

                      SORTER( IREC ) = IREC

                    ELSE
                      READ( IDEV, *, END=33, IOSTAT=IOS ) IDUM,
     &                    ( RTAB1( J1, MOVETO( IREC )), J1=1,NCOL1 ),
     &                    ( RTAB2( J2, MOVETO( IREC )), J2=1,NCOL2 ) 

                    ENDIF


                    IF( IOS .NE. 0 ) THEN

                        EFLAG = .TRUE.
                        WRITE( MESG, 94010 )
     &                     'I/O Error:', IOS, ', while reading ' //
     &                     BUFFER( 1:TRIMLEN( BUFFER ) ) //
     &                     ' file record ', IREC
                        CALL M3MESG( MESG )

                    ENDIF

                ENDIF

            GO TO 11

   33       IF( EFLAG )
     &      CALL M3EXIT( 'GETFACRR', 0, 0, 'Error in file format', 2 )

            IF( READTHRU .EQ. 1 ) THEN

                N = IREC - 1

                CALL SORTI1( N, SORTER, INDEXA )

                DO 41 K = 1, N
                    INDEX( K ) = INDEXA( SORTER( K ) )
                    MOVETO( SORTER( K ) ) = K
   41           CONTINUE

                REWIND( IDEV )

            ENDIF

   44   CONTINUE

        RETURN

94010   FORMAT( 10 (A, :, 1X, I5, :, 1X ) )

        END

        SUBROUTINE GETFACI( NAME, DEFAULT, NMAX, N, 
     &                        INDEX, NCOL1, ITAB1  )

C***********************************************************************
C  program body starts at line
C
C  DESCRIPTION:
C      Program prompts, opens, and reads an ascii file with an array 
C      of integer indices, and then each an integer, then character, then
C      real fields of indeterminant size.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Prototype 12/95 by MRH
C
C***********************************************************************

        IMPLICIT NONE

        INCLUDE 'PARMS3.EXT'     !  I/O API parameters
        INCLUDE 'IODECL3.EXT'    !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'     !  I/O API file description data structures

C.........  EXTERNAL
        INTEGER  GETEFILE
        INTEGER  TRIMLEN
        EXTERNAL GETEFILE, TRIMLEN

C.........  Program arguments
        INTEGER       N, NCOL1, NMAX

        INTEGER       INDEX( NMAX )        ! Sorted profile ID codes
        INTEGER       ITAB1( NCOL1, NMAX ) ! Sorted profile of integers

        CHARACTER*(*) NAME, DEFAULT

C.........  Parameters
        INTEGER  CMAX
        PARAMETER  ( CMAX = 200 )   ! Constant maximum thanks to FORTRAN

C.........  Local variables
        INTEGER  INDEXA( CMAX )     ! Unsorted profile ID codes
        INTEGER  MOVETO( CMAX )     ! Index for rearranging profiles
        INTEGER  SORTER( CMAX )     ! Index for sorting profile ID codes

        INTEGER IDEV, IREC, J1, K   ! Counters and pointers
        INTEGER IDUM                ! Dummy integer code
        INTEGER IPDUM               ! Dummy integer profile
        INTEGER IOS                 ! Read status
        INTEGER READTHRU            ! Loop counter for reading through file
        LOGICAL EFLAG               ! Error flag
        DATA    EFLAG  / .FALSE. /

        CHARACTER*256 BUFFER
        CHARACTER*256 MESG

C***********************************************************************
C   begin body of program GETFACI

C......... Set buffer for writing messages
        BUFFER = NAME

C.........  Open input file
        IDEV = GETEFILE( DEFAULT, .TRUE., .TRUE., 'GETFACI' )
 
        IF ( IDEV .LT. 0 ) THEN     !  failure to open
 
            MESG = 'Could not open input file "' //
     &             DEFAULT( 1 : TRIMLEN( DEFAULT ) ) // '".'
            CALL M3EXIT( 'GETFACI', 0, 0, MESG, 2 )
 
        ELSE
 
            CALL M3MSG2( 'NOTE: ' // BUFFER ( 1:TRIMLEN( BUFFER  ) ) // 
     &                   ' file "' // DEFAULT( 1:TRIMLEN( DEFAULT ) ) //
     &                   '" opened successfully.' )
 
        END IF      !  if getefile() failed

C.........  First time through file, establish sorted order, second time
C.........  through file, store profiles in that sorted order
        DO 44 READTHRU = 1, 2

            IREC = 0
   11       CONTINUE

                IREC = IREC + 1

                IF( IREC .GT. NMAX ) THEN

                    WRITE( MESG, 94010 ) 
     &                   'Maximum dimension =', NMAX,
     &                   'was exceeded while reading '// 
     &                   BUFFER( 1:TRIMLEN( BUFFER ) ) // ' file.'
                    CALL M3EXIT( 'GETFACI', 0, 0, MESG, 2 )

                ELSEIF( IREC .GT. CMAX ) THEN

                    WRITE( MESG, 94010 ) 
     &                   'Code parameter CMAX (=', CMAX,
     &                   ') was exceeded while reading '// 
     &                   BUFFER( 1:TRIMLEN( BUFFER ) )  // ' file.'
                    CALL  M3MESG( MESG )
                    CALL  M3EXIT( 'GETFACI', 0, 0, 
     &                      'Change CMAX in GETFACI and recompile', 2 )

                ELSE

                    IF( READTHRU .EQ. 1 ) THEN
                      READ( IDEV, *, END=33, IOSTAT=IOS ) 
     &                    INDEXA( IREC ), ( IPDUM, J1=1,NCOL1 )

                      SORTER( IREC ) = IREC

                    ELSE
                      READ( IDEV, *, END=33, IOSTAT=IOS ) IDUM,
     &                    ( ITAB1( J1, MOVETO( IREC )), J1=1,NCOL1 )

                    ENDIF

                    IF( IOS .NE. 0 ) THEN

                        EFLAG = .TRUE.
                        WRITE( MESG, 94010 )
     &                     'I/O Error:', IOS, ', while reading ' //
     &                     BUFFER( 1:TRIMLEN( BUFFER ) ) //
     &                     ' file record ', IREC
                        CALL M3MESG( MESG )

                    ENDIF

                ENDIF

            GO TO 11

   33       IF( EFLAG )
     &      CALL M3EXIT( 'GETFACI', 0, 0, 'Error in file format', 2 )

            IF( READTHRU .EQ. 1 ) THEN

                N = IREC - 1

                CALL SORTI1( N, SORTER, INDEXA )

                DO 41 K = 1, N
                    INDEX( K ) = INDEXA( SORTER( K ) )
                    MOVETO( SORTER( K ) ) = K
   41           CONTINUE

                REWIND( IDEV )

            ENDIF

   44   CONTINUE

        RETURN

94010   FORMAT( 10 (A, :, 1X, I5, :, 1X ) )

        END

        SUBROUTINE GETFACICR( NAME, DEFAULT, NMAX, N, 
     &                        INDEX, NCOL1, ITAB1, NCOL2, CTAB2,
     &                        NCOL3, RTAB3 )

C***********************************************************************
C  program body starts at line
C
C  DESCRIPTION:
C      Program prompts, opens, and reads an ascii file with an array 
C      of integer indices, and then each an integer, then character, then
C      real fields of indeterminant size.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Prototype 12/95 by MRH
C
C***********************************************************************

        IMPLICIT NONE

        INCLUDE 'PARMS3.EXT'     !  I/O API parameters
        INCLUDE 'IODECL3.EXT'    !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'     !  I/O API file description data structures

C.........  EXTERNAL
        INTEGER  GETEFILE
        INTEGER  TRIMLEN
        EXTERNAL GETEFILE, TRIMLEN

C.........  Program arguments
        INTEGER       N, NCOL1, NCOL2, NCOL3, NMAX

        INTEGER       INDEX( NMAX )
        INTEGER       ITAB1( NCOL1, NMAX )

        REAL          RTAB3( NCOL3, NMAX )
 
        CHARACTER*(*) CTAB2( NCOL2, NMAX )
        CHARACTER*(*) NAME, DEFAULT

C.........  Parameters
        INTEGER  CMAX
        PARAMETER  ( CMAX = 200 )   ! Constant maximum thanks to FORTRAN

C.........  Local variables
        INTEGER  INDEXA( CMAX )     ! Unsorted profile ID codes
        INTEGER  MOVETO( CMAX )     ! Index for rearranging profiles
        INTEGER  SORTER( CMAX )     ! Index for sorting profile ID codes

        INTEGER IDEV, IREC, J1, J2, J3, K   ! Counters and pointers
        INTEGER IDUM                ! Dummy integer code
        INTEGER IPDUM               ! Dummy integer profile
        INTEGER IOS                 ! Read status
        INTEGER READTHRU            ! Loop counter for reading through file
        REAL    RPDUM               ! Dummy real profile
        LOGICAL EFLAG               ! Error flag
        DATA    EFLAG  / .FALSE. /

        CHARACTER*256 BUFFER
        CHARACTER*256 MESG

C***********************************************************************
C   begin body of program GETFACICR

C......... Set buffer for writing messages
        BUFFER = NAME

C.........  Open input file
        IDEV = GETEFILE( DEFAULT, .TRUE., .TRUE., 'GETFACICR' )
 
        IF ( IDEV .LT. 0 ) THEN     !  failure to open
 
            MESG = 'Could not open input file "' //
     &             DEFAULT( 1 : TRIMLEN( DEFAULT ) ) // '".'
            CALL M3EXIT( 'GETFACICR', 0, 0, MESG, 2 )
 
        ELSE
 
            CALL M3MSG2( 'NOTE: ' // BUFFER ( 1:TRIMLEN( BUFFER  ) ) // 
     &                   ' file "' // DEFAULT( 1:TRIMLEN( DEFAULT ) ) //
     &                   '" opened successfully.' )
 
        END IF      !  if getefile() failed

C.........  First time through file, establish sorted order, second time
C.........  through file, store profiles in that sorted order
        DO 44 READTHRU = 1, 2

            IREC = 0
   11       CONTINUE

                IREC = IREC + 1

                IF( IREC .GT. NMAX ) THEN

                    WRITE( MESG, 94010 ) 
     &                   'Maximum dimension =', NMAX,
     &                   'was exceeded while reading '//
     &                   BUFFER( 1:TRIMLEN( BUFFER ) ) // ' file.'
                    CALL M3EXIT( 'GETFACICR', 0, 0, MESG, 2 )

                ELSEIF( IREC .GT. CMAX ) THEN

                    WRITE( MESG, 94010 ) 
     &                   'Code parameter CMAX (=', CMAX,
     &                   ') was exceeded while reading '// 
     &                   BUFFER( 1:TRIMLEN( BUFFER ) )  // ' file.'
                    CALL  M3MESG( MESG )
                    CALL  M3EXIT( 'GETFACICR', 0, 0, 
     &                     'Change CMAX in GETFACICR and recompile', 2 )

                ELSE

C..................... For character dummy read, need to use actual array
C..................... to get the length correct.  The first read will be
C..................... overwritten by the second.
                    IF( READTHRU .EQ. 1 ) THEN
                      READ( IDEV, *, END=33, IOSTAT=IOS )
     &                      INDEXA( IREC ), ( IPDUM, J1=1,NCOL1 ),
     &                           ( CTAB2( J2,IREC ), J2=1,NCOL2 ),
     &                                      ( RPDUM, J3=1,NCOL3 )

                      SORTER( IREC ) = IREC

                    ELSE

                      READ( IDEV, *, END=33, IOSTAT=IOS ) IDUM,
     &                    ( ITAB1( J1, MOVETO( IREC )), J1=1,NCOL1 ),
     &                    ( CTAB2( J2, MOVETO( IREC )), J2=1,NCOL2 ),
     &                    ( RTAB3( J3, MOVETO( IREC )), J3=1,NCOL3 )

                    ENDIF

                    IF( IOS .NE. 0 ) THEN

                        EFLAG = .TRUE.
                        WRITE( MESG, 94010 )
     &                     'I/O Error:', IOS, ', while reading ' //
     &                     BUFFER( 1:TRIMLEN( BUFFER ) ) //
     &                 ' file record ', IREC
                        CALL M3MESG( MESG )

                    ENDIF

                ENDIF

            GO TO 11

   33       IF( EFLAG )
     &      CALL M3EXIT( 'GETFACICR', 0, 0, 'Error in file format', 2 )

            IF( READTHRU .EQ. 1 ) THEN

                N = IREC - 1

                CALL SORTI1( N, SORTER, INDEXA )

                DO 41 K = 1, N
                    INDEX( K ) = INDEXA( SORTER( K ) )
                    MOVETO( SORTER( K ) ) = K
   41           CONTINUE

                REWIND( IDEV )

            ENDIF

   44   CONTINUE

        RETURN

94010   FORMAT( 10 (A, :, 1X, I5, :, 1X ) )

        END

        SUBROUTINE GETFACIR( NAME, DEFAULT, NMAX, N, 
     &                      INDEX, NCOL1, ITAB1, NCOL2, RTAB2 )

C***********************************************************************
C  program body starts at line 47
C
C  DESCRIPTION:
C      Program prompts, opens, and reads an ascii file with an array 
C      of integer indices, and two real tables of indeterminant size
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Prototype 12/95 by MRH
C
C***********************************************************************

        IMPLICIT NONE

        INCLUDE 'PARMS3.EXT'     !  I/O API parameters
        INCLUDE 'IODECL3.EXT'    !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'     !  I/O API file description data structures

C.........  External
        INTEGER  GETEFILE
        INTEGER  TRIMLEN
        EXTERNAL GETEFILE, TRIMLEN

C.........  Program arguments
        INTEGER N, NCOL1, NCOL2, NMAX
        INTEGER INDEX( NMAX )

        INTEGER ITAB1( NCOL1, NMAX )
        REAL    RTAB2( NCOL2, NMAX )
 
        CHARACTER*(*) NAME, DEFAULT

C.........  Parameters
        INTEGER  CMAX
        PARAMETER  ( CMAX = 200 )   ! Constant maximum thanks to FORTRAN

C.........  Local variables
        INTEGER  INDEXA( CMAX )     ! Unsorted profile ID codes
        INTEGER  MOVETO( CMAX )     ! Index for rearranging profiles
        INTEGER  SORTER( CMAX )     ! Index for sorting profile ID codes

        INTEGER IDEV, IREC, J1, J2, K   ! Counters and pointers
        INTEGER IDUM                ! Dummy integer code
        INTEGER IPDUM               ! Dummy integer profile
        INTEGER IOS                 ! Read status
        INTEGER READTHRU            ! Loop counter for reading through file
        REAL    RPDUM               ! Dummy real profile
        LOGICAL EFLAG               ! Error flag
        DATA    EFLAG  / .FALSE. /

        CHARACTER*256 BUFFER
        CHARACTER*256 MESG

C***********************************************************************
C   begin body of program GETFACIR

C......... Set buffer for writing messages
        BUFFER = NAME

C.........  Open input file
        IDEV = GETEFILE( DEFAULT, .TRUE., .TRUE., 'GETFACIR' )
 
        IF ( IDEV .LT. 0 ) THEN     !  failure to open
 
            MESG = 'Could not open input file "' //
     &             DEFAULT( 1 : TRIMLEN( DEFAULT ) ) // '".'
            CALL M3EXIT( 'GETFACIR', 0, 0, MESG, 2 )
 
        ELSE
 
            CALL M3MSG2( 'NOTE: ' // BUFFER ( 1:TRIMLEN( BUFFER  ) ) // 
     &                   ' file "' // DEFAULT( 1:TRIMLEN( DEFAULT ) ) //
     &                   '" opened successfully.' )
 
        END IF      !  if getefile() failed

C.........  First time through file, establish sorted order, second time
C.........  through file, store profiles in that sorted order
        DO 44 READTHRU = 1, 2

            IREC = 0
   11       CONTINUE

                IREC = IREC + 1

                IF( IREC .GT. NMAX ) THEN

                    WRITE( MESG, 94010 ) 
     &                   'Maximum dimension =', NMAX,
     &                   'was exceeded while reading '//
     &                   BUFFER( 1:TRIMLEN( BUFFER ) ) // ' file.'
                    CALL  M3EXIT( 'GETFACIR', 0, 0, MESG, 2 )

                ELSEIF( IREC .GT. CMAX ) THEN

                    WRITE( MESG, 94010 ) 
     &                   'Code parameter CMAX (=', CMAX,
     &                   ') was exceeded while reading '// 
     &                   BUFFER( 1:TRIMLEN( BUFFER ) )  // ' file.'
                    CALL  M3MESG( MESG )
                    CALL  M3EXIT( 'GETFACIR', 0, 0, 
     &                      'Change CMAX in GETFACIR and recompile', 2 )

                ELSE

                    IF( READTHRU .EQ. 1 ) THEN
                      READ( IDEV, *, END=33, IOSTAT=IOS ) 
     &                    INDEXA( IREC ), ( IPDUM, J1=1,NCOL1 ),
     &                                    ( RPDUM, J2=1,NCOL2 )

                      SORTER( IREC ) = IREC

                    ELSE
                      READ( IDEV, *, END=33, IOSTAT=IOS ) IDUM,
     &                    ( ITAB1( J1, MOVETO( IREC )), J1=1,NCOL1 ),
     &                    ( RTAB2( J2, MOVETO( IREC )), J2=1,NCOL2 ) 

                    ENDIF

                    IF( IOS .NE. 0 ) THEN

                        EFLAG = .TRUE.
                        WRITE( MESG, 94010 )
     &                     'I/O Error:', IOS, ', while reading ' //
     &                     BUFFER( 1:TRIMLEN( BUFFER ) ) //
     &                     ' file record ', IREC
                        CALL M3MESG( MESG )

                    ENDIF

                ENDIF

            GO TO 11

   33       IF( EFLAG ) 
     &      CALL M3EXIT( 'GETFACIR', 0, 0, 'Error in file format', 2 )

            IF( READTHRU .EQ. 1 ) THEN

                N = IREC - 1

                CALL SORTI1( N, SORTER, INDEXA )

                DO 41 K = 1, N
                    INDEX( K ) = INDEXA( SORTER( K ) )
                    MOVETO( SORTER( K ) ) = K
   41           CONTINUE

                REWIND( IDEV )

            ENDIF

   44   CONTINUE

        RETURN

94010   FORMAT( 10 (A, :, 1X, I5, :, 1X ) )

        END

        SUBROUTINE GETIIR3D( NAME, DEFAULT, NMAX, N, 
     &                       INDEX, NCOL1, ITAB1, N2D, NCOL2,
     &                       ITAB2, NCOL3, RTAB3 )

C***********************************************************************
C  program body starts at line 47
C
C  DESCRIPTION:
C      Program prompts, opens, and reads an ascii file with an array 
C      of integer indices, an integer 3D table and a real 3D table of 
C      indeterminant sizes
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Prototype 12/95 by MRH
C
C***********************************************************************

        IMPLICIT NONE

        INCLUDE 'PARMS3.EXT'     !  I/O API parameters
        INCLUDE 'IODECL3.EXT'    !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'     !  I/O API file description data structures

C.........  External
        INTEGER  GETEFILE
        INTEGER  TRIMLEN
        EXTERNAL GETEFILE, TRIMLEN

C.........  Program arguments
        INTEGER N, N2D, NCOL1, NCOL2, NCOL3, NMAX
        INTEGER INDEX( NMAX )
        INTEGER ITAB1( NCOL1, NMAX )
        INTEGER ITAB2( NCOL2, N2D, NMAX )
        REAL    RTAB3( NCOL3, N2D, NMAX )
 
        CHARACTER*(*) NAME, DEFAULT

C.........  Parameters
        INTEGER  CMAX
        PARAMETER  ( CMAX = 200 )   ! Constant maximum thanks to FORTRAN

C.........  Local variables
        INTEGER  INDEXA( CMAX )     ! Unsorted profile ID codes
        INTEGER  MOVETO( CMAX )     ! Index for rearranging profiles
        INTEGER  SORTER( CMAX )     ! Index for sorting profile ID codes

        INTEGER I, IDEV, IREC, J2, J3, K   ! Counters and pointers
        INTEGER IDUM                ! Dummy integer code
        INTEGER IPDUM2              ! Dummy integer profile
        INTEGER IPDUM               ! Dummy integer profile
        INTEGER IOS                 ! Read status
        INTEGER READTHRU            ! Loop counter for reading through file
        REAL    RPDUM               ! Dummy real profile
        LOGICAL EFLAG               ! Error flag
        DATA    EFLAG  / .FALSE. /

        CHARACTER*256 BUFFER
        CHARACTER*256 MESG

C***********************************************************************
C   begin body of program GETIIR3D

C......... Set buffer for writing messages
        BUFFER = NAME

C.........  Open input file
        IDEV = GETEFILE( DEFAULT, .TRUE., .TRUE., 'GETIIR3D' )
 
        IF ( IDEV .LT. 0 ) THEN     !  failure to open
 
            MESG = 'Could not open input file "' //
     &             DEFAULT( 1 : TRIMLEN( DEFAULT ) ) // '".'
            CALL M3EXIT( 'GETIIR3D', 0, 0, MESG, 2 )
 
        ELSE
 
            CALL M3MSG2( 'NOTE: ' // BUFFER ( 1:TRIMLEN( BUFFER  ) ) // 
     &                   ' file "' // DEFAULT( 1:TRIMLEN( DEFAULT ) ) //
     &                   '" opened successfully.' )
 
        END IF      !  if getefile() failed

C.........  First time through file, establish sorted order, second time
C.........  through file, store profiles in that sorted order
        DO 44 READTHRU = 1, 2

            IREC = 0
   11       CONTINUE

                IREC = IREC + 1

                IF( IREC .GT. NMAX ) THEN

                    WRITE( MESG, 94010 ) 
     &                   'Maximum dimension =', NMAX,
     &                   'was exceeded while reading '//
     &                   BUFFER( 1:TRIMLEN( BUFFER ) ) // ' file.'
                    CALL  M3EXIT( 'GETIIR3D', 0, 0, MESG, 2 )

                ELSEIF( IREC .GT. CMAX ) THEN

                    WRITE( MESG, 94010 ) 
     &                   'Code parameter CMAX (=', CMAX,
     &                   ') was exceeded while reading '// 
     &                   BUFFER( 1:TRIMLEN( BUFFER ) )  // ' file.'
                    CALL  M3MESG( MESG )
                    CALL  M3EXIT( 'GETIIR3D', 0, 0, 
     &                      'Change CMAX in GETIIR3D and recompile', 2 )

                ELSE

                    IF( READTHRU .EQ. 1 ) THEN
                      READ( IDEV, *, END=33, IOSTAT=IOS ) 
     &                    INDEXA( IREC ), IPDUM

                      IF( IPDUM .GT. N2D ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG, 94010 )
     &                     'Dimension 2: ', N2D,
     &                     'too small while reading ' //
     &                     BUFFER( 1:TRIMLEN( BUFFER ) ) //
     &                     ' file record ', IREC
                        CALL M3MESG( MESG )

                      ENDIF

                      DO 22 I = 1, IPDUM
                          READ( IDEV, *, END=29, IOSTAT=IOS )
     &                                  ( IPDUM2, J2=1,NCOL2 ),
     &                                  ( RPDUM , J3=1,NCOL3 )
   22                 CONTINUE

                      SORTER( IREC ) = IREC

                    ELSE

                      READ( IDEV, *, END=33, IOSTAT=IOS )
     &                    IDUM, ITAB1( 1, MOVETO( IREC ) )

                      DO 23 I = 1, ITAB1( 1, MOVETO( IREC ) )

                        READ( IDEV, *, END=29, IOSTAT=IOS )
     &                      ( ITAB2( J2,I,MOVETO(IREC) ), J2=1,NCOL2 ),
     &                      ( RTAB3( J3,I,MOVETO(IREC) ), J3=1,NCOL3 )

   23                 CONTINUE

                    ENDIF

                    IF( IOS .NE. 0 ) THEN

                        EFLAG = .TRUE.
                        WRITE( MESG, 94010 )
     &                     'I/O Error:', IOS, ', while reading ' //
     &                     BUFFER( 1:TRIMLEN( BUFFER ) ) //
     &                     ' file record ', IREC
                        CALL M3MESG( MESG )

                    ENDIF

                ENDIF

            GO TO 11

C.............  Report error if file ends prematurely
   29       CONTINUE
            WRITE( MESG,94010 ) BUFFER( 1:TRIMLEN( BUFFER ) ) // 
     &                          'file ended prematurely at record', IREC
            CALL M3EXIT( 'GETIIR3D', 0, 0, 'Error in file format', 2 )

   33       IF( EFLAG ) 
     &      CALL M3EXIT( 'GETIIR3D', 0, 0, 'Error in file format', 2 )

            IF( READTHRU .EQ. 1 ) THEN

                N = IREC - 1

                CALL SORTI1( N, SORTER, INDEXA )

                DO 41 K = 1, N
                    INDEX( K ) = INDEXA( SORTER( K ) )
                    MOVETO( SORTER( K ) ) = K
   41           CONTINUE

                REWIND( IDEV )

            ENDIF

   44   CONTINUE

        RETURN

94010   FORMAT( 10 (A, :, 1X, I5, :, 1X ) )

        END

        SUBROUTINE GETFACD( NAME, DEFAULT, NMAX, N, 
     &                      INDEX, ISIZE, NCOL, RTAB1 )

C***********************************************************************
C  program body starts at line
C
C  DESCRIPTION:
C      Program prompts, opens, and reads an ascii file with an array 
C      of integer indices, an array of actual counts of real values, and 
C      a real table of indeterminant number of rows and columns
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Prototype 12/95 by MRH
C
C***********************************************************************

        IMPLICIT NONE

        INCLUDE 'PARMS3.EXT'     !  I/O API parameters
        INCLUDE 'IODECL3.EXT'    !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'     !  I/O API file description data structures

C.........  External
        INTEGER  GETEFILE
        INTEGER  TRIMLEN
        EXTERNAL GETEFILE, TRIMLEN

C.........  Program arguments
        INTEGER N, NCOL, NMAX
        INTEGER INDEX( NMAX )       ! Sorted profile ID codes
        INTEGER ISIZE( NMAX )       ! Actual number of RTAB1 columns

        REAL    RTAB1( NCOL, NMAX ) ! Sorted profile of real numbers
 
        CHARACTER*(*) NAME, DEFAULT

C.........  Parameters
        INTEGER  CMAX
        PARAMETER  ( CMAX = 200 )   ! Constant maximum thanks to FORTRAN

C.........  Local variables
        INTEGER  INDEXA( CMAX )     ! Unsorted profile ID codes
        INTEGER  MOVETO( CMAX )     ! Index for rearranging profiles
        INTEGER  SORTER( CMAX )     ! Index for sorting profile ID codes

        INTEGER IDEV, IREC, J1, K   ! Counters and pointers
        INTEGER IDUM                ! Dummy integer code
        INTEGER IOS                 ! Read status
        INTEGER NTMP                ! Actual number of table entries
        INTEGER READTHRU            ! Loop counter for reading through file
        REAL    RPDUM               ! Dummy real profile
        LOGICAL EFLAG               ! Error flag
        DATA    EFLAG  / .FALSE. /

        CHARACTER*256 BUFFER
        CHARACTER*256 MESG

C***********************************************************************
C   begin body of program GETFACD

C......... Set buffer for writing messages
        BUFFER = NAME

C.........  Open input file
        IDEV = GETEFILE( DEFAULT, .TRUE., .TRUE., 'GETFACD' )
 
        IF ( IDEV .LT. 0 ) THEN     !  failure to open
 
            MESG = 'Could not open input file "' //
     &             DEFAULT( 1 : TRIMLEN( DEFAULT ) ) // '".'
            CALL M3EXIT( 'GETFACD', 0, 0, MESG, 2 )
 
        ELSE
 
            CALL M3MSG2( 'NOTE: ' // BUFFER ( 1:TRIMLEN( BUFFER  ) ) // 
     &                   ' file "' // DEFAULT( 1:TRIMLEN( DEFAULT ) ) //
     &                   '" opened successfully.' )
 
        END IF      !  if getefile() failed

C.........  First time through file, establish sorted order, second time
C.........  through file, store profiles in that sorted order
        DO 44 READTHRU = 1, 2

            IREC = 0
   11       CONTINUE

                IREC = IREC + 1

                IF( IREC .GT. NMAX ) THEN

                    WRITE( MESG, 94010 ) 
     &                   'Maximum dimension =', NMAX,
     &                   'was exceeded while reading '// 
     &                   BUFFER( 1:TRIMLEN( BUFFER ) )  // ' file.'
                    CALL  M3EXIT( 'GETFACD', 0, 0, MESG, 2 )

                ELSEIF( IREC .GT. CMAX ) THEN

                    WRITE( MESG, 94010 ) 
     &                   'Code parameter CMAX (=', CMAX,
     &                   ') was exceeded while reading '// 
     &                   BUFFER( 1:TRIMLEN( BUFFER ) )  // ' file.'
                    CALL  M3MESG( MESG )
                    CALL  M3EXIT( 'GETFACD', 0, 0, 
     &                      'Change CMAX in GETFACD and recompile', 2 )

                ELSE

                    IF( READTHRU .EQ. 1 ) THEN
                      READ( IDEV, *, END=33, IOSTAT=IOS ) 
     &                      INDEXA( IREC ), NTMP, ( RPDUM, J1=1,NTMP )

                      SORTER( IREC ) = IREC

C.......................  Initialize parameter array to zero
                      DO 22 K = 1, NCOL
                          RTAB1( K, IREC ) = 0.E0
   22                 CONTINUE

                    ELSE
                      READ( IDEV, *, END=33, IOSTAT=IOS ) IDUM, NTMP,
     &                    ( RTAB1( J1, MOVETO( IREC )), J1=1,NTMP )

                      ISIZE( MOVETO( IREC ) ) = NTMP

                    ENDIF
    
                    IF( IOS .NE. 0 ) THEN

                        EFLAG = .TRUE.
                        WRITE( MESG, 94010 )
     &                     'I/O Error:', IOS, ', while reading ' //
     &                     BUFFER( 1:TRIMLEN( BUFFER ) ) // 
     &                     ' file record ', IREC
                        CALL M3MESG( MESG )

                    ENDIF

                ENDIF

            GO TO 11

   33       IF( EFLAG ) 
     &      CALL M3EXIT( 'GETFACD', 0, 0, 'Error in file format', 2 )

            IF( READTHRU .EQ. 1 ) THEN

                N = IREC - 1

                CALL SORTI1( N, SORTER, INDEXA )

                DO 41 K = 1, N
                    INDEX( K ) = INDEXA( SORTER( K ) )
                    MOVETO( SORTER( K ) ) = K
   41           CONTINUE

                REWIND( IDEV )

            ENDIF

   44   CONTINUE

        RETURN

94010   FORMAT( 10 (A, :, 1X, I5, :, 1X ) )

        END


