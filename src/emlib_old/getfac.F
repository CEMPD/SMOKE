
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
        INTEGER  PROMPTFFILE
        INTEGER  TRIMLEN
        EXTERNAL PROMPTFFILE, TRIMLEN

C.........  Program arguments
        INTEGER N, NCOL, NMAX
        INTEGER INDEX( NMAX )

        REAL    RTAB1( NCOL, NMAX )
 
        CHARACTER*(*) NAME, DEFAULT

C.........  Local variables
        INTEGER IDEV, IOS, IREC, J1 
        LOGICAL EFLAG
        DATA    EFLAG  / .FALSE. /

        CHARACTER*256 BUFFER
        CHARACTER*256 MESG

C***********************************************************************
C   begin body of program GETFACR

C......... Set buffer for writing messages
        BUFFER = NAME

C.........  Prompt for input file

        IDEV = PROMPTFFILE(
     &         'Enter logical name for '// NAME( 1:TRIMLEN( NAME ) )
     &         //' input file', .TRUE., .TRUE., DEFAULT, 'GETFACR' )

        IREC = 0
   11   CONTINUE

            IREC = IREC + 1

            IF( IREC .GT. NMAX ) THEN

                WRITE( MESG, 94010 ) 
     &               'Maximum dimension =', NMAX,
     &               'was exceeded while reading '// 
     &               BUFFER( 1:TRIMLEN( BUFFER ) )  // ' file.'
                CALL  M3EXIT( 'GETFACR', 0, 0, MESG, 2 )

            ELSE

                READ( IDEV, *, END = 33, IOSTAT=IOS ) 
     &                INDEX( IREC ), ( RTAB1( J1, IREC ), J1=1,NCOL )

                IF( IOS .NE. 0 ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 )
     &                 'I/O Error:', IOS, ', while reading ' //
     &                 BUFFER( 1:TRIMLEN( BUFFER ) ) // 
     &                 ' file record ', IREC
                    CALL M3MESG( MESG )

                ENDIF

            ENDIF

        GO TO 11

   33   IF( EFLAG ) THEN

            CALL M3EXIT( 'GETFACR', 0, 0, 'Error in file format', 2 )

        ENDIF

        N = IREC - 1

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
        INTEGER  PROMPTFFILE
        INTEGER  TRIMLEN
        EXTERNAL PROMPTFFILE, TRIMLEN

C.........  Program arguments
        INTEGER N, NCOL1, NCOL2, NMAX
        INTEGER INDEX( NMAX )

        REAL    RTAB1( NCOL1, NMAX )
        REAL    RTAB2( NCOL2, NMAX )
 
        CHARACTER*(*) NAME, DEFAULT

C.........  Local variables
        INTEGER IDEV, IOS, IREC, J1, J2
        LOGICAL EFLAG
        DATA    EFLAG  / .FALSE. /

        CHARACTER*256 BUFFER
        CHARACTER*256 MESG

C***********************************************************************
C   begin body of program GETFACRR

C......... Set buffer for writing messages
        BUFFER = NAME

C.........  Prompt for input file

        IDEV = PROMPTFFILE(
     &         'Enter logical name for '// NAME( 1:TRIMLEN( NAME ) )
     &         //' input file', .TRUE., .TRUE., DEFAULT, 'GETFACRR' )

        IREC = 0
   11   CONTINUE

            IREC = IREC + 1

            IF( IREC .GT. NMAX ) THEN

                WRITE( MESG, 94010 ) 
     &               'Maximum dimension =', NMAX,
     &               'was exceeded while reading '// 
     &               BUFFER( 1:TRIMLEN( BUFFER ) ) // ' file.'
                CALL  M3EXIT( 'GETFACRR', 0, 0, MESG, 2 )

            ELSE

                READ( IDEV, *, END = 33, IOSTAT=IOS ) INDEX( IREC ),
     &                ( RTAB1( J1, IREC ), J1=1,NCOL1 ),
     &                ( RTAB2( J2, IREC ), J2=1,NCOL2 ) 

                IF( IOS .NE. 0 ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 )
     &                 'I/O Error:', IOS, ', while reading ' //
     &                 BUFFER( 1:TRIMLEN( BUFFER ) ) //
     &                 ' file record ', IREC
                    CALL M3MESG( MESG )

                ENDIF

            ENDIF

        GO TO 11

   33   IF( EFLAG ) THEN

            CALL M3EXIT( 'GETFACRR', 0, 0, 'Error in file format', 2 )

        ENDIF

        N = IREC - 1

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
        INTEGER  PROMPTFFILE
        INTEGER  TRIMLEN
        EXTERNAL PROMPTFFILE, TRIMLEN

C.........  Program arguments
        INTEGER       N, NCOL1, NMAX

        INTEGER       INDEX( NMAX )
        INTEGER       ITAB1( NCOL1, NMAX )

        CHARACTER*(*) NAME, DEFAULT

C.........  Local variables
        INTEGER IDEV, IOS, IREC, J1
        LOGICAL EFLAG
        DATA    EFLAG  / .FALSE. /

        CHARACTER*256 BUFFER
        CHARACTER*256 MESG

C***********************************************************************
C   begin body of program GETFACI

C......... Set buffer for writing messages
        BUFFER = NAME

C.........  Prompt for input file

        IDEV = PROMPTFFILE(
     &         'Enter logical name for '// NAME( 1:TRIMLEN( NAME ) )
     &         //' input file', .TRUE., .TRUE., DEFAULT, 'GETFACI' )

        IREC = 0
   11   CONTINUE

            IREC = IREC + 1

            IF( IREC .GT. NMAX ) THEN

                WRITE( MESG, 94010 ) 
     &               'Maximum dimension =', NMAX,
     &               'was exceeded while reading '// 
     &               BUFFER( 1:TRIMLEN( BUFFER ) ) // ' file.'
                CALL M3EXIT( 'GETFACI', 0, 0, MESG, 2 )

            ELSE

                READ( IDEV, *, END = 33, IOSTAT=IOS ) 
     &              INDEX( IREC ),
     &              ( ITAB1( J1, IREC ), J1=1,NCOL1 )

                IF( IOS .NE. 0 ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 )
     &                 'I/O Error:', IOS, ', while reading ' //
     &                 BUFFER( 1:TRIMLEN( BUFFER ) ) //
     &                 ' file record ', IREC
                    CALL M3MESG( MESG )

                ENDIF

            ENDIF

        GO TO 11

   33   IF( EFLAG ) THEN

            CALL M3EXIT( 'GETFACI', 0, 0, 'Error in file format', 2 )

        ENDIF

        N = IREC - 1

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
        INTEGER  PROMPTFFILE
        INTEGER  TRIMLEN
        EXTERNAL PROMPTFFILE, TRIMLEN

C.........  Program arguments
        INTEGER       N, NCOL1, NCOL2, NCOL3, NMAX

        INTEGER       INDEX( NMAX )
        INTEGER       ITAB1( NCOL1, NMAX )

        REAL          RTAB3( NCOL3, NMAX )
 
        CHARACTER*(*) CTAB2( NCOL2, NMAX )
        CHARACTER*(*) NAME, DEFAULT

C.........  Local variables
        INTEGER IDEV, IOS, IREC, J1, J2, J3
        LOGICAL EFLAG
        DATA    EFLAG  / .FALSE. /

        CHARACTER*256 BUFFER
        CHARACTER*256 MESG

C***********************************************************************
C   begin body of program GETFACICR

C......... Set buffer for writing messages
        BUFFER = NAME

C.........  Prompt for input file

        IDEV = PROMPTFFILE(
     &         'Enter logical name for '// NAME( 1:TRIMLEN( NAME ) )
     &         //' input file', .TRUE., .TRUE., DEFAULT, 'GETFACICR' )

        IREC = 0
   11   CONTINUE

            IREC = IREC + 1

            IF( IREC .GT. NMAX ) THEN

                WRITE( MESG, 94010 ) 
     &               'Maximum dimension =', NMAX,
     &               'was exceeded while reading '//
     &               BUFFER( 1:TRIMLEN( BUFFER ) ) // ' file.'
                CALL M3EXIT( 'GETFACICR', 0, 0, MESG, 2 )

            ELSE

                READ( IDEV, *, END = 33, IOSTAT=IOS ) 
     &              INDEX( IREC ),
     &              ( ITAB1( J1, IREC ), J1=1,NCOL1 ),
     &              ( CTAB2( J2, IREC ), J2=1,NCOL2 ),
     &              ( RTAB3( J3, IREC ), J3=1,NCOL3 )

                IF( IOS .NE. 0 ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 )
     &                 'I/O Error:', IOS, ', while reading ' //
     &                 BUFFER( 1:TRIMLEN( BUFFER ) ) //
     &                 ' file record ', IREC
                    CALL M3MESG( MESG )

                ENDIF

            ENDIF

        GO TO 11

   33   IF( EFLAG ) THEN

            CALL M3EXIT( 'GETFACICR', 0, 0, 'Error in file format', 2 )

        ENDIF

        N = IREC - 1

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
        INTEGER  PROMPTFFILE
        INTEGER  TRIMLEN
        EXTERNAL PROMPTFFILE, TRIMLEN

C.........  Program arguments
        INTEGER N, NCOL1, NCOL2, NMAX
        INTEGER INDEX( NMAX )

        INTEGER ITAB1( NCOL1, NMAX )
        REAL    RTAB2( NCOL2, NMAX )
 
        CHARACTER*(*) NAME, DEFAULT

C.........  Local variables
        INTEGER IDEV, IOS, IREC, J1, J2
        LOGICAL EFLAG
        DATA    EFLAG  / .FALSE. /

        CHARACTER*256 BUFFER
        CHARACTER*256 MESG

C***********************************************************************
C   begin body of program GETFACIR

C......... Set buffer for writing messages
        BUFFER = NAME

C.........  Prompt for input file

        IDEV = PROMPTFFILE(
     &         'Enter logical name for '// NAME( 1:TRIMLEN( NAME ) )
     &         //' input file', .TRUE., .TRUE., DEFAULT, 'GETFACIR' )

        IREC = 0
   11   CONTINUE

            IREC = IREC + 1

            IF( IREC .GT. NMAX ) THEN

                WRITE( MESG, 94010 ) 
     &               'Maximum dimension =', NMAX,
     &               'was exceeded while reading '//
     &               BUFFER( 1:TRIMLEN( BUFFER ) ) // ' file.'
                CALL  M3EXIT( 'GETFACIR', 0, 0, MESG, 2 )

            ELSE

                READ( IDEV, *, END = 33, IOSTAT=IOS ) INDEX( IREC ),
     &                ( ITAB1( J1, IREC ), J1=1,NCOL1 ),
     &                ( RTAB2( J2, IREC ), J2=1,NCOL2 ) 

                IF( IOS .NE. 0 ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 )
     &                 'I/O Error:', IOS, ', while reading ' //
     &                 BUFFER( 1:TRIMLEN( BUFFER ) ) //
     &                 ' file record ', IREC
                    CALL M3MESG( MESG )

                ENDIF

            ENDIF

        GO TO 11

   33   IF( EFLAG ) THEN

            CALL M3EXIT( 'GETFACR', 0, 0, 'Error in file format', 2 )

        ENDIF

        N = IREC - 1

        RETURN

94010   FORMAT( 10 (A, :, 1X, I5, :, 1X ) )

        END
