
        SUBROUTINE WRINVPOL( CATEGORY, INPATH, INNAM, NREC, 
     &                       NPVAR, VBUF, SRCID, POLBUF, EFLAG )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine writes the inventory pollutant-based variables
C      and/or VMT to the I/O API files in sparse storage format
C
C  PRECONDITIONS REQUIRED:
C      Number of sources NSRC defined correctly
C      Pollutant count IPCNT and output names POLNAM defined correctly
C      Output array POLVAL populated
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutine, BLDENAMS
C      Functions: I/O API subroutine
C
C  REVISION  HISTORY:
C      Created 11/98 by M. Houyoux
C     09/2025 by HT UNC-IE:  Use M3UTILIO
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************
        USE M3UTILIO

C...........   MODULES for public variables
        USE MODFILESET

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
c       INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI function declarations

C.........  EXTERNAL FUNCTIONS
c       CHARACTER(2)  CRLF
c       LOGICAL       SETENVVAR

c       EXTERNAL      CRLF, SETENVVAR

C.............  Subroutine arguments 
        CHARACTER(*), INTENT (IN) :: CATEGORY        ! source category
        CHARACTER(*), INTENT (IN) :: INPATH          ! path for output file
        CHARACTER(*), INTENT (IN) :: INNAM           ! name for output file
        INTEGER     , INTENT (IN) :: NREC            ! size of output arrays
        INTEGER     , INTENT (IN) :: NPVAR           ! number variables per pol/act
        CHARACTER(*), INTENT (IN) :: VBUF            ! names of pols/act
        INTEGER     , INTENT (IN) :: SRCID( NREC )   ! source IDs
        REAL        , INTENT (IN) :: POLBUF( NREC,NPVAR )  ! pol or act data
        LOGICAL     , INTENT(OUT) :: EFLAG           ! true: error found

C.............  Arrays for sparse output
        INTEGER, ALLOCATABLE :: INTBUF ( : )

C.............  Local variables
        INTEGER      J, L, N
        INTEGER      IOS      ! i/o status

        LOGICAL       :: MSGFLAG  = .FALSE.  ! true: need to output error message

        CHARACTER(16)      :: FNAME = 'IOAPI_DAT' ! tmp logical name for outputs
        CHARACTER(256)     :: MESG                ! message buffer
        CHARACTER(PHYLEN3) :: FPHYS = ' '         ! full physical file name

        CHARACTER(16) :: PROGNAME = 'WRINVPOL' !  program name

C***********************************************************************
C   begin body of program WRINVPOL

C........  Allocate memory for output arrays
        ALLOCATE( INTBUF( NREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INTBUF', PROGNAME )

        FPHYS = TRIM( INPATH ) // TRIM( INNAM )

C.........  Set output logical file name
        IF( .NOT. SETENVVAR( FNAME, FPHYS ) ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: Could not set logical file name ' //
     &             'for "' // TRIM( VBUF ) // '" file:'// 
     &             CRLF()// BLANK10// TRIM( FPHYS )
            CALL M3MSG2( MESG )

C.........  Open I/O API file
        ELSE IF( .NOT. OPENSET( FNAME, FSNEW3, PROGNAME )) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: Could not open I/O API "' // TRIM( VBUF )//
     &             '" file for file name:' // CRLF() // BLANK10 // 
     &             TRIM( FPHYS )
            CALL M3MSG2( MESG )

        END IF

C.........  Write source ID index
        IF( .NOT. WRITESET( FNAME, 'SRCID', ALLFILES, 
     &                      0, 0, SRCID               )) THEN
            EFLAG = .TRUE.
            MESG = 'Could not write "SRCID" to file:'//
     &              CRLF()//BLANK10// TRIM( FPHYS )
            CALL M3MSG2( MESG )
        END IF

C.........  Write all variables for each pollutant or activity
        DO J = 1, NPVAR

            MSGFLAG = .FALSE.
            IF( VTYPESET( J+1 ) .EQ. M3REAL ) THEN
                IF( .NOT. WRITESET( FNAME, VNAMESET( J+1 ), 
     &              ALLFILES, 0, 0, POLBUF( 1,J )   ) ) MSGFLAG = .TRUE.

            ELSE
                INTBUF( 1:NREC ) = INT( POLBUF( 1:NREC, J ) )  ! array
                IF( .NOT. WRITESET( FNAME, VNAMESET( J+1 ), 
     &              ALLFILES, 0, 0, INTBUF          ) ) MSGFLAG = .TRUE.
            END IF

            IF( MSGFLAG ) THEN
                EFLAG = .TRUE.
                MESG = 'Could not write "'// TRIM( VNAMESET(J+1) )
     &                 // '" to file:'// CRLF()// BLANK10// 
     &                 TRIM( FPHYS )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 ) 
            END IF

        END DO

C........  Close output file for this variable
        IF( .NOT. CLOSESET( FNAME ) ) THEN
            MESG = 'Could not close file:'//CRLF()//BLANK10//
     &                     TRIM( FPHYS )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C...............  Deallocate memory for temporary write arrays
        DEALLOCATE( INTBUF )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93100   FORMAT( I2 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94100   FORMAT( 9( A, I2 ) )

        END SUBROUTINE WRINVPOL
