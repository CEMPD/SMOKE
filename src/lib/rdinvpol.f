
        SUBROUTINE RDINVPOL( FILNAM, NCNT, VCNT, VNAMES, VTYPES, SRCID, 
     &                       POLDAT, EFLAG )

C***********************************************************************
C  subroutine body starts at line 82
C
C  DESCRIPTION:
C      Reads inventory pollutant-specific data for variables listed in VNAMES
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created ??/???? by ???
C       09/2025 by HT UNC-IE:  Use M3UTILIO
C
C**************************************************************************
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

C.........  MODULES for public variables

C.........  This module is required by the FileSetAPI
        USE MODFILESET

        IMPLICIT NONE

C...........   INCLUDE FILES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables
c       INCLUDE 'IODECL3.EXT'   !  I/O API function declarations

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
c       CHARACTER(2) CRLF
c       EXTERNAL     CRLF

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(IOFLEN3), INTENT (IN) :: FILNAM     ! Logical file name
        INTEGER     , INTENT (IN) :: NCNT             ! Number of records
        INTEGER     , INTENT (IN) :: VCNT             ! No. vars other than SRCID
        CHARACTER(IOVLEN3), INTENT (IN) :: VNAMES( VCNT )   ! Variable names
        INTEGER     , INTENT (IN) :: VTYPES( VCNT )   ! Variable types
        INTEGER     , INTENT(OUT) :: SRCID ( NCNT )   ! Data
        REAL        , INTENT(OUT) :: POLDAT( NCNT,VCNT ) ! Data
        LOGICAL     , INTENT(OUT) :: EFLAG            ! true: error found

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: IREAD ( : )  ! integer read array

C...........   Other local variables

        INTEGER         K, LV, V  ! counters and indices
        INTEGER         IOS       ! i/o status

        CHARACTER(IOVLEN3)   VARBUF
        CHARACTER(256)       MESG 

        CHARACTER(16) :: PROGNAME = 'RDINVPOL' ! program name

C***********************************************************************
C   begin body of subroutine RDINVPOL

        EFLAG = .FALSE.

        IF( .NOT. READSET( FILNAM, 'SRCID', ALLAYS3, 1,
     &                     0, 0, SRCID                    ) ) THEN
            EFLAG = .TRUE.
            MESG = 'Error reading "SRCID" from file: ' // CRLF()//
     &             BLANK10// TRIM( FILNAM )
            CALL M3MSG2( MESG )

        END IF

C.........  Read variables by type, and store as REAL in POLDAT
        DO V = 1, VCNT

            VARBUF = VNAMES( V )
            LV = LEN_TRIM( VARBUF )

            IF( VTYPES( V ) .EQ. M3INT ) THEN

C.................  If memory is not allocated for integer read array, then
C                   allocate it
                IF( .NOT. ALLOCATED( IREAD ) ) THEN
                    ALLOCATE( IREAD( NCNT ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'IREAD', PROGNAME )
                END IF

                IF( .NOT. READSET( FILNAM, VARBUF, ALLAYS3,
     &                             ALLFILES, 0, 0, IREAD )) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not read "' //
     &                      VARBUF( 1:LV ) // '" from file.'
                    CALL M3MSG2( MESG )

                ELSE

                    POLDAT( :,V ) = REAL( IREAD )   ! Array

                END IF

            ELSE IF( .NOT. READSET( FILNAM, VARBUF, ALLAYS3,
     &                              ALLFILES, 0, 0, POLDAT(1,V) ) ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not read "' //
     &                 VARBUF( 1:LV ) // '" from file.'
                CALL M3MSG2( MESG )

            END IF

        END DO

C.........  Deallocate local memory
        IF( ALLOCATED( IREAD ) ) DEALLOCATE( IREAD )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

        END

