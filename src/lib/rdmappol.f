
        SUBROUTINE RDMAPPOL( NSRC, NVARS, NPVAR, VARNAMS, OUTVALS )

C***********************************************************************
C  program body starts at line 
C
C  DESCRIPTION:
C     Reads in pollutant or activity data and associated 
C     variables from map-formatted or old-format inventory files.  
C     Stores the data by source instead of by sparse record.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Version ??/???? by ???
C       09/2025 by HT UNC-IE:  Use M3UTILIO
C
C************************************************************************
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
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NMAP, MAPNAM, MAPFIL

        USE MODFILESET

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
c       INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
c       CHARACTER(2)     CRLF
c       INTEGER          INDEX1
c       LOGICAL          SETENVVAR

c       EXTERNAL    CRLF, INDEX1, SETENVVAR

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NSRC                 ! no. inven sources        
        INTEGER     , INTENT (IN) :: NVARS                ! number of pol/act
        INTEGER     , INTENT (IN) :: NPVAR                ! number of vars per pol/act
        CHARACTER(*), INTENT (IN) :: VARNAMS( NVARS )     ! vars to map
        REAL        , INTENT(OUT) :: OUTVALS( NSRC, NVARS*NPVAR )! output data

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: SRCID( :   )  ! source IDs in sparse pol files
        REAL   , ALLOCATABLE :: ETMP ( :,: )  ! emissions in sparse pol files

C...........   Other local variables
        INTEGER          I, J, L, M, N, S, V   ! counters and indices

        INTEGER          ADJINDX      ! adjustment index for special average day read
        INTEGER          IOS          ! i/o status
        INTEGER          NSPARSE      ! number of rows in sparse pol input file
        INTEGER          NSRCLOC      ! number of rows in sparse pol input file

        LOGICAL       :: EFLAG = .FALSE.  ! true: error found
        LOGICAL       :: OFLAG = .FALSE.  ! true: old inventory format
        LOGICAL       :: SFLAG = .FALSE.  ! true: error found on read of inven data

        CHARACTER(16)   :: RNAME = 'IOAPI_DAT' ! logical name for reading pols
        CHARACTER(256)     MESG   ! message buffer
        CHARACTER(IOVLEN3) VBUF   ! tmp variable name buffer

        CHARACTER(16) :: PROGNAME = 'RDMAPPOL'   !  program name

C***********************************************************************
C   begin body of subroutine RDMAPPOL

C........  Initialize output arrays to zero
        OUTVALS = 0.   ! array

        DO V = 1, NVARS

            ADJINDX = 2
            VBUF = VARNAMS( V )

C............  Check if average day prefix is in the name
            IF( VARNAMS( V )( 1:3 ) .EQ. AVEDAYRT ) THEN
                L = LEN_TRIM( VARNAMS( V ) )
                VBUF = VARNAMS( V )( 4:L )  
                ADJINDX = 3

C...............  Give error if reading more than a single variable
                IF( NPVAR .GT. 1 ) THEN
                    MESG = 'INTERNAL ERROR: Cannot specify average ' //
     &                     'day value and have NPVAR > 1'
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

            END IF

C............  Find variable in map
            M = INDEX1( VBUF, NMAP, MAPNAM )
            IF( M .LE. 0 ) THEN
                MESG = 'INTERNAL ERROR: Could not find variable "'//
     &                 TRIM( VBUF ) // '" in map file'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C............  Open physical file name for this pollutant or activity
C............  Also, get description
            CALL OPENPHYS( PROGNAME, RNAME, FSREAD3, MAPFIL( M ),
     &                     EFLAG )
            
            NSPARSE = NROWS3D

C............  Allocate memory for local read arrays
            IF( ALLOCATED( ETMP ) ) DEALLOCATE( ETMP, SRCID )
            ALLOCATE( ETMP( NSPARSE, NPVAR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ETMP', PROGNAME )
            ALLOCATE( SRCID( NSPARSE ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SRCID', PROGNAME )

C.............  Use subroutine to read in data. This subroutine
C               already handles integer and real
            J = ADJINDX
            CALL RDINVPOL( RNAME, NSPARSE, NPVAR, VNAMESET( J ), 
     &                     VTYPESET( J ), SRCID, ETMP, SFLAG     )

            IF( SFLAG ) EFLAG = .TRUE.
            IF( EFLAG ) CYCLE

C...........  Transfer sparse-storage pollutant to output arrays
            DO I = 1, NPVAR

                J = I + ( V-1 ) * NPVAR
                DO N = 1, NSPARSE
                    S = SRCID( N )
                    OUTVALS( S,J ) = ETMP( N,I )
                END DO
 
            END DO     ! End loop over variable per pol or act (if any)

C...............  Close output file for this variable
            IF( .NOT. CLOSESET( RNAME ) ) THEN
                MESG = 'Could not close file:'//CRLF()//BLANK10//
     &                     TRIM( MAPFIL( M ) )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        END DO         ! End loop over pols or acts

C.........  Abort if an error has been found
        IF( EFLAG ) THEN

            MESG = 'Problem reading pollutant data'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

        IF( ALLOCATED( SRCID ) ) DEALLOCATE( SRCID, ETMP )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000       FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

94100   FORMAT( 1X,  A, 1X,  A, I7 )

        END SUBROUTINE RDMAPPOL


