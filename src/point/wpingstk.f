
        SUBROUTINE WPINGSTK( FNAME, SDATE, STIME )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C       Opens the output files for the Elevpoint program
C       
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 8/99 by M Houyoux
C
C************************************************************************
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
C Last updated: %G 
C  
C***********************************************************************

C...........   MODULES for public variables
C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2     CRLF
        EXTERNAL        CRLF

C..........    Subroutine arguments and their descriptions
        CHARACTER(*), INTENT (IN) :: FNAME   ! i/o api inventory file
        INTEGER     , INTENT (IN) :: SDATE   ! Julian start date
        INTEGER     , INTENT (IN) :: STIME   ! start time

C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: SCCSW  = '@(#)$Id$'  ! SCCS string with vers no.

C...........   Local variables that are dimensioned by module variables
C...........   These are for sorting groups and outputting in sorted order
        INTEGER         LOCIDX( NGROUP )
        INTEGER         LOCGID( NGROUP )
        INTEGER         LOCCNT( NGROUP )
        INTEGER         LOCCOL( NGROUP )
        INTEGER         LOCROW( NGROUP )

        REAL            LOCDM ( NGROUP )
        REAL            LOCFL ( NGROUP )
        REAL            LOCHT ( NGROUP )
        REAL            LOCLAT( NGROUP )
        REAL            LOCLON( NGROUP )
        REAL            LOCTK ( NGROUP )
        REAL            LOCVE ( NGROUP )
        REAL            LOCXL ( NGROUP )
        REAL            LOCYL ( NGROUP )

C...........   Other local variables
        INTEGER         I, J, L      ! indices and counters
         
        LOGICAL      :: EFLAG    = .FALSE.  !  true: error found

        CHARACTER*300   MESG

        CHARACTER*16 :: PROGNAME = 'WPINGSTK'   !  subroutine name

C***********************************************************************
C   begin body of subroutine WPINGSTK

C.........  Store sorted information
        DO I = 1, NGROUP
            J = GRPIDX( I )

            LOCGID( I ) = GRPGIDA( J )
            LOCCNT( I ) = GRPCNT( J )
            LOCCOL( I ) = GRPCOL( J )
            LOCROW( I ) = GRPROW( J )
            LOCDM ( I ) = GRPDM ( J )
            LOCFL ( I ) = GRPFL ( J )
            LOCHT ( I ) = GRPHT ( J )
            LOCLAT( I ) = GRPLAT( J )
            LOCLON( I ) = GRPLON( J )
            LOCTK ( I ) = GRPTK ( J )
            LOCVE ( I ) = GRPVE ( J )
            LOCXL ( I ) = GRPXL ( J )
            LOCYL ( I ) = GRPYL ( J )
            
        END DO

        MESG = 'Error writing to output file "' // FNAME // '"'

        IF ( .NOT. WRITE3( FNAME, 'ISTACK', SDATE, STIME, LOCGID )) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( FNAME,'LATITUDE',SDATE,STIME,LOCLAT )) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( FNAME,'LONGITUDE',SDATE,STIME,LOCLON )) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( FNAME, 'STKDM', SDATE, STIME, LOCDM ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( FNAME, 'STKHT', SDATE, STIME, LOCHT ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( FNAME, 'STKTK', SDATE, STIME, LOCTK ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( FNAME, 'STKVE', SDATE, STIME, LOCVE ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( FNAME, 'STKFLW', SDATE, STIME, LOCFL ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( FNAME, 'STKCNT', SDATE, STIME, LOCCNT )) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( FNAME, 'ROW', SDATE, STIME, LOCROW ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( FNAME, 'COL', SDATE, STIME, LOCCOL ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( FNAME, 'XLOCA', SDATE, STIME, LOCXL ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( FNAME, 'YLOCA', SDATE, STIME, LOCYL ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )


C...........   Formatted file I/O formats............ 93xxx

93010   FORMAT( 4 I10, F10.2 )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I8, :, 2X  ) )

94015   FORMAT( A, 1X, I6, 1X, A, 1X, I5.5, 1X, A, 1X, I8.8, 1X,
     &          A, 1X, I6, 1X, A, 1X, I6,   1X, A )

94020   FORMAT( 5X, 'H[m]:', 1X, F6.2, 1X, 'D[m]:'  , 1X, F4.2, 1X,
     &              'T[K]:', 1X, F7.1, 1X, 'V[m/s]:', 1X, F10.1 )

        END SUBROUTINE WPINGSTK

