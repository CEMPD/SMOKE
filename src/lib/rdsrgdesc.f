
        SUBROUTINE RDSRGDESC( FDEV )

C**************************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine allocates memory for surrogate code, region, 
C      description and file list from SRGDESC file.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 8/2005 by B. Baek 
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
C**************************************************************************

C...........   Modules for public variables

C...........   This module contains the gridding surrogates desciption files
        USE MODSURG, ONLY: NSRGS, IDXSRGA, SRGLIST, SRGFCOD, SRGFREG,
     &                     SRGFDES, SRGFNAM, NTSRGDSC
     
        IMPLICIT NONE

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER        GETFLINE
        LOGICAL        BLKORCMT
        INTEGER        PROMPTFFILE
        INTEGER        STR2INT
        
        EXTERNAL       BLKORCMT, GETFLINE, PROMPTFFILE, STR2INT

C...........   Subroutine arguments
        INTEGER, INTENT ( IN )  :: FDEV       ! file unit number

C...........   Local variables
        INTEGER         I, J, K, N, NT           ! indices and counters

        INTEGER         MDEV                  ! for surrogate files
        INTEGER         ENDLEN                ! end length for reading descriptn
        INTEGER         IOS                   ! i/o status
        INTEGER         NSRGALL               ! No. entries in surrgoates file
        INTEGER      :: IREC = 0              ! record number
        INTEGER      :: NLINES = 0            ! number of lines in input SRGDESC file
        INTEGER         SSC                   ! temporary spatial surrogate code
        INTEGER         LSSC                  ! previous temporary SSC

        CHARACTER(256)  LINE                  ! Read buffer for a line
        CHARACTER(300)  MESG                  ! Message buffer
        CHARACTER(60)   SEGMENT( 4 )           ! line parsing array

        CHARACTER(16) :: PROGNAME = 'RDSRGDESC'    !  program name

C***********************************************************************
C   Begin body of subroutine RDSRGDESC

        REWIND( FDEV )  ! In case of multiple calls
        
C.........  Get the number of lines in the surrogate description file
        NLINES = GETFLINE( FDEV, 'SRGDESC Descriptions' )

C.........  Allocate memory for the SRGDESC descriptions and initialize
        ALLOCATE( IDXSRGA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDXSRGA', PROGNAME )
        ALLOCATE( SRGFREG( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRGFREG', PROGNAME )
        ALLOCATE( SRGFCOD( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRGFCOD', PROGNAME )
        ALLOCATE( SRGFDES( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRGFDES', PROGNAME )
        ALLOCATE( SRGFNAM( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRGFNAM', PROGNAME )
        
C.........  Read surrogate files in SRGDESC file and store
        N = 0
        IREC = 0       
        DO I = 1, NLINES

            READ ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010)
     &                'I/O error', IOS, 'reading SRGDESC '//
     &                'surrogate description file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Left adjust line
            LINE = ADJUSTL( LINE )

C.............  Skip blank and comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Get line
            CALL PARSLINE( LINE, 4, SEGMENT )

            N = N + 1      ! actual line #s after skipping blank and comments
            
            IDXSRGA( N ) = N
            SRGFREG( N ) = ADJUSTL( SEGMENT( 1 ) )
            SRGFCOD( N ) = STR2INT( SEGMENT( 2 ) )
            SRGFDES( N ) = ADJUSTL( SEGMENT( 3 ) )
            SRGFNAM( N ) = ADJUSTL( SEGMENT( 4 ) )

999     END DO  

        NT = N
        NTSRGDSC = NT

C.........  Sort surrogates by county code & cell & surrogate code
        CALL SORTI2( NTSRGDSC, IDXSRGA, SRGFCOD, SRGFNAM )
        
C.........  Count the surrogate codes and no of surrogate files
        LSSC = -1
        NSRGS = 0
        
        DO I = 1, NT

            J   = IDXSRGA( I )
            SSC = SRGFCOD( J )

            IF( SSC .NE. LSSC ) THEN
                NSRGS = NSRGS + 1
                LSSC = SSC
            END IF
               
        END DO

C.........  Allocate memory for derived surrogates tables
        ALLOCATE( SRGLIST( NSRGS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRGLIST', PROGNAME )

C.........  Store the surrgage codes list and names of surrogate file
        LSSC = -1
        NSRGS = 0
        DO I = 1, NT
       
            J   = IDXSRGA( I )
            SSC = SRGFCOD( J )
            
            IF( SSC .NE. LSSC ) THEN
                NSRGS = NSRGS + 1
                SRGLIST( NSRGS ) = SSC
                LSSC = SSC
            END IF

        END DO
        
C.........  Successful completion
        RETURN

C.........  Unexpected end of file
998     MESG = 'INTERNAL ERROR : Unexpected end of surrogate file '
        CALL M3MSG2( MESG )

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

C******************  FORMAT  STATEMENTS   ******************************

93000   FORMAT( A )

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDSRGDESC
                             