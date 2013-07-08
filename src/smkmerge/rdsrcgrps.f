
        SUBROUTINE RDSRCGRPS( SGDEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C     This subroutine reads the source grouping file and assigns sources
C     to groups.
C
C  PRECONDITIONS REQUIRED:
C     File unit SGDEV is open
C     State and country codes have been loaded by INITSTCY
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 5/2013 by Dongmei Yang
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2013, Environmental Modeling for Policy Development
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

C.........  MODULES for public variables
        USE MODMERGE, ONLY: NSRCGRP, IGRPNUM, IFIPGRP, EMGGRD, NSGOUTPUT,
     &                      GRPCNT,
     &                      AFLAG, BFLAG, MFLAG, PFLAG,
     &                      NASRC, NMSRC, NPSRC,
     &                      ANGMAT, AGMATX, MNGMAT, MGMATX, PGMATX

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVIFIP, INVIFIP

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: AICNY, MICNY, PICNY

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NGRID

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV, ONLY: NGROUP
        
        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER         GETFLINE
        INTEGER         STR2INT, FIND1, FINDC
        INTEGER         ENVINT
        LOGICAL         BLKORCMT, CHKINT
 
        EXTERNAL  GETFLINE, STR2INT, FIND1, ENVINT, BLKORCMT, FINDC,
     &            CHKINT

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: SGDEV           ! file unit number

C...........   Local parameters
        INTEGER, PARAMETER :: MXCOL = 4

C...........   Array for parsing list-formatted inputs
        CHARACTER(50)          SEGMENT( MXCOL )

C...........   Local allocatable arrays
        INTEGER,            ALLOCATABLE :: INDEXA  ( : ) ! sorting index
        INTEGER,            ALLOCATABLE :: IGRPNUMA( : ) ! unsorted group number
        CHARACTER(FIPLEN3), ALLOCATABLE :: CGRPFIPA( : ) ! unsorted FIPS
        CHARACTER(FIPLEN3), ALLOCATABLE :: CGRPFIP ( : ) ! sorted FIPS

C...........   Other local variables
        INTEGER         I, J, N, INDX   !  counters and indices
        INTEGER         MXERR   !  max no. errors of each type
        INTEGER         MXWARN  !  max no. warnings of each type
        INTEGER         NLINES  !  number of lines
        INTEGER         INUM    !  source group number
        INTEGER         IFIP    !  integer FIPS code
        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  record counter

        LOGICAL      :: EFLAG = .FALSE.   !  true: error found

        CHARACTER(FIPLEN3) CFIP     !  character FIPS code
        CHARACTER(FIPLEN3) CSTA     !  state code
        CHARACTER(FIPLEN3) PREVFIP  !  previous FIPS code from sorted list
        CHARACTER(512)     LINE     !  line buffer
        CHARACTER(512)     MESG     !  message buffer

        CHARACTER(16) :: PROGNAME = 'RDSRCGRPS' ! program name

***********************************************************************
C   begin body of subroutine RDSRCGRPS

C.........  Write status message
        MESG = 'Reading source groups file...'
        CALL M3MSG2( MESG )

C.........  Get the number of lines in the file
        NLINES = GETFLINE( SGDEV, 'Source groups file' )

C.........  Allocate memory for unsorted data
        ALLOCATE( INDEXA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDEXA', PROGNAME )
        ALLOCATE( IGRPNUMA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IGRPNUMA', PROGNAME )
        ALLOCATE( CGRPFIPA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CGRPFIPA', PROGNAME )

C.........  Read lines and store unsorted data
        IREC   = 0
        N      = 0

        DO I = 1, NLINES

            READ( SGDEV, 93000, END=999, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'I/O error', IOS,
     &              'reading source groups file at line', IREC
                CALL M3MSG2( MESG )
                CYCLE
            END IF

C.............  Skip blank lines or comments
            IF( BLKORCMT( LINE ) ) CYCLE

            CALL PARSLINE( LINE, MXCOL, SEGMENT )

C.............  Check group number
            IF( CHKINT( SEGMENT( 1 ) ) ) THEN
                INUM = STR2INT( SEGMENT( 1 ) )

C.................  Reserve group number zero
                IF( INUM == 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Group number 0 ' //
     &                     'at line', IREC, 'of source grouping file.'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

            ELSE
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Bad group number at line',
     &                 IREC, 'of source grouping file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Check FIPS code
            IF( SEGMENT( 2 ) == ' ' ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Missing FIPS code at line',
     &                 IREC, 'of source grouping file.'
                CALL M3MESG( MESG )
                CYCLE

            ELSE IF( CHKINT( SEGMENT( 2 ) ) ) THEN
                CFIP = SEGMENT( 2 )
                IFIP = STR2INT( CFIP )

C.................  Standardize character version of FIPS code                
                CALL FLTRNEG( CFIP )
                CALL PADZERO( CFIP )

C.................  Check if FIPS code matches inventory
                IF( CFIP( 4:6 ) /= '000' ) THEN
                    J = FIND1( IFIP, NINVIFIP, INVIFIP )
                    IF( J .LE. 0 ) THEN
                        WRITE( MESG,94010 ) 
     &                    'WARNING: State/county FIPS code "' // CFIP //
     &                    '" at line', IREC, 'does not match inventory.'
                        CALL M3MSG2( MESG )
                        CYCLE
                    END IF
                END IF
                        
            ELSE
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Bad FIPS code at line',
     &                 IREC, 'of source grouping file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Increment count of valid lines and check it
            N = N + 1
            IF( N .GT. NLINES ) CYCLE  ! Ensure no overflow

C.............  Store fields from source groups file
            INDEXA  ( N ) = N
            IGRPNUMA( N ) = INUM
            CGRPFIPA( N ) = CFIP

        END DO      ! End of loop on I for reading in source groups file

C.........  Write warning if no groups match inventory
        IF( .NOT. EFLAG .AND. N .EQ. 0 ) THEN
            MESG = 'WARNING: All sources assigned to default group.'
            CALL M3MSG2( MESG )
        END IF

C.........  Check for errors reading file, and abort
        IF( EFLAG ) THEN
            MESG = 'Problem reading source groups file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        MESG = 'Processing source groups file...'
        CALL M3MSG2( MESG )

C.........  Save total number of source groups; add additional spot
C           for default group
        NSRCGRP = N + 1

C.........  Allocate memory for sorted groups
        ALLOCATE( IGRPNUM( NSRCGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IGRPNUM', PROGNAME )
        ALLOCATE( CGRPFIP( N ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CGRPFIP', PROGNAME )

C.........  Sort source groups by FIPS code
        CALL SORTIC( N, INDEXA, CGRPFIPA )

        PREVFIP = ' '
        DO I = 1, N

            J = INDEXA( I )
            
C.............  Check for duplicate FIPS codes
            CFIP = CGRPFIPA( J )
            IF( CFIP .EQ. PREVFIP ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Duplicate FIPS code "' // CFIP // 
     &                 '" in source grouping file.'
                CALL M3MSG2( MESG )
                CYCLE
            END IF

            IGRPNUM( I ) = IGRPNUMA( J )
            CGRPFIP( I ) = CFIP
            
            PREVFIP = CFIP

        END DO

C.........  Check for sorting errors
        IF( EFLAG ) THEN
            MESG = 'Problem with source groups file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
        IGRPNUM( NSRCGRP ) = 0

C.........  Assign FIPS from inventory to source groups
        ALLOCATE( IFIPGRP( NINVIFIP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IFIPGRP', PROGNAME )

        DO I = 1, NINVIFIP
        
            IFIP = INVIFIP( I )
            WRITE( CFIP, '(I5.5)' ) IFIP
            CALL PADZERO( CFIP )  ! pad with zeros
            CSTA = CFIP( 1:STALEN3 ) // '000'
            
            INDX = FINDC( CFIP, N, CGRPFIP )

C.............  If county-level FIPS isn't found, check state level
            IF( INDX .LT. 0 ) THEN
                INDX = FINDC( CSTA, N, CGRPFIP )
            END IF
            
            IF( INDX .GT. 0 ) THEN
                IFIPGRP( I ) = INDX
            ELSE

C.................  If no state or county match, assign to default group
                IFIPGRP( I ) = NSRCGRP
            END IF
        
        END DO

C.........  Determine number of source group / grid cell interactions
        NSGOUTPUT = 0

        ALLOCATE( GRPCNT( NGRID, NSRCGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPCNT', PROGNAME )
        GRPCNT = 0  ! array

        IF( AFLAG ) THEN
            CALL SRCGRPCNT( NASRC, ANGMAT, AGMATX( 1 ), 
     &                      AGMATX( NGRID+1 ), AICNY )
        END IF
        
        IF( BFLAG ) THEN
C.........  TODO: add biogenic handling
        END IF
        
        IF( MFLAG ) THEN
            CALL SRCGRPCNT( NMSRC, MNGMAT, MGMATX( 1 ), 
     &                      MGMATX( NGRID+1 ), MICNY )
        END IF
        
        IF( PFLAG ) THEN
            CALL SRCGRPCNT( NPSRC, NPSRC, PGMATX( 1 ), 
     &                      PGMATX( NGRID+1 ), PICNY )

C.............  Increment number of output records to account
C               for elevated sources
            NSGOUTPUT = NSGOUTPUT + NGROUP
        END IF

C.........  Allocate memory for emissions data
        ALLOCATE( EMGGRD( NGRID, NSRCGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMGGRD', PROGNAME )
        EMGGRD = 0.  ! array

C.........  Deallocate local memory
        DEALLOCATE( INDEXA, IGRPNUMA, CGRPFIPA, CGRPFIP )

        RETURN

C.........  Unexpected end of file
999     MESG = 'INTERNAL ERROR: Unexpected end of source grouping file'
        CALL M3MSG2( MESG )

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDSRCGRPS

