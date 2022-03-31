
        SUBROUTINE RDSRCGRPS( SGDEV, SPCFLAG, TSFLAG )

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
        USE MODMERGE, ONLY: NSRCGRP, IGRPNUM, ISRCGRP,
     &                      EMGGRD, EMGGRDSPC, EMGGRDSPCT, NSGOUTPUT,
     &                      GRPCNT, NGRPS, IUGRPNUM, IUGRPIDX, SUBSECFLAG,
     &                      AFLAG, BFLAG, MFLAG, PFLAG,
     &                      NASRC, NMSRC, NPSRC,
     &                      ANGMAT, AGMATX, MNGMAT, MGMATX, PGMATX,
     &                      NMSPC, NSTEPS,
     &                      AENAME, ASDEV, MENAME, MSDEV, 
     &                      PENAME, PSDEV

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVIFIP, INVCFIP, NINVSCC, INVSCC

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NGRID

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV, ONLY: NGROUP, GROUPID, NELEVGRPS, ELEVGRPID,
     &                     ELEVSTKGRP, ELEVSRCGRP, ELEVSTKCNT, EMELEVGRP

C...........   This module contains the gridding surrogates tables
        USE MODSURG, ONLY: NSRGFIPS, SRGFIPS, NCELLS, FIPCELL

C...........   This module contains the inventory arrays
        USE MODSOURC, ONLY: CIFIP, CSCC, CSOURC

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER         GETFLINE
        INTEGER         STR2INT, FIND1, FINDC, INDEX1
        INTEGER         ENVINT
        LOGICAL         BLKORCMT, CHKINT
        INTEGER         PROMPTFFILE
        LOGICAL         USEEXPGEO
 
        EXTERNAL  GETFLINE, STR2INT, FIND1, ENVINT, BLKORCMT, FINDC,
     &            CHKINT, INDEX1, PROMPTFFILE, USEEXPGEO

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: SGDEV           ! file unit number
        LOGICAL, INTENT (IN) :: SPCFLAG         ! include species in emissions array
        LOGICAL, INTENT (IN) :: TSFLAG          ! include timesteps in emissions array

C...........   Local parameters
        INTEGER, PARAMETER :: MXCOL = 5
        INTEGER, PARAMETER :: SRCLEN = FIPLEN3+SCCLEN3+PLTLEN3+CHRLEN3

C...........   Array for parsing list-formatted inputs
        CHARACTER(50)          SEGMENT( MXCOL )

C.........  Array that contains the names of the inventory variables to read
        CHARACTER(IOVLEN3) IVARNAMS( MXINVARR )

C...........   Local allocatable arrays
        INTEGER,            ALLOCATABLE :: INDEXA  ( : ) ! sorting index
        INTEGER,            ALLOCATABLE :: IGRPNUMA( : ) ! unsorted group number
        CHARACTER(SRCLEN),  ALLOCATABLE :: CGRPSRCA( : ) ! unsorted group characteristics
        CHARACTER(SRCLEN),  ALLOCATABLE :: CGRPSRC ( : ) ! sorted group characteristics
        CHARACTER(16),      ALLOCATABLE :: COMBOGRPS(: ) ! combo groups by source

C...........   Other local variables
        INTEGER         I, J, K, N, INDX, C, F, GIDX, CNT, IGRP, PGRP   !  counters and indices
        INTEGER         MXERR   !  max no. errors of each type
        INTEGER         MXWARN  !  max no. warnings of each type
        INTEGER         NLINES  !  number of lines
        INTEGER         INUM    !  source group number
        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  record counter
        INTEGER         NSRC    !  number of sources or (for biogenics) FIPS
        INTEGER         NINVARR !  number inventory variables to read
        INTEGER         RDEV    !  report output file
        INTEGER         GRPNUM  !  group number for report
        INTEGER         MXGRPNUM ! maximum group number

        LOGICAL      :: EFLAG = .FALSE.   !  true: error found

        CHARACTER(FIPLEN3) CFIP     !  character FIPS code
        CHARACTER(FIPLEN3) CSTA     !  state code
        CHARACTER(SCCLEN3) TSCC     !  SCC code
        CHARACTER(SCCLEN3) PSCC     !  previous SCC code
        CHARACTER(PLTLEN3) CPLTID   !  facility / plant ID code
        CHARACTER(CHRLEN3) CPNTID   !  unit / point ID code
        CHARACTER(SRCLEN)  CSRC     !  source characteristics
        CHARACTER(SRCLEN)  PREVCSRC !  previous source info from sorted list
        CHARACTER(16)      COMBOGRP !  stack group ID + source group index
        CHARACTER(512)     LINE     !  line buffer
        CHARACTER(512)     MESG     !  message buffer

        CHARACTER(16) :: PROGNAME = 'RDSRCGRPS' ! program name

***********************************************************************
C   begin body of subroutine RDSRCGRPS
        
C.........  Read source characteristics needed for matching from inventory

C.........  Check if Movesmrg has already loaded the data it needs
        IF( MFLAG .AND. ASSOCIATED( CSCC ) ) THEN
            NSRC = NMSRC

        ELSE
            IVARNAMS( 1 ) = 'CIFIP'
            IVARNAMS( 2 ) = 'CSCC'
            NINVARR = 2
            IF( AFLAG ) THEN
                NSRC = NASRC
                CALL RDINVCHR( 'AREA', AENAME, ASDEV, NSRC, NINVARR, IVARNAMS )
            
            ELSE IF( BFLAG ) THEN
                NSRC = NSRGFIPS
    
            ELSE IF( MFLAG ) THEN
                NSRC = NMSRC
                CALL RDINVCHR( 'MOBILE', MENAME, MSDEV, NSRC, NINVARR, IVARNAMS )
                
            ELSE IF( PFLAG ) THEN
                NSRC = NPSRC
                IVARNAMS( 3 ) = 'CSOURC'
                NINVARR = 3
                CALL RDINVCHR( 'POINT', PENAME, PSDEV, NSRC, NINVARR, IVARNAMS )
                
            END IF

C.............  Build list of unique SCCs
            IF( .NOT. BFLAG ) THEN

C.................  Create sorting index
                ALLOCATE( INDEXA( NSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INDEXA', PROGNAME )
        
                DO I = 1, NSRC
                    INDEXA( I ) = I
                END DO

C.................  Sort SCCs        
                CALL SORTIC( NSRC, INDEXA, CSCC )

C.................  Count number of unique SCCs
                NINVSCC = 0
                PSCC = ' '
                DO I = 1, NSRC
                
                    TSCC = CSCC( INDEXA( I ) )
                    IF( TSCC .NE. PSCC ) NINVSCC = NINVSCC + 1
                    PSCC = TSCC
        
                END DO

C.................  Build list of unique SCCs        
                ALLOCATE( INVSCC( NINVSCC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INVSCC', PROGNAME )
        
                NINVSCC = 0
                PSCC = ' '
                DO I = 1, NSRC
                
                    TSCC = CSCC( INDEXA( I ) )
                    IF( TSCC .NE. PSCC ) THEN
                        NINVSCC = NINVSCC + 1
                        INVSCC( NINVSCC ) = TSCC
                    END IF
                    PSCC = TSCC
                
                END DO
                
                DEALLOCATE( INDEXA )
            
            END IF  ! build SCC list if not biogenics
        
        END IF  ! need to load inventory data

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
        ALLOCATE( CGRPSRCA( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CGRPSRCA', PROGNAME )

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

            IF( .NOT. CHKINT( SEGMENT( 2 ) ) ) CYCLE  ! skip header lines

            TSCC = SEGMENT( 3 )
            CALL FLTRNEG( TSCC )

            CPLTID = SEGMENT( 4 )
            CALL FLTRNEG( CPLTID )

            CPNTID = SEGMENT( 5 )
            CALL FLTRNEG( CPNTID )

C.............  Skip records that don't apply to current category
            IF( BFLAG .AND. 
     &          ( TSCC /= ' ' .OR. 
     &            CPLTID /= ' ' .OR.
     &            CPNTID /= ' ' ) ) CYCLE
            IF( ( AFLAG .OR. MFLAG ) .AND.
     &          ( CPLTID /= ' ' .OR.
     &            CPNTID /= ' ' ) ) CYCLE

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
            IF( .NOT. USEEXPGEO() .AND.
     &          .NOT. CHKINT( SEGMENT( 2 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Bad FIPS code at line',
     &                 IREC, 'of source grouping file.'
                CALL M3MESG( MESG )
                CYCLE
                        
            ELSE
                CFIP = SEGMENT( 2 )

C.................  Standardize character version of FIPS code                
                CALL FLTRNEG( CFIP )
                CALL PADZERO( CFIP )

C.................  Check if FIPS code matches inventory or 
C                   surrogates (for biogenics)
                IF( USEEXPGEO() .OR. CFIP( FIPEXPLEN3+4:FIPEXPLEN3+6 ) /= '000' ) THEN
                    IF( BFLAG ) THEN
                        J = FINDC( CFIP, NSRGFIPS, SRGFIPS )
                    ELSE
                        J = FINDC( CFIP, NINVIFIP, INVCFIP )
                    END IF

                    IF( J .LE. 0 ) THEN
                        WRITE( MESG,94010 ) 
     &                    'WARNING: State/county FIPS code "' // CFIP //
     &                    '" at line', IREC, 'does not match inventory.'
                        CALL M3MSG2( MESG )
                        CYCLE
                    END IF
                END IF
            END IF

C.............  Standardize SCC code
            CALL PADZERO( TSCC )

C.............  Check SCC code
            IF( .NOT. BFLAG ) THEN

C.................  Check if SCC matches inventory
                IF( TSCC .NE. REPEAT( '0', SCCLEN3 ) ) THEN
                    J = FINDC( TSCC, NINVSCC, INVSCC )
                    
                    IF( J .LE. 0 ) THEN
                        WRITE( MESG,94010 )
     &                    'WARNING: SCC "' // TSCC // '" at line',
     &                    IREC, 'does not match inventory.'
                        CALL M3MSG2( MESG )
                        CYCLE
                    END IF
                END IF

            END IF

C.............  Increment count of valid lines and check it
            N = N + 1
            IF( N .GT. NLINES ) CYCLE  ! Ensure no overflow

C.............  Store fields from source groups file
            INDEXA  ( N ) = N
            IGRPNUMA( N ) = INUM
            CGRPSRCA( N ) = CFIP // TSCC // 
     &                      ADJUSTR( CPLTID ) // ADJUSTR( CPNTID )

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
        ALLOCATE( CGRPSRC( N ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CGRPSRC', PROGNAME )

C.........  Sort source groups by FIPS code
        CALL SORTIC( N, INDEXA, CGRPSRCA )

        PREVCSRC = ' '
        DO I = 1, N

            J = INDEXA( I )

C.............  Check for duplicate group info
            CSRC = CGRPSRCA( J )
            IF( CSRC .EQ. PREVCSRC ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Duplicate source group "' // TRIM( CSRC ) //
     &                 '" in source grouping file.'
                CALL M3MSG2( MESG )
                CYCLE
            END IF

            IGRPNUM( I ) = IGRPNUMA( J )
            CGRPSRC( I ) = CSRC
            
            PREVCSRC = CSRC

        END DO

C.........  Sort source groups by group ID
        CALL SORTI1( N, INDEXA, IGRPNUMA )

C.........  Unique list of groups
        NGRPS = 0
        PGRP = -1
        DO I = 1, N
            J = INDEXA( I )
            IGRP = IGRPNUMA( J )
            IF( PGRP .NE. IGRP ) THEN
                NGRPS = NGRPS + 1
            END IF
            PGRP = IGRP
        END DO
        MXGRPNUM = PGRP

        ALLOCATE( IUGRPNUM( NGRPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IUGRPNUM', PROGNAME )
        ALLOCATE( IUGRPIDX( MXGRPNUM ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IUGRPIDX', PROGNAME )
        IUGRPNUM = 0
        IUGRPIDX = 0

        NGRPS = 0
        PGRP = -1
        DO I = 1, N
            J = INDEXA( I )
            IGRP = IGRPNUMA( J )
            IF( PGRP .NE. IGRP ) THEN
                NGRPS = NGRPS + 1
                IUGRPNUM( NGRPS ) = IGRP
                IUGRPIDX( IGRP ) = NGRPS
            END IF
            PGRP = IGRP
        END DO

C.........  Check for sorting errors
        IF( EFLAG ) THEN
            MESG = 'Problem with source groups file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Set default source group        
        IGRPNUM( NSRCGRP ) = 0

C.........  Output source groups report
        RDEV = PROMPTFFILE(
     &           'Enter name for SOURCE GROUPS REPORT file', 
     &           .FALSE., .TRUE., 'SRCGRP_REPORT', PROGNAME )

C.........  Assign sources to source groups
        ALLOCATE( ISRCGRP( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISRCGRP', PROGNAME )
        ISRCGRP = 0
        
        DO I = 1, NSRC
            
            IF( BFLAG ) THEN
                CFIP = SRGFIPS( I )
            ELSE
                CFIP = CIFIP( I )
            END IF

            CSTA = CFIP( 1:STALEN3 ) // '000'
            
            INDX = 0
            CSRC = ' '
            
            IF( PFLAG ) THEN

C.................  Match plant characteristics
                TSCC = CSCC( I )
                CALL PADZERO( TSCC )
                
                CPLTID = CSOURC( I )( PLTPOS3:PLTPOS3+PLTLEN3 )
                CPNTID = CSOURC( I )( CH1POS3:CH1POS3+CHRLEN3 )

C.................  full FIPS, SCC, plant, point
                CSRC = CFIP // TSCC // CPLTID // CPNTID
                INDX = FINDC( CSRC, N, CGRPSRC )
                
C.................  full FIPS, SCC, plant
                IF( INDX .LT. 0 ) THEN
                    CSRC = CFIP // TSCC // CPLTID
                    INDX = FINDC( CSRC, N, CGRPSRC )
                END IF
                
C.................  full FIPS, blank SCC, plant, point
                TSCC = REPEAT( '0', SCCLEN3 )
                IF( INDX .LT. 0 ) THEN
                    CSRC = CFIP // TSCC // CPLTID // CPNTID
                    INDX = FINDC( CSRC, N, CGRPSRC )
                END IF
                
C.................  full FIPS, blank SCC, plant
                IF( INDX .LT. 0 ) THEN
                    CSRC = CFIP // TSCC // CPLTID
                    INDX = FINDC( CSRC, N, CGRPSRC )
                END IF

            END IF
            
            IF( INDX .LE. 0 .AND. .NOT. BFLAG ) THEN

C.................  Match SCC
                TSCC = CSCC( I )
                CALL PADZERO( TSCC )
                
C.................  full FIPS, SCC
                CSRC = CFIP // TSCC
                INDX = FINDC( CSRC, N, CGRPSRC )
                
C.................  state, SCC
                IF( INDX .LT. 0 .AND. .NOT. USEEXPGEO() ) THEN
                    CSRC = CSTA // TSCC
                    INDX = FINDC( CSRC, N, CGRPSRC )
                END IF
                
C.................  blank FIPS, SCC
                IF( INDX .LT. 0 ) THEN
                    CSRC = REPEAT( '0', FIPLEN3 ) // TSCC
                    INDX = FINDC( CSRC, N, CGRPSRC )
                END IF
                
            END IF
            
            IF( INDX .LE. 0 ) THEN

C.................  full FIPS
                CSRC = CFIP // REPEAT( '0', SCCLEN3 )
                INDX = FINDC( CSRC, N, CGRPSRC )
                
C.................  state
                IF( INDX .LT. 0 .AND. .NOT. USEEXPGEO() ) THEN
                    CSRC = CSTA // REPEAT( '0', SCCLEN3 )
                    INDX = FINDC( CSRC, N, CGRPSRC )
                END IF

C.................  user-assigned default group
                IF( INDX .LT. 0 ) THEN
                    CSRC = REPEAT( '0', FIPLEN3) // REPEAT( '0', SCCLEN3 )
                    INDX = FINDC( CSRC, N, CGRPSRC )
                END IF
                
            END IF
            
            IF( INDX .GT. 0 ) THEN
                ISRCGRP( I ) = INDX
            ELSE

C.................  If no match, assign to default group
                ISRCGRP( I ) = NSRCGRP
            END IF

C.............  Add source to report file
            GRPNUM = IGRPNUM( ISRCGRP( I ) )

            IF( AFLAG .OR. MFLAG ) THEN

                WRITE( RDEV,'(I8,1X,A,1X,A,1X,I8)' ) 
     &            I, CIFIP( I ), CSCC( I ), GRPNUM
                
            ELSE IF( BFLAG ) THEN
            
                WRITE( RDEV,'(A,1X,I8)' ) SRGFIPS( I ), GRPNUM
            
            ELSE IF( PFLAG ) THEN
            
                CPLTID = CSOURC( I )( PLTPOS3:PLTPOS3+PLTLEN3 )
                CPNTID = CSOURC( I )( CH1POS3:CH1POS3+CHRLEN3 )

                WRITE( RDEV,'(I8,1X,A,1X,A,1X,A,1X,A,1X,I8)' ) 
     &            I, CIFIP( I ), CSCC( I ), CPLTID, CPNTID, GRPNUM

            END IF
        
        END DO

C.........  Reassign elevated stack groups
        IF( PFLAG .AND. NGROUP .NE. 0 ) THEN

            ALLOCATE( ELEVGRPID( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ELEVGRPID', PROGNAME )
            ALLOCATE( COMBOGRPS( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'COMBOGRPS', PROGNAME )

            ELEVGRPID = 0

            K = 0
            DO I = 1, NSRC

C.................  Skip sources that aren't elevated
                IF( GROUPID( I ) == 0 ) CYCLE
            
C.................  Build combo group ID using stack group number and 
C                   source group index for current source
                WRITE( COMBOGRP, '(I8.8, I8.8)' ) GROUPID( I ), ISRCGRP( I )

C.................  Check if combo group has already been assigned
                INDX = INDEX1( COMBOGRP, K, COMBOGRPS )
                IF( INDX > 0 ) THEN

C.....................  Assign source to existing group
                    ELEVGRPID( I ) = INDX
                
                ELSE

C.....................  Create new combo group
                    K = K + 1
                    COMBOGRPS( K ) = COMBOGRP
                    ELEVGRPID( I ) = K

                END IF
            
            END DO

            NELEVGRPS = K
            
            DEALLOCATE( COMBOGRPS )

C.............  Build lists mapping new elevated groups to original stack groups
C               and source groups
            ALLOCATE( ELEVSTKGRP( NELEVGRPS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ELEVSTKID', PROGNAME )
            ALLOCATE( ELEVSRCGRP( NELEVGRPS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ELEVSRCGRP', PROGNAME )
            ALLOCATE( ELEVSTKCNT( NELEVGRPS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ELEVSTKCNT', PROGNAME )
            ELEVSTKCNT = 0
            
            DO I = 1, NSRC

C.................  Skip sources that aren't elevated
                IF( GROUPID( I ) == 0 ) CYCLE
            
                INDX = ELEVGRPID( I )
                ELEVSTKGRP( INDX ) = GROUPID( I )
                ELEVSRCGRP( INDX ) = ISRCGRP( I )
                ELEVSTKCNT( INDX ) = ELEVSTKCNT( INDX ) + 1
            
            END DO

C.............  Allocate storage for summed emissions per group            
            ALLOCATE( EMELEVGRP( NELEVGRPS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EMELEVGRP', PROGNAME )
            EMELEVGRP = 0.

        END IF

C.........  Determine number of source group / grid cell interactions
        NSGOUTPUT = 0

        ALLOCATE( GRPCNT( NGRID, NSRCGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPCNT', PROGNAME )
        GRPCNT = 0  ! array

        IF( AFLAG ) THEN
            CALL SRCGRPCNT( NSRC, ANGMAT, AGMATX( 1 ), 
     &                      AGMATX( NGRID+1 ) )
        END IF
        
        IF( BFLAG ) THEN

C.............  Loop through FIPS codes from surrogates
            DO F = 1, NSRGFIPS

                GIDX = ISRCGRP( F )

C.................  Loop through grid cells for FIPS code
                DO N = 1, NCELLS( F )

                    C = FIPCELL( N,F )

                    CNT = GRPCNT( C, GIDX )
                    IF( CNT == 0 ) THEN
                        NSGOUTPUT = NSGOUTPUT + 1
                    END IF
                    GRPCNT( C, GIDX ) = CNT + 1

                END DO  ! end loop over grid cells for FIPS

            END DO  ! end loop over FIPS codes

        END IF
        
        IF( MFLAG ) THEN
            CALL SRCGRPCNT( NSRC, MNGMAT, MGMATX( 1 ), 
     &                      MGMATX( NGRID+1 ) )
        END IF
        
        IF( PFLAG ) THEN
            CALL SRCGRPCNT( NSRC, NSRC, PGMATX( 1 ), 
     &                      PGMATX( NGRID+1 ) )

C.............  Increment number of output records to account
C               for elevated sources
            NSGOUTPUT = NSGOUTPUT + NELEVGRPS
        END IF

C.........  Upate no of srcgrp when SUBSECFLAG output is set to Y
        IF( SUBSECFLAG ) NSRCGRP = NGRPS

C.........  Allocate memory for emissions data
        ALLOCATE( EMGGRD( NGRID, NSRCGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMGGRD', PROGNAME )
        EMGGRD = 0.  ! array

C.........  Allocate additional arrays for Movesmrg processing
        IF( SPCFLAG ) THEN
            IF( TSFLAG ) THEN
                ALLOCATE( EMGGRDSPCT( NGRID, NSRCGRP, NMSPC, NSTEPS ), STAT=IOS )
                CALL CHECKMEM( IOS, 'EMGGRDSPCT', PROGNAME )
                EMGGRDSPCT = 0.  ! array
            ELSE
                ALLOCATE( EMGGRDSPC( NGRID, NSRCGRP, NMSPC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'EMGGRDSPC', PROGNAME )
                EMGGRDSPC = 0.  ! array
            END IF
        END IF

C.........  Deallocate local memory
        DEALLOCATE( INDEXA, IGRPNUMA, CGRPSRCA, CGRPSRC )

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

