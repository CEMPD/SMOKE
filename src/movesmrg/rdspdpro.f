
        SUBROUTINE RDSPDPRO( SPDEV )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C       Reads hourly speed data from SPDPRO file
C
C  PRECONDITIONS REQUIRED:
C       SPDEV must be opened
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     09/10: Created by C. Seppanen
C     08/15: Modified by B.H. Baek 
C
C***********************************************************************
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***********************************************************************

C.........  MODULES for public variables
C.........  This module contains data structures and flags specific to Movesmrg
        USE M3UTILIO

        USE MODMVSMRG, ONLY: SPDPRO

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVIFIP, INVCFIP, NINVSCC, INVSCC

C.........  This module is for mobile-specific data
        USE MODMOBIL, ONLY: NSCCMAP, SCCMAPFLAG, SCCMAPLIST, EXCLSCCFLAG

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL       BLKORCMT
        LOGICAL       CHKINT
        LOGICAL       CHKREAL
C       LOGICAL       ENVYN
C       INTEGER       INDEX1 
        INTEGER       GETFLINE
C       INTEGER       FINDC
C       INTEGER       STR2INT
C       INTEGER       PROMPTFFILE
C       INTEGER       ENVINT
C       REAL          STR2REAL
C       CHARACTER(2)  CRLF
        
C        EXTERNAL BLKORCMT, CHKINT, CHKREAL, FINDC, GETFLINE, ENVINT, 
C     &           STR2INT, STR2REAL, CRLF, ENVYN, INDEX1, PROMPTFFILE
        EXTERNAL     BLKORCMT, CHKINT, CHKREAL, GETFLINE

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: SPDEV             ! SPDPRO file unit no.

C...........   Local allocatable arrays

C...........   Local arrays
        CHARACTER(20)  SEGMENT( 53 )          ! parsed input line

C...........   Other local variables
        INTEGER         I, J, JJ, KK     ! counters and indexes
        INTEGER         IOS         ! error status
        INTEGER      :: IREC = 0    ! record counter
        INTEGER         FIPIDX      ! current FIPS index
        INTEGER         NSCC        ! no of reference SCCs
        INTEGER         SCCIDX      ! current SCC index
        INTEGER         NLINES      ! number of lines
        INTEGER         MDEV        ! open SCCXREF input file
        INTEGER     NWARN1, NWARN2, MXWARN   ! number of warning messages
 
        REAL            SPDVAL      ! hourly speed value

        LOGICAL      :: EFLAG = .FALSE.   ! true: error found

        CHARACTER(1060)    LINE     ! line buffer
        CHARACTER(300)     MESG     ! message buffer
        CHARACTER(FIPLEN3) CFIP     ! current FIPS
        CHARACTER(SCCLEN3) SCC      ! current SCC
        CHARACTER(10)      KEYWORD  ! temperature keyword

        CHARACTER(16) :: PROGNAME = 'RDSPDPRO'   ! program name

C***********************************************************************
C   begin body of subroutine RDSPDPRO

C.........  Get maximum number of warnings
        MXWARN = ENVINT( WARNSET, ' ', 100, IOS )

C.........  Read cross-refrenced SCC input file
        MESG = 'Use referenced SCC activity inventory file'
        SCCMAPFLAG = ENVYN ( 'USE_REF_SCC_YN', MESG, .FALSE., IOS )

        IF( SCCMAPFLAG ) THEN
            MESG = 'Enter logical name for reference SCC input file'
            MDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'SCCXREF',
     &                      PROGNAME )
            CALL RDSCCMAP( MDEV )

            MESG = 'Exclude SCCs not found in SCCXREF input file'
            EXCLSCCFLAG = ENVYN ( 'EXCLUDE_REF_SCC_YN', MESG, .FALSE., IOS )

        END IF

C.........  Allocate storage based on number of FIPs and SCCs in inventory
        ALLOCATE( SPDPRO( NINVIFIP, NINVSCC, 2, 24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPDPRO', PROGNAME )
        SPDPRO = BADVAL3   ! array

C.........  Get the number of lines in the file
        NLINES = GETFLINE( SPDEV, 'SPDPRO file' )

C.........  Read through file and store hourly data
        NWARN1 = 0
        NWARN2 = 0
        DO I = 1, NLINES
        
          READ( SPDEV, 93000, END=999, IOSTAT=IOS ) LINE
            
          IREC = IREC + 1
            
          IF( IOS .NE. 0 ) THEN
              EFLAG = .TRUE.
              WRITE( MESG, 94010 ) 'I/O error', IOS,
     &          'reading hourly speed file at line', IREC
              CALL M3MESG( MESG )
              CYCLE
          END IF

C...........  Skip blank or comment lines
          IF( BLKORCMT( LINE ) ) CYCLE

C...........  Parse line into fields
          CALL PARSLINE( LINE, 53, SEGMENT )

C...........  SCC mapping loop based on SCCXREF reference input file
          SCC = ADJUSTL( SEGMENT( 2 ) )
          CALL PADZERO( SCC )

          KK = 0
          NSCC = 0
          IF( SCCMAPFLAG ) THEN
              KK   = INDEX1( SCC, NSCCMAP, SCCMAPLIST( :,1 ) )
              IF( KK > 0 ) THEN
                  NSCC = STR2INT( SCCMAPLIST( KK,3 ) )
              ELSE
                  IF( EXCLSCCFLAG ) CYCLE    ! drop SCCs not listed in SCCXREF file
              END IF
          END IF

          DO JJ = 0, NSCC

            IF( JJ > 0 .AND. KK > 0 ) IREC = IREC + 1
            IF( SCCMAPFLAG .AND. KK > 0 ) SCC = SCCMAPLIST( KK+JJ,2 )

C.............  Find SCC in inventory list
            SCCIDX = FINDC( SCC, NINVSCC, INVSCC )
            IF( SCCIDX .LE. 0 ) THEN
                WRITE( MESG, 94010 ) 'NOTE: Skipping line',
     &            IREC, 'of hourly speed file because SCC ' //
     &            SCC // ' is not in the inventory.'
                NWARN1 = NWARN1 + 1
                IF( NWARN1 < MXWARN ) CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Find county in inventory list
            CFIP = ADJUSTL( SEGMENT( 1 ) )
            CALL PADZERO( CFIP )
            FIPIDX = FINDC( CFIP, NINVIFIP, INVCFIP )

            IF( FIPIDX .LE. 0 ) THEN
                WRITE( MESG, 94010 ) 'NOTE: Skipping line',
     &            IREC, 'of hourly speed file because FIPS code '//
     &            CFIP //' is not in the inventory.'
                NWARN2 = NWARN2 + 1
                IF( NWARN2 < MXWARN ) CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Check weekday keyword
            KEYWORD = ADJUSTL( SEGMENT( 3 ) )
            IF( KEYWORD .NE. 'WEEKDAY' ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Missing expected ' //
     &            'keyword WEEKDAY at line', IREC,
     &            'of hourly speed file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Check weekday speed values
            DO J = 1, 24
                IF( .NOT. CHKREAL( SEGMENT( 3 + J ) ) ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 ) 'ERROR: Bad weekday ' //
     &                'hourly speed value at line', IREC,
     &                'of hourly speed file.'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF
                
                SPDVAL = STR2REAL( SEGMENT( 3 + J ) )                
                SPDPRO( FIPIDX, SCCIDX, 2, J ) = SPDVAL
            END DO

C.............  Check weekend keyword
            KEYWORD = ADJUSTL( SEGMENT( 28 ) )
            IF( KEYWORD .NE. 'WEEKEND' ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Missing expected ' //
     &            'keyword WEEKEND at line', IREC,
     &            'of hourly speed file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Check weekend speed values
            DO J = 1, 24
                IF( .NOT. CHKREAL( SEGMENT( 28 + J ) ) ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 ) 'ERROR: Bad weekend ' //
     &                'hourly speed value at line', IREC,
     &                'of hourly speed file.'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF
                
                SPDVAL = STR2REAL( SEGMENT( 28 + J ) )
                SPDPRO( FIPIDX, SCCIDX, 1, J ) = SPDVAL
            END DO
            
          END DO    ! SCCXREF loop

        END DO

        CLOSE( SPDEV )
        
        IF( EFLAG ) THEN
            MESG = 'Problem found in hourly speed file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
        RETURN

999     MESG = 'End of file'
        MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of SPDPRO' // CRLF() // BLANK5 //
     &         'input file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )  
      
C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
        
        END SUBROUTINE RDSPDPRO

