
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
C***********************************************************************

C.........  MODULES for public variables
C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: SPDPRO

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVIFIP, INVIFIP, NINVSCC, INVSCC

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL       BLKORCMT
        LOGICAL       CHKINT
        LOGICAL       CHKREAL
        INTEGER       GETFLINE
        INTEGER       FIND1
        INTEGER       FINDC
        INTEGER       STR2INT
        REAL          STR2REAL
        CHARACTER(2)  CRLF
        
        EXTERNAL BLKORCMT, CHKINT, CHKREAL, FIND1, GETFLINE, 
     &           STR2INT, STR2REAL, CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: SPDEV             ! SPDPRO file unit no.

C...........   Local allocatable arrays

C...........   Local arrays
        CHARACTER(20)  SEGMENT( 53 )          ! parsed input line

C...........   Other local variables
        INTEGER         I, J        ! counters and indexes
        INTEGER         IOS         ! error status
        INTEGER      :: IREC = 0    ! record counter
        INTEGER         FIPIDX      ! current FIPS index
        INTEGER         SCCIDX      ! current SCC index
        INTEGER         NLINES      ! number of lines
        INTEGER         CNTY        ! current FIPS code
        
        REAL            SPDVAL      ! hourly speed value

        LOGICAL      :: EFLAG = .FALSE.   ! true: error found

        CHARACTER(1060)     LINE     ! line buffer
        CHARACTER(300)     MESG     ! message buffer
        CHARACTER(SCCLEN3) SCC      ! current SCC
        CHARACTER(10)      KEYWORD  ! temperature keyword

        CHARACTER(16) :: PROGNAME = 'RDSPDPRO'   ! program name

C***********************************************************************
C   begin body of subroutine RDSPDPRO

C.........  Allocate storage based on number of FIPs and SCCs in inventory
        ALLOCATE( SPDPRO( NINVIFIP, NINVSCC, 2, 24 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPDPRO', PROGNAME )
        SPDPRO = BADVAL3   ! array

C.........  Get the number of lines in the file
        NLINES = GETFLINE( SPDEV, 'SPDPRO file' )

C.........  Read through file and store hourly data
        DO I = 1, NLINES
        
            READ( SPDEV, 93000, END=999, IOSTAT=IOS ) LINE
            
            IREC = IREC + 1
            
            IF( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'I/O error', IOS,
     &            'reading hourly speed file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Skip blank or comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Parse line into fields
            CALL PARSLINE( LINE, 53, SEGMENT )

C.............  Convert FIP to integer
            IF( .NOT. CHKINT( SEGMENT( 1 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Bad FIPS code ' //
     &            'at line', IREC, 'of hourly speed file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Find county in inventory list
            CNTY = STR2INT( ADJUSTR( SEGMENT( 1 ) ) )
            FIPIDX = FIND1( CNTY, NINVIFIP, INVIFIP )
            
            IF( FIPIDX .LE. 0 ) THEN
                WRITE( MESG, 94010 ) 'NOTE: Skipping line',
     &            IREC, 'of hourly speed file because FIPS code',
     &            CNTY, 'is not in the inventory.'
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Find SCC in inventory list
            SCC = ADJUSTL( SEGMENT( 2 ) )
            SCCIDX = FINDC( SCC, NINVSCC, INVSCC )
            
            IF( SCCIDX .LE. 0 ) THEN
                WRITE( MESG, 94010 ) 'NOTE: Skipping line',
     &            IREC, 'of hourly speed file because SCC ' //
     &            SCC // ' is not in the inventory.'
                CALL M3MESG( MESG )
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
                
                SPDVAL = STR2REAL( ADJUSTR( SEGMENT( 3 + J ) ) )                
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
                
                SPDVAL = STR2REAL( ADJUSTR( SEGMENT( 28 + J ) ) )
                SPDPRO( FIPIDX, SCCIDX, 1, J ) = SPDVAL
            END DO
            
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

