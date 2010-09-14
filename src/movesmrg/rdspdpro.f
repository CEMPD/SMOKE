
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
        CHARACTER(20)  SEGMENT( 50 )          ! parsed input line

C...........   Other local variables
        INTEGER         I, J        ! counters and indexes
        INTEGER         IOS         ! error status
        INTEGER      :: IREC = 0    ! record counter
        INTEGER         DAYIDX      ! weekday vs. weekend index
        INTEGER         FIPIDX      ! current FIPS index
        INTEGER         HOURIDX     ! current hour index
        INTEGER         SCCIDX      ! current SCC index
        INTEGER         NLINES      ! number of lines
        INTEGER         CNTY        ! current FIPS code
        
        REAL            SPDVAL      ! hourly speed value

        LOGICAL      :: EFLAG = .FALSE.   ! true: error found

        CHARACTER(500)     LINE     ! line buffer
        CHARACTER(300)     MESG     ! message buffer
        CHARACTER(SCCLEN3) SCC      ! current SCC

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

C.............  Parse line into 50 segments
            CALL PARSLINE( LINE, 50, SEGMENT )

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

C.............  Check temperature values
            DO J = 1, 48
                IF( .NOT. CHKREAL( SEGMENT( 2 + J ) ) ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 ) 'ERROR: Bad hourly ' //
     &                'speed value at line', IREC,
     &                'of hourly speed file.'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF
                
                SPDVAL = STR2REAL( ADJUSTR( SEGMENT( 2 + J ) ) )
                
                DAYIDX = 2
                HOURIDX = J
                IF( J > 24 ) THEN
                    DAYIDX = 1
                    HOURIDX = J - 24
                END IF
                
                SPDPRO( FIPIDX, SCCIDX, DAYIDX, HOURIDX ) = SPDVAL
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

