
        SUBROUTINE RDCFPRO( CFDEV )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C       Reads control factor data from CFPRO file
C
C  PRECONDITIONS REQUIRED:
C       CFDEV must be opened
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     06/12: Created by Dongmei Yang
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
C.........  This module contains the major data structure and control flags
        USE MODMERGE, ONLY:  NIPPA, EANAM

C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: CFPRO

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
        INTEGER       INDEX1
        INTEGER       FIND1
        INTEGER       FINDC
        INTEGER       STR2INT
        INTEGER       ENVINT
        REAL          STR2REAL
        CHARACTER(2)  CRLF
        
        EXTERNAL BLKORCMT, CHKINT, CHKREAL, FIND1, GETFLINE, 
     &           STR2INT, STR2REAL, CRLF, INDEX1, ENVINT

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: CFDEV             ! CFPRO file unit no.

C...........   Local allocatable arrays

C...........   Local arrays
        CHARACTER(20)  SEGMENT( 10 )          ! parsed input line

C...........   Other local variables
        INTEGER         I, J, K, L, M, L1, L2    ! counters and indexes
        INTEGER         IOS         ! error status
        INTEGER      :: IREC = 0    ! record counter
        INTEGER         NFIPS       ! total matched FIPS
        INTEGER         NSCCS       ! total matched SCC
        INTEGER         NPOLS       ! total matched POL 
        INTEGER         NMONS       ! total matched POL 
        INTEGER         POLIDX      ! current POL index
        INTEGER         SCCIDX      ! current POL index
        INTEGER         NLINES      ! number of lines
        INTEGER         CNTY        ! current FIPS code
        INTEGER         MON         ! current Month
        INTEGER         MXWARN      !  maximum number of warnings
        INTEGER      :: NWARN = 0   !  current number of warnings

        
        REAL            CFVAL       ! control factor value
        REAL            OLDVAL      ! duplicate control factor value

        LOGICAL      :: EFLAG = .FALSE.   ! true: error found

        CHARACTER(500)     LINE     ! line buffer
        CHARACTER(300)     MESG     ! message buffer
        CHARACTER(SCCLEN3) SCC      ! current SCC
        CHARACTER(IOVLEN3) POLNAME  ! current pollutant name 
        CHARACTER(IOVLEN3) EPOLNAM  ! current mode+pollutant name 
        CHARACTER(3)       MODNAME  ! current mode name

        INTEGER     NLFIPS(NINVIFIP)      ! FIPS matched
        INTEGER     NLSCCS(NINVSCC)       ! SCC matched
        INTEGER     NLPOLS(NIPPA)        ! POL matched
        INTEGER     NLMONS(12)        ! POL matched

        CHARACTER(16) :: PROGNAME = 'RDCFPRO'   ! program name

C***********************************************************************
C   begin body of subroutine RDCFPRO

C.........  Get maximum number of warnings
        MXWARN = ENVINT( WARNSET, ' ', 100, IOS );

C.........  Allocate storage based on number of FIPs and SCCs in inventory
        ALLOCATE( CFPRO( NINVIFIP, NINVSCC, NIPPA, 12 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CFPRO', PROGNAME )
        CFPRO = 1.0   ! array
        
C.........  Get the number of lines in the file
        NLINES = GETFLINE( CFDEV, 'CFPRO file' )

C.........  Read through file and store hourly data
        DO I = 1, NLINES
            NFIPS = 0
            NSCCS = 0
            NPOLS = 0
            NMONS = 0
        
            READ( CFDEV, 93000, END=999, IOSTAT=IOS ) LINE
            
            IREC = IREC + 1
            
            IF( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'I/O error', IOS,
     &            'reading control factor file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Skip blank or comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Parse line into fields
            CALL PARSLINE( LINE, 6, SEGMENT )

C.............  Convert FIP to integer
            K = 0
            IF( STR2INT( SEGMENT( 1 ) ) == 0  ) THEN 
                NFIPS = NINVIFIP
                DO J = 1, NINVIFIP
                    NLFIPS(J) = J
                END DO
            ELSE IF ( .NOT. CHKINT( SEGMENT( 1 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94020 ) 'ERROR: Bad FIPS code ', TRIM(SEGMENT( 1 )) //
     &            ' at line ', IREC, ' of control factor file.'
                CALL M3MESG( MESG )
                CYCLE
            ELSE         ! FIPS is integer value 
                L2 = LEN_TRIM(SEGMENT( 1 ))
                CNTY = STR2INT( SEGMENT( 1 ) )
                L1 = FIND1( CNTY, NINVIFIP, INVIFIP )
                IF( L1 > 0 ) THEN 
                    NFIPS = 1
                    NLFIPS(1) = L1
                ELSE     ! FIPS is country/state  
                    CNTY = CNTY/1000
                    DO J =1, NINVIFIP 
                        IF ( CNTY == INVIFIP(J)/1000 ) THEN
                            K = K+1
                            NLFIPS(K) = J
                        END IF
                    END DO
                    NFIPS = K
                END IF
            END IF

            IF( NFIPS == 0 ) THEN
                WRITE( MESG, 94010 ) 'WARNING: Skipping line', 
     &            IREC, ' of control factor file because FIPS code', 
     &            TRIM(SEGMENT( 1 )), ' is not in the inventory.'
                CALL M3WARN(PROGNAME, 0, 0, MESG )
                CYCLE
            END IF

C.............  Find SCC in inventory list
            K = 0 
            IF ( SEGMENT( 2 ) == ' ' ) THEN 
                NSCCS = NINVSCC
                DO J = 1, NINVSCC
                    NLSCCS(J) = J
                END DO
            ELSE
                SCC = ADJUSTL( SEGMENT( 2 ) )
                SCCIDX = FINDC( SCC, NINVSCC, INVSCC )
                IF( SCCIDX .LE. 0 ) THEN
                    WRITE( MESG, 94010 ) 'WARNING: Skipping ' //
     &                "line ", IREC, ' of control factor file because SCC ', 
     &                 SCC, ' is not in the inventory.'
                    CALL M3WARN(PROGNAME, 0, 0, MESG )
                    CYCLE
                END IF
                NSCCS =1
                NLSCCS(1) = SCCIDX
            END IF

C.............  Check pollutant name and mode
            POLNAME = ADJUSTL( SEGMENT( 3 ) )
            MODNAME = ADJUSTL( SEGMENT( 4 ) )
            K = 0
            IF( SEGMENT( 3 ) == ' ' .AND. SEGMENT( 4 ) == ' ' ) THEN
                NPOLS  = NIPPA
                DO J = 1, NIPPA
                    NLPOLS(J) =  J
                END DO
            ELSE IF( SEGMENT( 3 ) == ' ' ) THEN
                L1 = LEN_TRIM(SEGMENT( 4 ))
                DO J = 1, NIPPA
                    IF ( MODNAME == EANAM(J)(1:L1) ) THEN
                        K = K +1
                        NLPOLS(K) = J
                    END IF
                END DO
                NPOLS = K
                IF (  NPOLS .EQ. 0 ) THEN 
                    WRITE( MESG, 94010 ) 'WARNING: Skipping ' //
     &                "line ", IREC, ' of control factor file because ',
     &                TRIM(MODNAME), ' is not in the mode list.'
                    CALL M3WARN(PROGNAME, 0, 0, MESG )
                    CYCLE
                END IF
            ELSE IF( SEGMENT( 4 ) .EQ. ' ' ) THEN
                DO J = 1, NIPPA
                    L1 = INDEX(EANAM(J), ETJOIN)
                    IF ( POLNAME == EANAM(J)(L1+2:) ) THEN
                         K = K +1
                         NLPOLS(K) = J
                    END IF
                END DO
                NPOLS = K
                IF (  NPOLS .EQ. 0 ) THEN 
                    WRITE( MESG, 94010 ) 'WARNING: Skipping ' //
     &                "line ", IREC, ' of control factor file because ',
     &                TRIM(POLNAME), ' is not in the pollutant list.'
                    CALL M3WARN(PROGNAME, 0, 0, MESG )
                    CYCLE
                END IF
            ELSE 
                L1 =  LEN_TRIM(POLNAME)
                L2 =  LEN_TRIM(MODNAME)
                EPOLNAM = MODNAME(1:L2)//ETJOIN//POLNAME(1:L1)
                POLIDX =  INDEX1( EPOLNAM, NIPPA, EANAM )
                IF ( POLIDX > 0) THEN
                    NPOLS = 1
                    NLPOLS(1) = POLIDX
                END IF
                IF (  NPOLS .EQ. 0 ) THEN 
                    WRITE( MESG, 94010 ) 'WARNING: Skipping ' //
     &                "line ", IREC, ' of control factor file because ',
     &                EPOLNAM, ' is not in the pollutant list.'
                    CALL M3WARN(PROGNAME, 0, 0, MESG )
                    CYCLE
                END IF
            END IF


C.............  Check month values
            NMONS = 0
            IF( STR2INT( SEGMENT( 5 ) ) ) THEN 
                NMONS = 12
                DO J = 1, NMONS
                    NLMONS(J) = J
                END DO
            ELSE IF ( .NOT. CHKINT( SEGMENT( 6 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94020 ) 'ERROR: Bad month code ', TRIM(SEGMENT( 5 ))//
     &            ' at line', IREC, ' of control factor file.'
                CALL M3MESG( MESG )
                CYCLE
            ELSE         ! month is integer value 
                MON = STR2INT( SEGMENT( 5 ) )
                IF( MON < 0  .OR. MON > 12 ) THEN 
                    EFLAG = .TRUE.
                    WRITE( MESG, 94020 ) 'ERROR: Bad month code ', TRIM(SEGMENT( 5 ))//
     &                ' at line', IREC, ' of control factor file.'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF
                NMONS = 1
                NLMONS(1) = MON
            END IF
                
C.............  check and get up control factor values
            IF ( .NOT. CHKREAL( SEGMENT( 6 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Bad contol factor value ' //
     &            'at line', IREC, 'of control factor file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF
 
            CFVAL = STR2REAL( SEGMENT( 6 ) )                
            DO L = 1, NFIPS
            DO J = 1, NSCCS
            DO K = 1, NPOLS
            DO M = 1, NMONS
                OLDVAL = CFPRO(NLFIPS(L), NLSCCS(J), NLPOLS(K), NLMONS(M))
                IF ( OLDVAL < 1.0 ) THEN
                IF( NWARN < MXWARN ) THEN
                    WRITE( MESG, 94050 ) 'WARNING: Duplicate entry: Overwriting '//
     &                  'factor from ', OLDVAL, ' to ', CFVAL,
     &                  CRLF() // BLANK10 // ' FIPS:', INVIFIP(NLFIPS(L)), 
     &                  ', SCC: ', INVSCC(NLSCCS(J)), ', POLL: ', EANAM(NLPOLS(K)),
     &                  ', MONTH: ', NLMONS(M)
                    CALL M3WARN(PROGNAME, 0, 0, MESG )
                    NWARN = NWARN + 1
               END IF
               END IF
               CFPRO(NLFIPS(L), NLSCCS(J), NLPOLS(K), NLMONS(M)) = CFVAL
            END DO
            END DO
            END DO
            END DO

        END DO  

        CLOSE( CFDEV )
        
        IF( EFLAG ) THEN
            MESG = 'Problem found in control factor file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
C       MESG = "test read control factor file "
C       CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        
        RETURN

999     MESG = 'End of file'
        MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of CFPRO' // CRLF() // BLANK5 //
     &         'input file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )  
      
C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, I8, 1X, A, A,: ) )
94020   FORMAT( 10( 2A, 1X, I8, A ,:) )
94050   FORMAT(  A,:,F5.2,:,A,:,F5.2,:,A,:,I6,5A,A,A, I2 )
        
        END SUBROUTINE RDCFPRO

