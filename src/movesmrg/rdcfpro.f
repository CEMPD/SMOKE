
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
        USE MODMERGE, ONLY:  NIPPA, EANAM, NMSPC, EMNAM, NSMATV,
     &                       TSVDESC

C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: CFPRO, EXPCFFLAG, REFCFFLAG

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVIFIP, INVIFIP, NINVSCC, INVSCC

C.........  This module is used for reference county information
        USE MODMBSET, ONLY: NREFC, MCREFIDX, NINVC, MCREFSORT

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
        INTEGER         SPCIDX      ! current species index
        INTEGER         SCCIDX      ! current scc index
        INTEGER         NLINES      ! number of lines
        INTEGER         CNTY        ! current FIPS code
        INTEGER         MON         ! current Month
        INTEGER         MXWARN      !  maximum number of warnings
        INTEGER      :: NWARN = 0   !  current number of warnings

        
        REAL            CFVAL       ! control factor value
        REAL            OLDVAL      ! duplicate control factor value

        LOGICAL      :: EFLAG = .FALSE.   ! true: error found
        LOGICAL      :: DUFLAG = .FALSE.   ! true: Duplicate found 
        LOGICAL, ALLOCATABLE :: CFLAG( :,:,:,: )    ! true: duplicate found

        CHARACTER(500)     LINE     ! line buffer
        CHARACTER(300)     MESG     ! message buffer
        CHARACTER(SCCLEN3) SCC      ! current SCC
        CHARACTER(IOVLEN3) POLNAM, SPCNAM  ! current pollutant-species name 

        INTEGER     NLFIPS(NINVIFIP)      ! FIPS matched
        INTEGER     NLSCCS(NINVSCC)       ! SCC matched
        INTEGER     NLPOLS(NIPPA+NMSPC)        ! POL matched
        INTEGER     NLMONS(12)        ! POL matched

        CHARACTER(16) :: PROGNAME = 'RDCFPRO'   ! program name

C***********************************************************************
C   begin body of subroutine RDCFPRO

C.........  Get maximum number of warnings
        MXWARN = ENVINT( WARNSET, ' ', 100, IOS );

C.........  Allocate storage based on number of FIPs and SCCs in inventory
        ALLOCATE( CFPRO( NINVIFIP, NINVSCC, NIPPA+NMSPC, 12 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CFPRO', PROGNAME )
        ALLOCATE( CFLAG( NINVIFIP, NINVSCC, NIPPA+NMSPC, 12 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CFLAG', PROGNAME )
        CFPRO = 1.0   ! array
        CFLAG = .FALSE.   ! array
 
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
            IF( SEGMENT( 1 ) == '' ) SEGMENT( 1 ) = '0'
            IF( STR2INT( SEGMENT( 1 ) ) == 0  ) THEN 

C.....................  State-level is not applicable when REF_CFPRO_YN is set to Y
                IF( REFCFFLAG ) THEN
                    WRITE( MESG, 94010 ) 'WARNING: Skipping line',
     &                  IREC, ' of control factor file because FIPS code '//
     &                  TRIM( SEGMENT(1) ) // ' is not a reference county'
                    CALL M3MESG( MESG )
                ELSE 
                    WRITE( MESG, 94010 )'WARNING: All counties will be '//
     &                  'controlled by zero or blank FIPS entry at line', IREC 
                    CALL M3MESG( MESG )
                END IF

                NFIPS = NINVIFIP
                DO J = 1, NINVIFIP
                    NLFIPS(J) = J
                END DO

            ELSE         ! FIPS is integer value 
                CNTY = STR2INT( SEGMENT( 1 ) )
                L1 = FIND1( CNTY, NINVIFIP, INVIFIP )

                IF( L1 > 0 ) THEN 
                    NFIPS = 1
                    NLFIPS(1) = L1

C.....................  Propagate reference-county-specific control factor to inventory counties.
                    IF( REFCFFLAG ) THEN
                        L2 = FIND1( CNTY, NREFC, MCREFIDX( :,1 ) )
                        IF( L2 < 1 ) THEN
                            WRITE( MESG, 94010 ) 'WARNING: Skipping line',
     &                          IREC, ' of control factor file because FIPS code ',
     &                          CNTY, ' is not a reference county'
                            CALL M3MESG( MESG )
                            CYCLE
                        END IF

C.........................  find ref county and apply CF to ref-inventory counties
                        DO J = 1, NINVC
                            IF( CNTY == MCREFSORT( J,2 ) ) THEN  ! found matched ref county
                                L1 = FIND1( MCREFSORT( J,1 ), NINVIFIP, INVIFIP )
                                IF( L1 > 0 ) THEN
                                    K = K + 1
                                    NLFIPS( K ) = L1
                                END IF
                            END IF
                        END DO
                        NFIPS = K
                    END IF
                ELSE     ! FIPS is country/state  

C.....................  State-level is not applicable when REF_CFPRO_YN is set to Y
                    IF( MOD( CNTY,1000 ) /= 0 ) THEN
                        WRITE( MESG, 94010 ) 'NOTE: Skipping line',
     &                  IREC, ' of control factor file because FIPS code '
     &                   //TRIM(SEGMENT(1))//' is not listed in the inventory'
                        CALL M3MESG( MESG )
                        CYCLE
                    END IF

                    CNTY = CNTY/1000     ! Convert to State ID
                    DO J = 1, NINVIFIP 
                        IF( CNTY == INVIFIP(J)/1000 ) THEN
                            K = K+1
                            NLFIPS(K) = J
                        END IF
                    END DO
                    NFIPS = K

                END IF

            END IF

            IF( NFIPS == 0 ) THEN
                WRITE( MESG, 94010 ) 'NOTE: Skipping line', 
     &            IREC, ' of control factor file because FIPS code '
     &            //TRIM(SEGMENT(1))//' is not listed in the inventory'
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Find SCC in inventory list
            K = 0 

            IF( SEGMENT( 2 ) == ' ' .OR. 
     &          SEGMENT( 2 ) == '0000000000' ) THEN 
                NSCCS = NINVSCC
                DO J = 1, NINVSCC
                    NLSCCS(J) = J
                END DO
            ELSE
                SCC = ADJUSTL( SEGMENT( 2 ) )
                SCCIDX = FINDC( SCC, NINVSCC, INVSCC )
                IF( SCCIDX .LE. 0 ) THEN
                    WRITE( MESG, 94010 ) 'NOTE: Skipping ' //
     &                "line ", IREC, ' of control factor file because SCC ' 
     &                //SCC// ' is not listed in the inventory.'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF
                NSCCS =1
                NLSCCS(1) = SCCIDX
            END IF

C.............  Check pollutant name and mode
            POLNAM = ADJUSTL( SEGMENT( 3 ) )
            K = 0
            IF( SEGMENT( 3 ) == ' ' ) THEN
                NPOLS  = NIPPA+NMSPC
                DO J = 1, NPOLS
                    NLPOLS( J ) =  J
                END DO
            ELSE 
                NPOLS = 0
                IF ( EXPCFFLAG ) THEN
                    POLIDX = INDEX1( POLNAM, NIPPA, EANAM )
                    SPCIDX = INDEX1( POLNAM, NMSPC, EMNAM )
                    NPOLS =1
                    IF( POLIDX > 0 ) THEN
                        NLPOLS( NPOLS ) = POLIDX
                    ELSE IF( SPCIDX > 0 ) THEN
                        NLPOLS( NPOLS ) = NIPPA + SPCIDX
                    ELSE
                        WRITE( MESG, 94010 ) 'NOTE: Skipping line at',
     &                    IREC, ' of control factor file because ' //
     &                    TRIM(POLNAM)// ' is not in the pollutant/species list.'
                        CALL M3MESG( MESG )
                        CYCLE
                    END IF
                ELSE        ! inv poll-specific control factors application
                    POLIDX = INDEX1( POLNAM, NIPPA, EANAM )
                    IF ( POLIDX > 0) THEN
                        NPOLS = NPOLS + 1 
                        NLPOLS( NPOLS ) = POLIDX
                        DO J = 1, NSMATV
                           L1 = INDEX( TSVDESC( J ), SPJOIN )
                           L2 = LEN_TRIM( TSVDESC( J ) )
                           IF( POLNAM == TSVDESC( J )( 1:L1-1 ) ) THEN
                               SPCNAM = TSVDESC( J )( L1+1:L2 )
                               SPCIDX = INDEX1( SPCNAM, NMSPC, EMNAM )
                               NPOLS = NPOLS + 1
                               NLPOLS( NPOLS ) = NIPPA + SPCIDX
                           END IF
                        END DO
                    ELSE
                        WRITE( MESG, 94010 ) 'NOTE: Skipping line at',
     &                    IREC, ' of control factor file because ' //
     &                    TRIM(POLNAM)// ' is not in the pollutant list.'
                        CALL M3MESG( MESG )
                        CYCLE
                    END IF
                END IF
            END IF

C.............  Check month values
            NMONS = 0
            IF( SEGMENT( 4 ) == ' ' ) SEGMENT( 4 ) = '0'
            IF( STR2INT( SEGMENT( 4 ) ) == 0 ) THEN 
                NMONS = 12
                DO J = 1, NMONS
                    NLMONS( J ) = J
                END DO
            ELSE IF ( .NOT. CHKINT( SEGMENT( 4 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Bad month format '
     &               //TRIM(SEGMENT(4))// ' at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            ELSE         ! month is integer value 
                MON = STR2INT( SEGMENT( 4 ) )
                IF( MON < 1  .OR. MON > 12 ) THEN 
                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 ) 'ERROR: Can not process month '
     &                   //TRIM(SEGMENT(4))// ' at line', IREC
                    CALL M3MESG( MESG )
                    CYCLE
                END IF
                NMONS = 1
                NLMONS( 1 ) = MON
            END IF

C.............  check no of fips/scc/poll/mon
            IF( NFIPS < 1 .OR. NSCCS < 1 .OR. NPOLS < 1 .OR. NMONS < 1 )  CYCLE
                
C.............  check and get up control factor values
            IF ( .NOT. CHKREAL( SEGMENT( 5 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Bad contol factor value'//
     &            ' at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF
 
            CFVAL = STR2REAL( SEGMENT( 5 ) )                

            DO L = 1, NFIPS
            DO J = 1, NSCCS
            DO K = 1, NPOLS
            DO M = 1, NMONS
                DUFLAG = CFLAG(NLFIPS(L), NLSCCS(J), NLPOLS(K), NLMONS(M))
                OLDVAL = CFPRO(NLFIPS(L), NLSCCS(J), NLPOLS(K), NLMONS(M))
                IF( DUFLAG .AND. OLDVAL .NE. CFVAL ) THEN
                    OLDVAL = CFPRO(NLFIPS(L), NLSCCS(J), NLPOLS(K), NLMONS(M))
                    WRITE( MESG, 94050 ) 'ERROR: Duplicate entry at line',
     &                  IREC,': Previous factor:', OLDVAL, ' vs New factor:', CFVAL,
     &                  CRLF() // BLANK10 // ' FIPS:', INVIFIP(NLFIPS(L)), 
     &                  ', SCC: ', INVSCC(NLSCCS(J)), ', POLL: ',
     &                  TRIM( EANAM(NLPOLS(K))), ', MONTH: ', NLMONS(M)
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 ) 
               END IF
               CFPRO(NLFIPS(L), NLSCCS(J), NLPOLS(K), NLMONS(M)) = CFVAL
               CFLAG(NLFIPS(L), NLSCCS(J), NLPOLS(K), NLMONS(M)) = .TRUE.
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

94010   FORMAT( 10( A, :, I8, :, 1X ) )
94050   FORMAT(  A,I8,A,F5.2,A,F5.2,A,I6.6,5A,I2 )
        
        END SUBROUTINE RDCFPRO

