
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
C     06/15: Updated by BH Baek
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
        CHARACTER(20)  SEGMENT( 10 )          ! parsed input line

        INTEGER,            ALLOCATABLE :: NLFIPS( : )     ! FIPS matched
        INTEGER,            ALLOCATABLE :: NLSCCS( : )     ! SCC matched
        INTEGER,            ALLOCATABLE :: NLPOLS( : )     ! POL matched
        INTEGER,            ALLOCATABLE :: NLMONS( : )     ! month matched
        INTEGER           , ALLOCATABLE :: INDXCF( : )     ! index for search matched x-ref
        
        CHARACTER(FPSLEN3), ALLOCATABLE :: CHRT02( : )
        CHARACTER(FPSLEN3), ALLOCATABLE :: CHRT03( : )
        CHARACTER(FPSLEN3), ALLOCATABLE :: CHRT04( : )
        CHARACTER(FPSLEN3), ALLOCATABLE :: CHRT05( : )
        CHARACTER(FPSLEN3), ALLOCATABLE :: CHRT06( : )
        CHARACTER(FPSLEN3+CASLEN3+2), ALLOCATABLE :: DUPCHECK( : )   ! store duplicate entries
        CHARACTER(300),     ALLOCATABLE :: CFIPSCCLIST( : )          ! store matched entries

C...........   Other local variables
        INTEGER         I, J, K, L, M, N, NS, NX, NSCC, L1, L2    ! counters and indexes
        INTEGER         IOS         ! error status
        INTEGER      :: IREC = 0    ! record counter
        INTEGER         NFIPS       ! total matched FIPS
        INTEGER         NSCCS       ! total matched SCC
        INTEGER         NPOLS       ! total matched POL 
        INTEGER         POLIDX      ! current POL index
        INTEGER         SPCIDX      ! current species index
        INTEGER         SCCIDX      ! current scc index
        INTEGER         NLINES      ! number of lines
        INTEGER         NXREF       ! max no of matached x-ref entries
        INTEGER         IFIP        ! current FIPS code
        INTEGER         MON, NMONS  ! current Month
        INTEGER         IDUM        ! tmp dummy integer
        INTEGER         MXWARN      ! maximum number of warnings
        INTEGER      :: NWARN = 0   ! current number of warnings
        INTEGER         ISTA, PSTA, NSTA
        INTEGER         NCHRT02, NCHRT03, NCHRT04, NCHRT05, NCHRT06
        INTEGER         FF, F1, F2, F3, F4, F5, F6
        
        REAL            CFVAL       ! control factor value

        LOGICAL, SAVE:: FIRSTIME = .TRUE. ! true: firstime
        LOGICAL      :: EFLAG = .FALSE.   ! true: error found
        LOGICAL      :: PFLAG, SKIPREC    ! true: skip pollutant 

        CHARACTER(300)     LINE     ! line buffer
        CHARACTER(300)     MESG     ! message buffer
        CHARACTER(FPSLEN3+CASLEN3+2) DUPTMP      ! tmp duplicate check buffer
        CHARACTER(FIPLEN3) CFIP, BLKFIP          ! tmp (character) FIPS code : CHRT03
        CHARACTER(SCCLEN3) CSCC, BLKSCC          ! current SCC : CHRT04
        CHARACTER(FPSLEN3) CFIPSCC, BLKFIPSCC    ! tmp FIPS code // SCC : CHRT06
        CHARACTER(IOVLEN3) POLNAM, SPCNAM, CBUF  ! current pollutant-species name 

        CHARACTER(16) :: PROGNAME = 'RDCFPRO'    ! program name

C***********************************************************************
C   begin body of subroutine RDCFPRO

        IF( FIRSTIME ) THEN
C.............  Get maximum number of warnings
            MXWARN = ENVINT( WARNSET, ' ', 100, IOS )

C.............  No of states
            NSTA = 0
            PSTA = 0
            DO I = 1 , NINVIFIP
                ISTA = INT( INVIFIP( I ) / 1000 ) * 1000
                IF( ISTA /= PSTA ) THEN
                    NSTA = NSTA + 1
                    PSTA = ISTA
                END IF
            END DO

C.............  Allocate locate FIPS/SCC hierarchy arrays
            NCHRT02 = NSTA               ! SCC=0,   FIP=state
            NCHRT03 = NINVIFIP           ! SCC=0,   FIP=all
            NCHRT04 = NINVSCC            ! SCC=all, FIP=0
            NCHRT05 = NSTA * NINVSCC     ! SCC=all, FIP=state
            NCHRT06 = NINVIFIP * NINVSCC ! SCC=all, FIP=all

            ALLOCATE( CHRT02( NCHRT02 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT02', PROGNAME )
            ALLOCATE( CHRT03( NCHRT03 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT03', PROGNAME )
            ALLOCATE( CHRT04( NCHRT04 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT04', PROGNAME )
            ALLOCATE( CHRT05( NCHRT05 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT05', PROGNAME )
            ALLOCATE( CHRT06( NCHRT06 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CHRT06', PROGNAME )
            CHRT02 = ''
            CHRT03 = ''
            CHRT04 = ''
            CHRT05 = ''
            CHRT06 = ''

C.............  Build FIPS and SCC hierarchy tables
            N = 0
            BLKFIP = '000000'
            BLKSCC = '0000000000'
            BLKFIPSCC = BLKFIP // BLKSCC
            DO I = 1, NINVSCC

                CHRT04( I ) = BLKFIP // INVSCC( I )
                NSTA = 0
                PSTA = -1 
                DO J = 1, NINVIFIP

                    N = N + 1
                    WRITE( CHRT03(J),'(I6.6,A)' ) INVIFIP(J), BLKSCC
                    WRITE( CHRT06(N),'(I6.6,A)' ) INVIFIP(J), INVSCC(I)
   
                    NS = 0
                    ISTA = INT( INVIFIP( J ) / 1000 ) * 1000
                    IF( ISTA /= PSTA ) THEN
                        NSTA = NSTA + 1
                        WRITE( CHRT02( NSTA ),'(I6.6,A)' ) ISTA, BLKSCC

                        DO K = 1, NINVSCC
                            NS = NS + 1
                            WRITE( CHRT05(NS),'(I6.6,A)' ) ISTA, INVSCC(K)
                        END DO
                        PSTA = ISTA
                    END IF

                END DO    ! FIPS loop

            END DO     ! SCC loop

C.............  Allocate storage based on number of FIPs and SCCs in inventory
            ALLOCATE( CFPRO( NINVIFIP, NINVSCC, NIPPA+NMSPC, 12 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CFPRO', PROGNAME )
            CFPRO = 1.0   ! array
            
C.............  Get the number of lines in the file
            NLINES = GETFLINE( CFDEV, 'CFPRO file' )
            ALLOCATE( CFIPSCCLIST( NLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CFIPSCCLIST', PROGNAME )
            ALLOCATE( INDXCF( NLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INDXCF', PROGNAME )
            ALLOCATE( DUPCHECK( NLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DUPCHECK', PROGNAME )
            ALLOCATE( NLFIPS( NINVIFIP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NLFIPS', PROGNAME )
            ALLOCATE( NLSCCS( NINVSCC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NLSCCS', PROGNAME )
            ALLOCATE( NLPOLS( NIPPA+NMSPC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NLPOLS', PROGNAME )
            ALLOCATE( NLMONS( 12 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NLMONS', PROGNAME )
            CFIPSCCLIST = ''
            DUPCHECK   = ''
            INDXCF = 0
            NLFIPS = 0
            NLSCCS = 0
            NLPOLS = 0
            NLMONS = 0

            FIRSTIME = .FALSE.

        END IF

C.........  Read through file and store hourly data
        NX = 0
        DO I = 1, NLINES

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
            
            CFIP   = TRIM( SEGMENT( 1 ) )
            CSCC   = TRIM( SEGMENT( 2 ) )
            POLNAM = ADJUSTL( SEGMENT( 3 ) )
            MON    = STR2INT( SEGMENT( 4 ) )
            IF( MON < 1 ) MON = 0

            CBUF = ' '
            CALL FLTRXREF( CFIP, CBUF, CSCC, POLNAM, CBUF,
     &                     IDUM, IDUM, IDUM, PFLAG, SKIPREC )

            CFIPSCC = CFIP // CSCC        ! CHRT09

            F2 = FINDC( CFIPSCC, NCHRT02, CHRT02 )  ! SCC=0,   FIP=state
            F3 = FINDC( CFIPSCC, NCHRT03, CHRT03 )  ! SCC=0,   FIP=all
            F4 = FINDC( CFIPSCC, NCHRT04, CHRT04 )  ! SCC=all, FIP=0
            F5 = FINDC( CFIPSCC, NCHRT05, CHRT05 )  ! SCC=all, FIP=state
            F6 = FINDC( CFIPSCC, NCHRT06, CHRT06 )  ! SCC=all, FIP=all

            FF = 0
            IF( F6 > 0 ) FF = 6
            IF( F5 > 0 ) FF = 5
            IF( F4 > 0 ) FF = 4
            IF( F3 > 0 ) FF = 3
            IF( F2 > 0 ) FF = 2
            IF( CFIPSCC == BLKFIPSCC ) FF = 1

            IF( FF < 1 ) CYCLE     ! skip unmatched entries
            NX = NX + 1

            WRITE( DUPTMP,'(A,I2.2,A)' ) CFIP//CSCC//POLNAM, MON
            L = FINDC( DUPTMP, NLINES, DUPCHECK )
            IF( L > 0 ) THEN
                WRITE( MESG,94010 ) 'ERROR: Duplicate entry found at', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ELSE
                DUPCHECK( NX ) = DUPTMP
            END IF    

            INDXCF( NX ) = FF
            CFIPSCCLIST( NX ) = TRIM( LINE )

C.............  No of CFPRO entries check
            IF( NX == 0 ) THEN
                MESG = 'ERROR: No applicable entry found from CFPRO input file'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 ) 
            END IF
            
        END DO

        NXREF = NX

C.........  Processing a list of valid CFPRO entries
        DO N = 1, 6     ! six hierarchy

          DO I = 1, NXREF
        
            IF( INDXCF( I ) /= N ) CYCLE   ! process low hierarchy first

            NFIPS = 0
            NSCCS = 0
            NPOLS = 0
       
C.............  Parse line into fields
            CALL PARSLINE( CFIPSCCLIST( I ), 6, SEGMENT )
            
            CFIP   = TRIM( SEGMENT( 1 ) )
            CSCC   = TRIM( SEGMENT( 2 ) )
            POLNAM = ADJUSTL( SEGMENT( 3 ) )

            CBUF = ' '
            CALL FLTRXREF( CFIP, CBUF, CSCC, POLNAM, CBUF,
     &                     IDUM, IDUM, IDUM, PFLAG, SKIPREC )

C.............  Convert FIP to integer
            K = 0
            IFIP = STR2INT( CFIP )
            IF( IFIP == 0  ) THEN 

C.....................  State-level is not applicable when REF_CFPRO_YN is set to Y
                IF( REFCFFLAG ) THEN
                    MESG = 'WARNING: FIPS code '//CFIP//' is not a reference county'
                    CALL M3MESG( MESG )
                END IF

                NFIPS = NINVIFIP
                DO J = 1, NINVIFIP
                    NLFIPS( J ) = J
                END DO

            ELSE         ! FIPS is integer value 
                L1 = FIND1( IFIP, NINVIFIP, INVIFIP )

                IF( L1 > 0 ) THEN 
                    NFIPS = 1
                    NLFIPS(1) = L1

C.....................  Propagate reference-county-specific control factor to inventory counties.
                    IF( REFCFFLAG ) THEN
                        L2 = FIND1( IFIP, NREFC, MCREFIDX( :,1 ) )
                        IF( L2 < 1 ) THEN
                            MESG = 'WARNING: FIPS code '//CFIP//' is not a reference county'
                            CALL M3MESG( MESG )
                            CYCLE
                        END IF

C.........................  find ref county and apply CF to ref-inventory counties
                        DO J = 1, NINVC
                            IF( IFIP == MCREFSORT( J,2 ) ) THEN  ! found matched ref county
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
                    IF( MOD( IFIP,1000 ) /= 0 ) THEN
                        MESG = 'WARNING: FIPS '//CFIP//' is not an inventory county'
                        CALL M3MESG( MESG )
                        CYCLE
                    END IF

                    IFIP = IFIP/1000     ! Convert to State ID
                    DO J = 1, NINVIFIP 
                        IF( IFIP == INVIFIP(J)/1000 ) THEN
                            K = K + 1
                            NLFIPS( K ) = J
                        END IF
                    END DO
                    NFIPS = K

                END IF

            END IF

            IF( NFIPS == 0 ) CYCLE

C.............  Find SCC in inventory list
            K = 0 
            IF( CSCC == ' ' .OR. CSCC == BLKSCC ) THEN 
                NSCCS = NINVSCC
                DO J = 1, NINVSCC
                    NLSCCS(J) = J
                END DO
            ELSE
                SCCIDX = FINDC( CSCC, NINVSCC, INVSCC )
                IF( SCCIDX < 1 ) THEN
                    MESG = 'WARNING: SCC '//CSCC//' is not in the inventory.'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF
                NSCCS =1
                NLSCCS(1) = SCCIDX
            END IF

C.............  Check pollutant name and mode
            K = 0
            IF( POLNAM == ' ' ) THEN
                NPOLS  = NIPPA+NMSPC
                DO J = 1, NPOLS
                    NLPOLS( J ) =  J
                END DO
            ELSE 
                NPOLS = 0
                IF ( EXPCFFLAG ) THEN
                    POLIDX = INDEX1( POLNAM, NIPPA, EANAM )
                    SPCIDX = INDEX1( POLNAM, NMSPC, EMNAM )
                    NPOLS = 1
                    IF( POLIDX > 0 ) THEN
                        NLPOLS( NPOLS ) = POLIDX
                    ELSE IF( SPCIDX > 0 ) THEN
                        NLPOLS( NPOLS ) = NIPPA + SPCIDX
                    ELSE
                        MESG = 'WARNING: '//TRIM(POLNAM)//' is not found in' //
     &                         ' both inventory pollutant and model species lists.'
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
                        MESG = 'WARNING: '//TRIM(POLNAM)//' is not found in' //
     &                         ' inventory pollutant list.'
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
                MESG = 'ERROR: Bad month format '//TRIM(SEGMENT(4))
                CALL M3MESG( MESG )
                CYCLE
            ELSE         ! month is integer value
                MON = STR2INT( SEGMENT( 4 ) )
                IF( MON < 1  .OR. MON > 12 ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Bad month format '//TRIM(SEGMENT(4))
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
                MESG = 'ERROR: Bad control factor value'// TRIM(SEGMENT(5))
                CALL M3MESG( MESG )
                CYCLE
            END IF
 
            CFVAL = STR2REAL( SEGMENT( 5 ) )

            DO L = 1, NFIPS
              DO J = 1, NSCCS
                DO K = 1, NPOLS
                  DO M = 1, NMONS
                    CFPRO( NLFIPS(L),NLSCCS(J),NLPOLS(K),NLMONS(M) ) = CFVAL
                  END DO
                END DO
              END DO
            END DO

          END DO  

        END DO  

        REWIND( CFDEV )
        
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

