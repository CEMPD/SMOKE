       SUBROUTINE SMKHAPSETUP
       
       USE ATHEAP, ONLY: ATHIV, ATHSNAME
       USE MODEMFAC, ONLY: NTOTHAPS, HAPNAMES, PMHAPS
       
       IMPLICIT NONE
       
       INCLUDE 'ATOX5.I'  ! NMAP
       INCLUDE 'IVTYPE.I' ! MAXVEH
       
       INTEGER :: K, L1   ! counters and indexes
       INTEGER :: IOS     ! I/O status
       INTEGER :: IMAP    ! pollutant number
       INTEGER :: IV      ! vehicle type (1 - 28)
       INTEGER :: JXV = 1 ! emission type (exhaust)
       INTEGER :: IRB = 4 ! ratio type (PM)
       
       LOGICAL, SAVE :: INITIAL = .TRUE.  ! true: first time through subroutine
       
       CHARACTER(16) :: TEMPNAME                  ! temporary pollutant name
       CHARACTER(16) :: PROGNAME = 'SMKHAPSETUP'  ! program name
       
C  Begin body of subroutine SMKHAPSETUP

C  Only need to set things up once, so skip if not first time
       IF( .NOT. INITIAL ) RETURN
       INITIAL = .FALSE.

       NTOTHAPS = NMAP

       IF( NTOTHAPS == 0 ) RETURN

       ALLOCATE( PMHAPS( NMAP ), STAT=IOS )
       CALL CHECKMEM( IOS, 'PMHAPS', PROGNAME ) 
       PMHAPS = 0

C  Count total number of pollutants and determine which are PM based
       DO IMAP = 1, NMAP
           DO IV = 1, MAXVEH
               IF( ATHIV( JXV, IRB, IMAP, IV ) ) THEN
                   NTOTHAPS = NTOTHAPS + 1
                   PMHAPS( IMAP ) = NTOTHAPS
                   EXIT
               END IF
           END DO
       END DO
      
       ALLOCATE( HAPNAMES( NTOTHAPS ), STAT=IOS )
       CALL CHECKMEM( IOS, 'HAPNAMES', PROGNAME )
       HAPNAMES = ' '

C  Store names of toxic pollutants
       DO IMAP = 1, NMAP
           
           TEMPNAME = ADJUSTL( ATHSNAME( IMAP ) )
           L1 = LEN_TRIM( TEMPNAME )
           
           IF( PMHAPS( IMAP ) > 0 ) THEN
               IF( L1 > 8 ) THEN
                   HAPNAMES( IMAP ) = TEMPNAME( 1:8 )
               ELSE
                   HAPNAMES( IMAP ) = TEMPNAME( 1:L1 )
               END IF
               TEMPNAME = HAPNAMES( IMAP )
               HAPNAMES( IMAP ) = TRIM( TEMPNAME ) // '25'
               HAPNAMES( PMHAPS( IMAP ) ) = TRIM( TEMPNAME ) // 'PMC'
           ELSE
               IF( L1 > 11 ) THEN
                   HAPNAMES( IMAP ) = TEMPNAME( 1:11 )
               ELSE
                   HAPNAMES( IMAP ) = TEMPNAME( 1:L1 )
               END IF
           END IF
       END DO
       
       END SUBROUTINE SMKHAPSETUP
       
