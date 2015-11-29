        SUBROUTINE RDPSPLIT( )

C.............  Allocate memory for reading stack splits file
            NSLINES = GETFLINE( TDEV, 'Stack splits' )

            ALLOCATE( SPTINDX( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTINDX', PROGNAME )
            ALLOCATE( SPTGIDA( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTGIDA', PROGNAME )
            ALLOCATE( SPTMMSA( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTMMSA', PROGNAME )
            ALLOCATE( SPTMPSA( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTMPSA', PROGNAME )
            ALLOCATE( SPTLON( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTLON', PROGNAME )
            ALLOCATE( SPTLAT( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTLAT', PROGNAME )
            ALLOCATE( SPTDM( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTDM', PROGNAME )
            ALLOCATE( SPTHT( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTHT', PROGNAME )
            ALLOCATE( SPTTK( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTTK', PROGNAME )
            ALLOCATE( SPTVE( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTVE', PROGNAME )
            ALLOCATE( SPTFL( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTFL', PROGNAME )
            ALLOCATE( SPTCSRCA( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTCSRCA', PROGNAME )
            ALLOCATE( SPTCSRC( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTCSRC', PROGNAME )
            ALLOCATE( FOUND  ( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'FOUND', PROGNAME )

C.............  Initialize
            SPTINDX = 1

C.............  Initialize status of PSPLIT entries found in inventory
            FOUND = .FALSE.    ! array 
        
            MESG = 'Reading stack splits file...'
            CALL M3MSG2( MESG )

C.............  Read stack splits file
            IREC = 0
            DO I = 1, NSLINES

                READ( TDEV, 93550, IOSTAT=IOS, END=999 ) 
     &                GID, CSWITCH1, CSWITCH2, COID, STID, CYID, PLT,
     &                CHAR1, CHAR2, LON, LAT, DM, HT, TK, VE, FL
                IREC = IREC + 1

C.................  Check read error status
                IF( IOS .GT. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'I/O error', IOS, 
     &                 'reading stack splits file at line', IREC
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

C.................  Store contents of splits file entry
                SPTINDX( I ) = I
                SPTGIDA( I ) = GID
                SPTMMSA( I ) = ( CSWITCH1 .NE. ' ' )
                SPTMPSA( I ) = ( CSWITCH2 .NE. ' ' )
                SPTLON ( I ) = LON
                SPTLAT ( I ) = LAT
                SPTDM  ( I ) = DM
                SPTHT  ( I ) = HT
                SPTTK  ( I ) = TK
                SPTVE  ( I ) = VE
                SPTFL  ( I ) = FL

                FIP = COID * 100000 + STID * 1000 + CYID
                WRITE( CFIP, FMTFIP ) FIP

                CSRC = ' '
                CALL BLDCSRC( CFIP, PLT, CHAR1, CHAR2, CHRBLNK3,
     &                        CHRBLNK3, CHRBLNK3, POLBLNK3, CSRC )

                SPTCSRCA( I ) = CSRC

            END DO    ! End loop on input file lines

            MESG = 'Processing splits data with inventory...'
            CALL M3MSG2( MESG )

C.............  Sort splits file source characteristics
            CALL SORTIC( NSLINES, SPTINDX, SPTCSRCA ) 

C.............  Store sorted splits file source characteristics for searching
            DO I = 1, NSLINES
                J = SPTINDX( I )
                SPTCSRC( I ) = SPTCSRCA( J )
            END DO

        END SUBROUTINE RDPSPLIT
