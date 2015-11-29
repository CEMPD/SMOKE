
        SUBROUTINE RDPGROUP( )

C.............  Get the number of lines in the stack groups file
            NGLINES = GETFLINE( GDEV, 'Stack groups' )

C.............  Scan the groups file to determine the number of PinG groups
            DO I = 1, NGLINES

                READ( GDEV, 93500, IOSTAT=IOS, END=999 ) 
     &                GID, CSWITCH1
                IREC = IREC + 1

C.................  Check read error status
                IF( IOS .GT. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'I/O error', IOS, 
     &                 'reading stack groups file at line', IREC
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

                IF( CSWITCH1 .NE. ' ' ) NGROUP = NGROUP + 1

            END DO

            REWIND( GDEV )

C.............  Allocate memory for stack groups based on the number of lines
            NGROUP = MAX( NGROUP, 1 )
            CALL ALLOCATE_GROUP_VARIABLES

            MESG = 'Reading stack split groups file...'
            CALL M3MSG2( MESG )

            IF( CFLAG ) THEN
                MESG = 'NOTE: Converting stack parameters from ' //
     &                 'English to metric units'
                CALL M3MSG2( MESG )
            END IF

C.............  Read stack groups file
            IREC = 0
            J    = 0
            DO I = 1, NGLINES

                READ( GDEV, 93500, IOSTAT=IOS, END=999 ) 
     &                GID, CSWITCH1, NPG, LON, LAT, DM, HT, TK, VE, FL
                IREC = IREC + 1

C.................  Check read error status
                IF( IOS .GT. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'I/O error', IOS, 
     &                 'reading stack groups file at line', IREC
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

C.................  Skip entry if not a plume-in-grid source
                IF( CSWITCH1 .EQ. ' ' ) CYCLE

C.................  Convert longitude to western hemisphere if needed
                IF( WFLAG .AND. LON .GT. 0 ) LON = -LON

C.................  Convert from English to metric, if needed...
                IF( CFLAG ) THEN
                    DM = DM * FT2M
                    HT = HT * FT2M
                    TK = ( TK - 32. ) * FTOC + CTOK
                    VE = VE * FT2M
                    FL = FL * FLWE2M
                END IF

C.................  When flow is not defined, set it with the vel & diam
                IF( FL .LE. 0. .AND. VE .GT. 0. ) THEN
                   FL = VE * PI * ( 0.25 * DM * DM )
                END IF

C.................  Store data
                J = J + 1
                IF( J .LE. NGROUP ) THEN
                    GRPIDX ( J ) = J
                    GRPGIDA( J ) = GID
                    GRPCNT ( J ) = NPG 
                    IF( LON .NE. 0. ) GRPLON ( J ) = LON
                    IF( LON .NE. 0. ) GRPXL  ( J ) = LON
                    IF( LAT .NE. 0. ) GRPLAT ( J ) = LAT
                    IF( LAT .NE. 0. ) GRPYL  ( J ) = LAT
                    IF( DM  .GT. 0. ) GRPDM  ( J ) = DM
                    IF( HT  .GT. 0. ) GRPHT  ( J ) = HT
                    IF( TK  .GT. 0. ) GRPTK  ( J ) = TK
                    IF( VE  .GT. 0. ) GRPVE  ( J ) = VE
                    IF( FL  .GT. 0. ) GRPFL  ( J ) = FL
                END IF

            END DO    ! End loop on input file lines

C.............  Abort if overflow
            IF( J .GT. NGROUP ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &                  'INTERNAL ERROR: Number of stack groups ' //
     &                  'J=', J, 
     &                  'exceeds dimension NGROUP=', NGROUP
                CALL M3MSG2( MESG ) 

            END IF

C.............  Sort stack group information
            CALL SORTI1( NGROUP, GRPIDX, GRPGIDA )

C.............  Store sorted stack groups for lookups in reading stack splits 
C               file
C.............  Note that GRPGID does not match with other GRP arrays
            N = 0
            PGID = -9
            DO I = 1, NGROUP

                J = GRPIDX( I )
                IF( GRPGIDA( J ) .NE. PGID ) THEN
                    N = N + 1
                    GRPGID( N ) = GRPGIDA( J )
                    PGID = GRPGIDA( J )
                END IF

            END DO

        END SUBROUTINE RDPGROUP
