
        SUBROUTINE ALOCPKTS( FDEV, WDEV, INYEAR, CPYEAR, PKTCNT, PKTBEG,
     &                       XRFCNT )

C***********************************************************************
C  subroutine body starts at line 107
C
C  DESCRIPTION:
C      This subroutine allocates memory for the valid control packets that are
C      in the input file
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Started 3/99 by M. Houyoux
C
C**************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains the control packet data and control matrices
        USE MODCNTRL, ONLY: CUTCTG, FACCTG, FACMACT, FACRACT, 
     &                      ICTLEQUP, ICTLSIC, FACCEFF, FACREFF, 
     &                      FACRLPN, IALWSIC, FACALW, EMCAPALW,
     &                      EMREPALW, IREASIC, EMREPREA, PRJFCREA,
     &                      MKTPNREA, CSCCREA, CSPFREA, IPRJSIC,
     &                      PRJFC, IEMSSIC, BASCEFF, BASREFF, BASRLPN, 
     &                      EMSCEFF, EMSREFF, EMSRLPN, EMSPTCF, EMSTOTL,
     &                      CTLRPLC, MACEXEFF, MACNWEFF, MACNWFRC,
     &                      CMACSRCTYP

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY:

C.........  This module contains the information about the source category
        USE MODINFO, ONLY:

        IMPLICIT NONE
        
C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS:
        CHARACTER*2   CRLF
        EXTERNAL      CRLF

C...........   SUBROUTINE ARGUMENTS:

        INTEGER     , INTENT (IN) :: FDEV      ! in file unit number
        INTEGER     , INTENT (IN) :: WDEV      ! errors/warning file
        INTEGER     , INTENT (IN) :: INYEAR    ! year to project from 
        INTEGER     , INTENT(OUT) :: CPYEAR    ! year to project to
        INTEGER     , INTENT(OUT) :: PKTCNT( NPACKET ) ! count of packet recs
        INTEGER     , INTENT(OUT) :: PKTBEG( NPACKET ) ! 1st line of pkt in file
        INTEGER     , INTENT(OUT) :: XRFCNT( NPACKET ) ! count of x-ref entries

C...........   Logical names and unit numbers

        INTEGER         PDEV      ! file unit no. for tmp PROJ file
        INTEGER         CDEV      ! file unit no. for tmp CTL file
        INTEGER         GDEV      ! file unit no. for tmp CTG file
        INTEGER         LDEV      ! file unit no. for tmp ALW file
        INTEGER         MDEV      ! file unit no. for tmp MACT file

C...........   Other local variables

        INTEGER         I, J, K, L, L1, L2      ! counters and indices

        INTEGER         IOS       ! i/o error status
        INTEGER         IREC      ! line number
        INTEGER         SYEAR     ! tmp start year for packet

        LOGICAL      :: EFLAG  = .FALSE.   ! error flag
        LOGICAL      :: INSIDE = .FALSE.   ! true: read loop inside a packet
        LOGICAL      :: VALID  = .FALSE.   ! true: current packet is valid
        LOGICAL      :: YRSPEC = .FALSE.   ! true: packet is year-specific

        CHARACTER*7               ACTION      ! buffer for PKTLOOP action
        CHARACTER*300             LINE        ! read buffer for a line
        CHARACTER*300             MESG        ! message buffer

        CHARACTER*16 :: PROGNAME = 'ALOCPKTS' ! program name

C***********************************************************************
C   Begin body of subroutine ALOCPKTS

C.........  Initialize packet count for all valid packets
        PKTCNT = 0   ! array

C.........  Write status message
        MESG = 'Scanning control/projection packets input file...'
        CALL M3MSG2( MESG )

C.........  Loop through control packets file and scan for packets. For each
C           packet type, count the number of entries.
C.........  For the projection packet, only pay attention to those that
C           have the same INYEAR.
C.........  Multiple packets of the same type are not permitted.  For
C           projections, this means that multiple packets for the same INYEAR 
C           are not permitted.

        IREC = 0
        DO                ! head of FDEV read loop

            READ( FDEV, 93000, END = 101, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF( IOS .GT. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'ERROR: I/O error', IOS,
     &              'reading control packets file at line', IREC
                CALL M3MESG( MESG )
                CYCLE

            END IF

            IF( LINE .EQ. ' ' ) CYCLE  ! Skip blank lines

C.............  If inside packet...
            IF( INSIDE ) THEN

                J = INDEX( LINE, '!' )
                IF ( J .LE. 0 ) J = LEN_TRIM( LINE )
                I = INDEX( LINE( 1:J ), '/END/' )
                J = INDEX( LINE( 1:J ), '/'     )

C.................  Check for /END/ of packet
                IF( I .GT. 0 ) THEN
                    INSIDE = .FALSE.

C.................  Encountered a new packet or a bad /END/
                ELSE IF( I .LE. 0 .AND. J .GT. 0 ) THEN 
                    WRITE( MESG,94010 ) 'Problem at line', IREC, 
     &                     'of control packets file.' // CRLF() // 
     &                     BLANK10 // 'Encountered a "/" before /END/.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.................  For valid packets, count records
                ELSEIF( VALID ) THEN

                    PKTCNT( K ) = PKTCNT( K ) + 1   ! Increment for packet K

                END IF

C.............  If outside packet, look for next packet header
            ELSE

                L  = LEN_TRIM( LINE )
                L1 = INDEX( LINE, '/' )
                L2 = INDEX( LINE( L1+1:L ), '/' ) + 1

C.................  Well-formed packet
                IF( L1. GT. 0 .AND. L2 .GT. L1 + 1 ) THEN

C.....................  Check if this packet is in list, with special treatment
C                       for reactivity and projection packets.  Return position
C                       K in parameter list of valid packet names, the control
C                       packet year (if applicable), and a year-specific flag.
                    CALL CHECK_PACKET( LINE, L1, L2, INYEAR, K, 
     &                                 CPYEAR, YRSPEC )

C.....................  Print message because this year-specific packet does
C                       not apply to the inventory file
                    IF( K .LE. 0 .AND. YRSPEC ) THEN
                        WRITE( MESG,94010 )
     &                         'WARNING: Packet ' // LINE( L1:L2 ) //
     &                         'does not apply to current' //
     &                         CRLF() // BLANK10 // 'inventory year (',
     &                         INYEAR, '). It will be skipped.'
                        CALL M3MSG2( MESG )
                        INSIDE = .TRUE.
                        VALID = .FALSE.

C.....................  Otherwise when K < 0, packet is not known
                    ELSE IF( K .LE. 0 ) THEN
                        MESG = 'WARNING: Packet ' // LINE( L1:L2 ) // 
     &                         ' is not recognized. It will be skipped.'
                        CALL M3MSG2( MESG )
                        INSIDE = .TRUE.
                        VALID = .FALSE.

C.....................  When K is found, make sure it's the first time
                    ELSE IF( PKTCNT( K ) .GT. 0 ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 )
     &                         'ERROR: Second packet '// LINE( L1:L2 )//
     &                         'found at line', IREC, 'of control file.'
                        CALL M3MSG2( MESG )
                        INSIDE = .TRUE.
                        VALID = .FALSE.

C.....................  Otherwise, store the line number this packet starts on
                    ELSE
                        INSIDE = .TRUE.
                        VALID = .TRUE.
                        PKTBEG( K ) = IREC

                    END IF

C.................  Badly-formed packet
                ELSE IF( L1 .GT. 0 .AND. L2 .LE. L1 ) THEN

                    WRITE( MESG,94010 ) 
     &                     'WARNING: Badly formed packet at line',
     &                     IREC, 'in control packets file.'
                    INSIDE = .TRUE.
                    VALID = .FALSE.
            
                END IF  ! End well-formed packet or not

            END IF      ! End inside a packet or not

        END DO          ! End loop to scan and count records

101     CONTINUE        ! Exit from read loop 

C.........  Abort of there was an error
        IF( EFLAG ) THEN
            MESG = 'Problem reading control packets file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Allocate memory for control packet inputs...
C.........  CTG packet
        J = PKTCNT( 1 ) 
        ALLOCATE( CUTCTG( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CUTCTG', PROGNAME )
        ALLOCATE( FACCTG( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FACCTG', PROGNAME )
        ALLOCATE( FACMACT( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FACMACT', PROGNAME )
        ALLOCATE( FACRACT( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FACRACT', PROGNAME )
        CUTCTG  = 0.  ! array
        FACCTG  = 0.  ! array
        FACMACT = 0.  ! array
        FACRACT = 0.  ! array

C.........  CONTROL packet
        J = PKTCNT( 2 ) 
        ALLOCATE( ICTLEQUP( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTLEQUP', PROGNAME )
        ALLOCATE( ICTLSIC( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTLSIC', PROGNAME )
        ALLOCATE( FACCEFF( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FACCEFF', PROGNAME )
        ALLOCATE( FACREFF( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FACREFF', PROGNAME )
        ALLOCATE( FACRLPN( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FACRLPN', PROGNAME )
        ALLOCATE( CTLRPLC( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CTLRPLC', PROGNAME )
        ICTLEQUP = 0   ! array
        ICTLSIC  = 0   ! array
        FACCEFF  = 0.  ! array
        FACREFF  = 0.  ! array
        FACRLPN  = 0.  ! array
        CTLRPLC = .FALSE.  ! array

C.........  ALLOWABLE packet
        J = PKTCNT( 3 ) 
        ALLOCATE( IALWSIC( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IALWSIC', PROGNAME )
        ALLOCATE( FACALW( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FACALW', PROGNAME )
        ALLOCATE( EMCAPALW( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMCAPALW', PROGNAME )
        ALLOCATE( EMREPALW( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMREPALW', PROGNAME )
        IALWSIC  = 0  ! array
        FACALW   = 0. ! array
        EMCAPALW = 0. ! array
        EMREPALW = 0. ! array

C.........  REACTIVITY packet
        J = PKTCNT( 5 ) 
        ALLOCATE( IREASIC( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IREASIC', PROGNAME )
        ALLOCATE( EMREPREA( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMREPREA', PROGNAME )
        ALLOCATE( PRJFCREA( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRJFCREA', PROGNAME )
        ALLOCATE( MKTPNREA( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MKTPNREA', PROGNAME )
        ALLOCATE( CSCCREA( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSCCREA', PROGNAME )
        ALLOCATE( CSPFREA( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPFREA', PROGNAME )
        IREASIC  = 0   ! array
        EMREPREA = 0.  ! array
        PRJFCREA = 0.  ! array
        MKTPNREA = 0.  ! array
        CSCCREA  = ' ' ! array
        CSPFREA  = ' ' ! array

C.........  PROJECTION packet
        J = PKTCNT( 6 )
        ALLOCATE( IPRJSIC( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IPRJSIC', PROGNAME )
        ALLOCATE( PRJFC( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRJFC', PROGNAME )
        IPRJSIC = 0   ! array
        PRJFC   = 0.  ! array

C.........  EMS-95 CONTROL packet
        J = PKTCNT( 7 ) 
        ALLOCATE( IEMSSIC( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICTLSIC', PROGNAME )
        ALLOCATE( BASCEFF( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BASCEFF', PROGNAME )
        ALLOCATE( BASREFF( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BASREFF', PROGNAME )
        ALLOCATE( BASRLPN( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BASRLPN', PROGNAME )
        ALLOCATE( EMSCEFF( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMSCEFF', PROGNAME )
        ALLOCATE( EMSREFF( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMSREFF', PROGNAME )
        ALLOCATE( EMSRLPN( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMSRLPN', PROGNAME )
        ALLOCATE( EMSPTCF( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMSPTCF', PROGNAME )
        ALLOCATE( EMSTOTL( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMSTOTL', PROGNAME )
        IEMSSIC  = 0    ! array
        BASCEFF  = 0.   ! array
        BASRLPN  = 100. ! array
        EMSCEFF  = 0.   ! array
        EMSRLPN  = 100. ! array
        EMSPTCF  = 1.   ! array - initialize because it will be used
        EMSTOTL  = 0.   ! array - initialize because it will be checked

C.........  MACT packet
        J = PKTCNT( 8 )
        ALLOCATE( CMACSRCTYP( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CMACSRCTYP', PROGNAME )
        ALLOCATE( MACEXEFF( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MACEXEFF', PROGNAME )
        ALLOCATE( MACNWEFF( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MACNWEFF', PROGNAME )
        ALLOCATE( MACNWFRC( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MACNWFRC', PROGNAME )
        CMACSRCTYP = '00'  ! array
        MACEXEFF = 0.  ! array
        MACNWEFF = 0.  ! array
        MACNWFRC = 0.  ! array

C.........  Make sure that at least one packet is defined
        J = 0
        DO I = 1, NPACKET
            J = J + PKTCNT( I ) 
        END DO 

        IF( J .EQ. 0 ) THEN
            MESG = 'No valid packets found in control packets file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Now determine memory needed the cross-referencing portion of the 
C           packets for each packet type.
C.........  For point sources, the memory needed is not necessarily equal to the
C           number of lines in the packet, because the SIC codes need to be 
C           expanded to SCC codes.  So, for point sources process the packets to
C           find SIC NE 0 and SCC EQ 0, and expand the memory requirements
C           accordingly.

        ACTION = 'COUNT'
        CALL PKTLOOP( FDEV, PDEV, CDEV, GDEV, LDEV, MDEV, WDEV, 
     &                CPYEAR, ACTION, BLANK5, PKTCNT, PKTBEG, 
     &                XRFCNT )

C.........  Rewind file
        REWIND( FDEV )
       
        RETURN
       
C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010  FORMAT( 10( A, :, I8, :, 1X ) )
        
C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram checks the packet header(s) to find
C               valid packets in the input file.
            SUBROUTINE CHECK_PACKET( LINE, L1, L2, INYEAR, PKTIDX, 
     &                               OUTYEAR, YFLAG )

C.............  Subprogram arguments
            CHARACTER(*), INTENT (IN) :: LINE     ! local size for dimensioning
            INTEGER     , INTENT (IN) :: L1       ! start of packet in LINE
            INTEGER     , INTENT (IN) :: L2       ! end of packet in LINE
            INTEGER     , INTENT (IN) :: INYEAR   ! year for valid packets
            INTEGER     , INTENT(OUT) :: PKTIDX   ! packet index in master list
            INTEGER     , INTENT(OUT) :: OUTYEAR  ! year generated by packet(s)
            LOGICAL     , INTENT(OUT) :: YFLAG    ! flag for year-specific pkt

C.............  Local variables
            INTEGER    I, J, L, S1, S2    ! indices
            INTEGER    IOS             ! i/o status

            INTEGER, SAVE :: SAVYEAR = 0  ! saved output year

            LOGICAL, SAVE :: SFLAG     ! true: there is a saved output year

            CHARACTER*40  :: BUFFER    ! packet buffer
            CHARACTER*300 :: MESG      ! message buffer

C----------------------------------------------------------------------

C.............  Initialize outputs
            OUTYEAR = 0
            YFLAG   = .FALSE.

C.............  If spaces are between slashes, reset comparison accordingly
            S1 = L1 + 1
            S2 = L2 - 1
            BUFFER = LINE( S1:S2 )
            J = INDEX( BUFFER, ' ' )
            IF( J .GT. 0 ) BUFFER = BUFFER( 1:J-1 )

C.............  Find packet in LINE in list of packets
            PKTIDX = 0
            DO I = 1, NPACKET

                L = LEN_TRIM( PKTLIST( I ) )
                IF( BUFFER .EQ. PKTLIST( I )( 1:L ) ) THEN
                    PKTIDX = I
                    EXIT
                END IF

            END DO

C.............  If packet is not recognized, return
            IF( PKTIDX .LE. 0 ) RETURN

C.............  Otherwise, continue...

C.............  Check for year information...
            L = L1 + LEN_TRIM( PKTLIST( PKTIDX ) ) + 1
            IF( L2 .GT. L ) THEN

                READ( LINE( L+1:L2 ), *, IOSTAT=IOS ) SYEAR, OUTYEAR

C.................  Misformated year information
                IF( IOS .GT. 0 ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                  'ERROR: I/O error', IOS,
     &                  'reading packet years at line', IREC  ! IREC from main
                    CALL M3MESG( MESG )
                    PKTIDX = -1
                    RETURN

C.................  Start year of packet does not apply (warning)
                ELSE IF( SYEAR .NE. INYEAR ) THEN

C NOTE: It is not this simple, because there might be multiple inventory years
C      in the input inventory.  Instead of comparing to the base year, compare
C      to a list of base years generated in GETSINFO.  The projection packet
C      data tables should then contain the start year.

                    PKTIDX = -1
                    YFLAG   = .TRUE.
                    RETURN

C.................  Start and end year are the same... (error)
                ELSE IF( SYEAR .EQ. OUTYEAR ) THEN
                    EFLAG = .TRUE.      ! EFLAG from main
                    WRITE( MESG,94010 )
     &                  'ERROR: I/O error', IOS,
     &                  'reading packet years at line', IREC  ! IREC from main
                    CALL M3MESG( MESG )
                    PKTIDX = -1
                    RETURN

C.................  Start year is fine, but output year disagrees with other pkt
                ELSE IF( SFLAG .AND. OUTYEAR .NE. SAVYEAR ) THEN ! (warning)
                    WRITE( MESG,94010 )
     &                     'WARNING: Output year', OUTYEAR, 'from ' // 
     &                     LINE( S1:S2 ) // ' packet, does not agree '//
     &                     'with previously encountered output year', 
     &                     SAVYEAR, CRLF() // BLANK10 // 
     &                     'Packet will be skipped.' 
                    CALL M3MESG( MESG )
                    PKTIDX = -1
                    RETURN

C.................  Packet is good
                ELSE

                    SFLAG   = .TRUE.    ! year-specific has been encountered
                    SAVYEAR = OUTYEAR   ! store saved year 
                    YFLAG   = .TRUE.

                END IF   ! End checks on year-specific records

            END IF       ! End year-specific processing

            OUTYEAR = SAVYEAR

            RETURN
       
C------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE CHECK_PACKET

       END SUBROUTINE ALOCPKTS
