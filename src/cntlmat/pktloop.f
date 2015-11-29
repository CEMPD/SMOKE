
        SUBROUTINE PKTLOOP( FDEV, PDEV, CDEV, GDEV, LDEV, MDEV, WDEV, 
     &                      CPYEAR, ACTION, ENAME, PKTCNT, PKTBEG, 
     &                      XRFCNT, LPTMP, LCTMP )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine loops through the packets and counts the number
C      of entries or processes the packets, depending on ACTION.
C
C  PRECONDITIONS REQUIRED:
C      Subroutine must be called with ACTION = 'COUNT' before it is called
C      with ACTION = 'PROCESS'
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Started 3/99 by M. Houyoux
C
C************************************************************************
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
C***************************************************************************

C.........  MODULES for public variables
C.........  This module is for cross reference tables
        USE MODXREF, ONLY: INDXTA, CSRCTA, CSCCTA, CMACTA, CISICA, ISPTA,
     &                     MPRNA

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, NIPPA

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'CPKTDAT.EXT'   !  control packet contents

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)   CRLF
        EXTERNAL       CRLF

C...........   SUBROUTINE ARGUMENTS:

        INTEGER     , INTENT(IN) :: FDEV      ! control packets file unit no.
        INTEGER     , INTENT(IN) :: PDEV      ! file unit no. for tmp PROJ file
        INTEGER     , INTENT(IN) :: CDEV      ! file unit no. for tmp CTL file 
        INTEGER     , INTENT(IN) :: GDEV      ! file unit no. for tmp CTG file
        INTEGER     , INTENT(IN) :: LDEV      ! file unit no. for tmp ALW file
        INTEGER     , INTENT(IN) :: MDEV      ! file unit no. for tmp MACT file
        INTEGER     , INTENT(IN) :: WDEV      ! warnings/errors file unit
        INTEGER     , INTENT(IN) :: CPYEAR    ! year to project to
        CHARACTER(*), INTENT(IN) :: ACTION    ! action to take for packets 
        CHARACTER(*), INTENT(IN) :: ENAME     ! inventory file name 
        INTEGER     , INTENT(IN) :: PKTCNT(  NPACKET ) ! count of packet recs
        INTEGER     , INTENT(IN) :: PKTBEG ( NPACKET ) ! 1st line of pkt in file
        INTEGER , INTENT(IN OUT) :: XRFCNT ( NPACKET ) ! count of x-ref recs
        LOGICAL    , INTENT(OUT) :: LPTMP     ! true: projection tmp file written
        LOGICAL    , INTENT(OUT) :: LCTMP     ! true: control tmp file written

C...........   Local arrays allocated to stack
        LOGICAL USEPOL( NIPPA )

C...........   Derived type local variables
        TYPE ( CPACKET ) PKTINFO     ! packet information

C...........   Other local variables
        INTEGER         I, J, K, L      ! counters and indices

        INTEGER         IOS       ! i/o error status
        INTEGER         IREC      ! line number
        INTEGER         IXSIC     ! index of SIC in master SIC list
        INTEGER         IXSCC     ! index of SCC in master SCC list or left-SCC
        INTEGER         JPOL      ! tmp index to master pollutant list
        INTEGER         JT        ! index to control data tables
        INTEGER         JX        ! index to ungrouped control x-ref tables
        INTEGER         XCNT      ! tmp counter for pt srcs x-ref table(s)

        LOGICAL      :: CFLAG  = .FALSE.   ! true: current line is comment
        LOGICAL      :: EFLAG  = .FALSE.   ! error flag
        LOGICAL      :: LTMP   = .FALSE.   ! tmp logical buffer
        LOGICAL      :: LPOLSPEC= .FALSE.  ! true: packet contains usable pol-specific entries
        LOGICAL      :: OFLAG  = .FALSE.   ! true: at least 1 packet was applied
        LOGICAL      :: SKIPPOL= .FALSE.   ! true: pol-spec entries skipped
        LOGICAL      :: SKIPREC= .FALSE.   ! true: packet entries skipped

        CHARACTER(5)       CPOS        ! char pollutant position in EINAM
        CHARACTER(256)     MESG        ! message buffer
        CHARACTER(MACLEN3) MACZERO     ! buffer for zero MACT code
        CHARACTER(SCCLEN3) SCCZERO     ! buffer for zero SCC
        CHARACTER(SICLEN3) SICZERO     ! buffer for zero SIC
        CHARACTER(SICLEN3) CSIC        ! buffer for char SIC
        CHARACTER(IOVLEN3) POLDUM      ! dummy pollutant variable

        CHARACTER(16) :: PROGNAME = 'PKTLOOP' ! program name

C***********************************************************************
C   Begin body of subroutine PKTLOOP

C.........  Set up zero strings for SCC code of zero and SIC code of zero
        SCCZERO = REPEAT( '0', SCCLEN3 )
        SICZERO = REPEAT( '0', SICLEN3 )
        MACZERO = REPEAT( '0', MACLEN3 )

C.........  Loop through packets
        DO K = 1, NPACKET

C.............  Check if packet is present
            IF( PKTCNT( K ) .EQ. 0 ) CYCLE    ! to next iteration

C.............  Reset status of pollutant-specific entries
            LPOLSPEC = .FALSE.

C.............  When cross-reference sizes have already been defined, allocate
C               unsorted x-ref data
            IF( ACTION .EQ. 'PROCESS' ) THEN

                MESG = 'Processing ' // 
     &                 PKTLIST( K )( 1:LEN_TRIM( PKTLIST( K ) ) ) //
     &                 ' packet...'
                CALL M3MSG2( MESG )
                
                J = XRFCNT( K )

C.................  Deallocate memory for ungrouped cross-reference information
                IF( ALLOCATED( INDXTA ) ) THEN
                    DEALLOCATE( INDXTA, ISPTA, MPRNA, CSCCTA, CMACTA, 
     &                          CISICA, CSRCTA )
                END IF

C.................  Allocate memory for ungrouped cross-reference information
                ALLOCATE( INDXTA( J ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INDXTA', PROGNAME )
                ALLOCATE( ISPTA( J ), STAT=IOS )
                CALL CHECKMEM( IOS, 'ISPTA', PROGNAME )
                ALLOCATE( MPRNA( J ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MPRNA', PROGNAME )
                ALLOCATE( CSCCTA( J ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSCCTA', PROGNAME )
                ALLOCATE( CMACTA( J ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CMACTA', PROGNAME )
                ALLOCATE( CISICA( J ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CISICA', PROGNAME )
                ALLOCATE( CSRCTA( J ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSRCTA', PROGNAME )

            ELSE IF( ACTION .EQ. 'COUNT' ) THEN

                MESG = 'Counting entries in ' // TRIM( PKTLIST( K ) ) //
     &                 ' packet...'
                CALL M3MSG2( MESG )

            END IF

C.............  Skip to first line of packet
            REWIND( FDEV )
            CALL SKIPL( FDEV, PKTBEG( K ) )

C.............  Initialize flag that indicates which pollutants appear in
C               a packet
            USEPOL = .FALSE.   ! array

C.............  Initialize various counters before following loop
            XCNT = 0            ! count of cross-reference records for this packet
            JX   = 0            ! control x-ref table counter
            JT   = 0            ! control packet table counter
            IREC = PKTBEG( K )
   
C.............  Loop through lines of current packet to read them
            I = 0
            DO

C.................  Exit if we've read all the lines in this packet
                IF( I == PKTCNT( K ) ) EXIT
            
C.................  Read packet information (populates CPKTDAT.EXT common)
                CALL RDPACKET( FDEV, PKTLIST( K ), 
     &                         USEPOL, IREC, PKTINFO, CFLAG, EFLAG )

C.................  Skip comment lines
                IF( CFLAG ) CYCLE
                I = I + 1

C.................  Format SIC to turn blank and negative values to 0000 - ensure
C                   SICs provided as 2-digits get left-justified.
                CALL FLTRNEG( PKTINFO%CSIC )     ! Filter 0 and -9 to blank
                CSIC = ADJUSTL( PKTINFO%CSIC )
                IF ( CSIC( 2:4 ) .EQ. '   ' ) THEN
                    PKTINFO%CSIC = CSIC( 1:1 ) // '000'
                ELSE IF ( CSIC( 3:4 ) .EQ. '  ' ) THEN
                    PKTINFO%CSIC = CSIC( 1:2 ) // '00'
                END IF
                CALL PADZERO( PKTINFO%CSIC )     ! Pad LHS with zeros

C................  Check if invalid SIC is provided in a cross-reference file.
C                  This reporting is part of the requirements for SMOKE enhancements 
C                  for toxics (EPA SMOKE/MPEI project, Task 6)
                IF( ACTION       .EQ. 'COUNT' .AND.
     &              PKTINFO%CSIC .NE. SICZERO .AND.
     &              PKTINFO%TSCC .EQ. SCCZERO       ) THEN

C.....................  Use FLTRXREF to check if SIC is in the inventory; pass 
C                       zero values for SCC, MACT, and pollutant so they don't 
C                       cause this record to be skipped
                    POLDUM = ' '                 
                    CALL FLTRXREF( PKTINFO%CFIP, PKTINFO%CSIC, SCCZERO, 
     &                             POLDUM, MACZERO, IXSIC, IXSCC, JPOL,
     &                             LTMP, SKIPREC )
     
                    IF( SKIPREC ) THEN
                        WRITE( MESG, 94010 ) 'WARNING: SIC "' //
     &                     TRIM( PKTINFO%CSIC ) // '" is in ' //
     &                     'cross-reference at line', I, 'of ' //
     &                     CRLF() // BLANK10 // TRIM( PKTLIST( K ) ) //
     &                     'packet, but it is not in the inventory.'
                        CALL M3MSG2( MESG )
                    END IF
                END IF

C.................  For non-point sectors, if plant field is filled in, then
C                   skip this record.
                IF ( CATEGORY /= 'POINT' .AND. 
     &               PKTINFO%PLT /= ' '        ) CYCLE
          
C.................  For CONTROL or MACT packet entries, check application control flag
                IF( PKTLIST( K ) == 'CONTROL' .OR.
     &              PKTLIST( K ) == 'MACT'         ) THEN
                    IF( PKTINFO%APPFLAG /= 'Y' ) CYCLE
                END IF                    

C.................  Post-process x-ref information to scan for '-9', pad
C                   with zeros, compare SCC version master list, compare
C                   SIC version to master list, and compare pollutant name 
C                   with master list.
                CALL FLTRXREF( PKTINFO%CFIP, PKTINFO%CSIC, 
     &                         SCCZERO, PKTINFO%CPOL, 
     &                         PKTINFO%CMCT, IXSIC, 
     &                         IXSCC, JPOL, LTMP, SKIPREC  )
                IF( SKIPREC ) CYCLE  ! Skip this record

                SKIPPOL = ( SKIPPOL .OR. LTMP )

C.................  Format SCC since this wasn't done in FLTRXREF so that
C                   partial-SCCs would not get filtered out.
                CALL FLTRNEG( PKTINFO%TSCC )     ! Filter 0 and -9 to blank
                CALL PADZERO( PKTINFO%TSCC )     ! Pad LHS with zeros

C.................  If SIC is defined, make sure SCC is not defined
                IF( PKTINFO%CSIC .NE. SICZERO .AND. 
     &              PKTINFO%TSCC .NE. SCCZERO       ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Both SCC and SIC ' //
     &                     'values are given at line', I, CRLF() // 
     &                     BLANK10 // 'of "' // TRIM( PKTLIST( K ) ) //
     &                     '" packet. Only the SCC will be used for ' //
     &                     'this cross-reference entry.'
                    CALL M3MSG2( MESG )
                    PKTINFO%CSIC = ' '
                END IF

                XCNT = XCNT + 1

C.................  Write pollutant sorted position to a character string
                WRITE( CPOS, '(I5)' ) JPOL
                PKTINFO%CPOS = CPOS

C.................  Set flag to indicate which packets have pol/act-specific 
C                   entries.
                IF( JPOL .GT. 0 ) LPOLSPEC = .TRUE. 

C.................  Process packet information
C.................  This routine increments JX and JT and stores the control
C                   data tables and ungrouped control x-ref tables
                IF( ACTION .EQ. 'PROCESS' ) THEN

                    CALL FILLCNTL( K, PKTCNT( K ), XRFCNT( K ), 
     &                             PKTINFO, JPOL, JT, JX )

                END IF

            END DO   ! End loop through lines of packet

C.............  Store memory needs for x-ref info for this packet
            IF( ACTION .EQ. 'COUNT' ) XRFCNT( K ) = XCNT  

            IF( ACTION .EQ. 'PROCESS' ) THEN

C.................  Check for overflow and write other end-of-loop messages
                CALL ERRPKTS( PKTLIST( K ), JT, JX, SKIPPOL,  
     &                        PKTCNT( K ), XRFCNT( K ), EFLAG )

                IF( EFLAG ) CYCLE  ! To next iteration (next packet)
        
C.................  Sort temporal cross-reference entries. Since CCOD was used
C                   in building CSRCTA, and CCOD will equal "0" when the x-ref 
C                   entry is not pollutant-specific, the non-pollutant-specific
C                   entries will always appear first.  This is necessary for the
C                   table-generating subroutines.
                CALL SORTIC( XRFCNT( K ), INDXTA, CSRCTA )

C.................  Group cross-reference information for current packet
                CALL XREFTBL( PKTLIST( K ), XRFCNT( K ) )

C.................  Match controls to sources and pollutants, as needed for 
C                   each packet type
                CALL PROCPKTS( PDEV, CDEV, GDEV, LDEV, MDEV, WDEV, 
     &                         CPYEAR, PKTLIST( K ), ENAME, LPOLSPEC, 
     &                         USEPOL, OFLAG, LPTMP, LCTMP )

            END IF  ! End process section

        END DO      ! End loop through packets

C.........  An error was found while reading one or more of the packets
        IF( EFLAG ) THEN
            MESG = 'Problem reading control packets file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF( ACTION .EQ. 'PROCESS' .AND. .NOT. OFLAG ) THEN

            MESG = 'No packets were applied to inventory! ' //
     &             CRLF() // BLANK10 // 
     &             'Input packets did not match inventory ' //
     &             CRLF() // BLANK10 //
     &             'or improper environment variable settings.'

            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C......... Rewind file
        REWIND( FDEV )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010  FORMAT( 10( A, :, I8, :, 1X ) )
       
        END SUBROUTINE PKTLOOP
