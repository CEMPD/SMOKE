
        SUBROUTINE PKTLOOP( FDEV, NSRC, CPYEAR, NPACKET, ACTION, ENAME, 
     &                      PKTCNT, PKTBEG, PKTLIST, XRFCNT )

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
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C.........  MODULES for public variables
C.........  This module is for cross reference tables
        USE MODXREF

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'CPKTDAT.EXT'   !  control packet contents

C...........   SUBROUTINE ARGUMENTS:

        INTEGER     , INTENT(IN) :: FDEV      ! control packets file unit no.
        INTEGER     , INTENT(IN) :: NSRC      ! number of sources
        INTEGER     , INTENT(IN) :: CPYEAR    ! year to project to
        INTEGER     , INTENT(IN) :: NPACKET   ! number of valid packets
        CHARACTER(*), INTENT(IN) :: ACTION    ! action to take for packets 
        CHARACTER(*), INTENT(IN) :: ENAME     ! inventory file name 
        INTEGER     , INTENT(IN) :: PKTCNT(  NPACKET ) ! count of packet recs
        INTEGER     , INTENT(IN) :: PKTBEG ( NPACKET ) ! 1st line of pkt in file
        CHARACTER(*), INTENT(IN) :: PKTLIST( NPACKET ) ! packet names 
        INTEGER , INTENT(IN OUT) :: XRFCNT ( NPACKET ) ! count of x-ref recs

C...........   Derived type local variables
        TYPE ( CPACKET ) PKTINFO     ! packet information

C...........   Other local variables
        INTEGER         I, J, K      ! counters and indices

c temp:
        INTEGER   :: IDUM =0
        INTEGER         IOS       ! i/o error status
        INTEGER         IREC      ! line number
        INTEGER         IXSIC     ! index of SIC in master SIC list
        INTEGER         IXSCC     ! index of SCC in master SCC list or left-SCC
        INTEGER         JPOL      ! tmp index to master pollutant list
        INTEGER         JT        ! index to control data tables
        INTEGER         JX        ! index to ungrouped control x-ref tables
        INTEGER         NEND      ! tmp loop end for expanding SIC-based rec
        INTEGER         NSTART    ! tmp loop start for expanding SIC-based rec
        INTEGER         XCNT      ! tmp counter for pt srcs x-ref table(s)

        LOGICAL      :: EFLAG  = .FALSE.   ! error flag
        LOGICAL      :: EXPAND = .FALSE.   ! true: expand SIC-based rec to SCCs
        LOGICAL      :: LTMP   = .FALSE.   ! tmp logical buffer
        LOGICAL      :: SKIPPOL= .FALSE.   ! true: pol-spec entries skipped
        LOGICAL      :: SKIPREC= .FALSE.   ! true: packet entries skipped

        CHARACTER*5     CPOS               ! char pollutant position in EINAM
        CHARACTER*300   MESG               ! message buffer
        CHARACTER(LEN=SCCLEN3) SCCZERO     ! buffer for zero SCC

        CHARACTER*16 :: PROGNAME = 'PKTLOOP' ! program name

C***********************************************************************
C   Begin body of subroutine PKTLOOP

C.........  Set up zero strings for FIPS code of zero and SCC code of zero
        SCCZERO = REPEAT( '0', SCCLEN3 )

C.........  Loop through packets
        DO K = 1, NPACKET

C.............  Check if packet is present
            IF( PKTCNT( K ) .EQ. 0 ) CYCLE    ! to next iteration

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
                    DEALLOCATE( INDXTA, ISPTA, MPRNA, CSCCTA, CSRCTA )
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
                ALLOCATE( CSRCTA( J ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSRCTA', PROGNAME )

            END IF

C.............  Skip to first line of packet
            REWIND( FDEV )
            CALL SKIPL( FDEV, PKTBEG( K ) )

C.............  Initialize various counters before following loop
            XCNT = 0            ! count of cross-reference records for this packet
            JX   = 0            ! control x-ref table counter
            JT   = 0            ! control packet table counter
            IREC = PKTBEG( K )
   
C.............  Loop through lines of current packet to read them
            DO I = 1, PKTCNT( K )

C.................  Read packet information (populates CPKTDAT.EXT common)
                CALL RDPACKET( FDEV, PKTLIST( K ), IREC, PKTINFO, 
     &                         EFLAG )

C.................  Post-process x-ref information to scan for '-9', pad
C                   with zeros, compare SCC version master list, compare
C                   SIC version to master list, and compare pollutant name 
C                   with master list.
                CALL FLTRXREF( PKTINFO%CFIP, PKTINFO%CSIC, 
     &                         PKTINFO%TSCC, PKTINFO%CPOL, IXSIC, 
     &                         IXSCC, JPOL, LTMP, SKIPREC  )
     
                IF( SKIPREC ) CYCLE  ! Skip this record

                SKIPPOL = ( SKIPPOL .OR. LTMP )

C.................  Initialize settings for no SIC expansion
                EXPAND = .FALSE.
                NSTART = IXSCC     ! Position of TSCC in INVSCC
                NEND   = IXSCC

C.................  When SIC is defined, but SCC is not defined, need to 
C                   count the additional records to insert.
                IF( IXSIC .GT. 0       .AND. 
     &              PKTINFO%TSCC  .EQ. SCCZERO       ) THEN
                
C.....................  Using position of SIC in list of inventory SICs, 
C                       extract start and end position of SCC in INVSCC
                    NSTART = IBEGSIC( IXSIC )
                    NEND   = IENDSIC( IXSIC )
                    EXPAND = .TRUE.

                    XCNT = XCNT + NEND - NSTART + 1

C.................  When SCC is defined, just add one to the count
                ELSE

                    XCNT = XCNT + 1

                END IF

C.................  Write pollutant sorted position to a character string
                WRITE( CPOS, '(I5)' ) JPOL
                PKTINFO%CPOS = CPOS

C.................  Process packet information
C.................  This routine increments JX and JT and stores the control
C                   data tables and ungrouped control x-ref tables
                IF( ACTION .EQ. 'PROCESS' ) THEN

                    CALL FILLCNTL( PKTLIST( K ), PKTCNT( K ),
     &                             XRFCNT( K ), NSTART, NEND,
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
                CALL PROCPKTS( NSRC, CPYEAR, PKTLIST( K ), ENAME )

            END IF  ! End process section

        END DO      ! End loop through packets

C.........  An error was found while reading one or more of the packets
        IF( EFLAG ) THEN
            MESG = 'Problem reading control packets file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C......... Rewind file
        REWIND( FDEV )
       
        END SUBROUTINE PKTLOOP
