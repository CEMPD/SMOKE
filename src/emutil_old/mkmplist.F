
         PROGRAM MKMPLIST

C***********************************************************************
C  program body starts at line 147
C
C  DESCRIPTION:
C       Program MKMPLIST creates a very basic default MPLIST file, which
C       can be used as a template for building Mobile source control 
C       strategies.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C	FIND1, GETYN, INDEX1, PROMPTFFILE
C
C  REVISION  HISTORY:
C       Prototype  9/96 by MRH
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
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

      IMPLICIT NONE

C...........   INCLUDES:
       
        INCLUDE 'MBDIMS3.EXT'
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        INTEGER         FIND1 
        LOGICAL         GETYN
        INTEGER         INDEX1
        INTEGER         PROMPTFFILE
        
        EXTERNAL  FIND1, GETYN, INDEX1, PROMPTFFILE

C...........   LOCAL VARIABLES and their descriptions:

C.........  FIPs for matching control strategies

        INTEGER     INDEXA ( MXFIP )  ! Index for sorting
        INTEGER     IFIP1  ( MXFIP )  ! Sorted   FIP code unique list
        INTEGER     IFIP2  ( MXFIP )  ! Sorted   FIP code for matching
        INTEGER     IFIP1A ( MXFIP )  ! Unsorted FIP code unique list
        INTEGER     IFIP2A ( MXFIP )  ! Unsorted FIP code for matching
        INTEGER     NFIP              ! Actual number of FIP entries
 
C.........  Variables for reading in OFNETSPD file

        REAL        SPDA  ( NRCLAS*MXFIP )  ! Unsorted speeds

        INTEGER     IFIPA ( NRCLAS*MXFIP )  ! Unsorted FIP 
        INTEGER     IRCLSA( NRCLAS*MXFIP )  ! Unsorted road class
        INTEGER     ISPIDA( NRCLAS*MXFIP )  ! Unsorted speeds indices to unique
        INTEGER     MFIPA ( NRCLAS*MXFIP )  ! Unsorted matching FIP
        INTEGER     SPDIDX( NRCLAS*MXFIP )  ! Sorting index

        INTEGER     NSPD                    ! Actual number of speed entries

C.........  Variables for creating MPLIST file


        INTEGER     IFIP   ( MXPLST )  ! Sorted Fip code x-ref list
        INTEGER     INDEXB ( MXPLST )  ! Index for re-sorting
        INTEGER     INDEXC ( MXPLST )  ! Index for re-sorting
        INTEGER     IRCLAS ( MXPLST )  ! Sorted Road Class x-ref list
        INTEGER     ISPID  ( MXPLST )  ! Sorted speed IDs pointing to SPDLIST
        INTEGER     IPSI   ( MXPLST )  ! Sorted parameter scheme indices
        INTEGER     NPLIST             ! Actual number of entries in MPLIST


C.........  Variables for default MFIPs
        INTEGER     DEFLST ( MXSTA )        ! MFIPs to be used as state default
        INTEGER     DFSTRCL( NRCLAS, MXSTA )! Speed index per default MFIP and RCL
        INTEGER     ICNTX                   ! First counter for default setting
        INTEGER     ICNTY                   ! Second counter for default setting
        INTEGER     FIPX                    ! First FIP for default setting
        INTEGER     FIPY                    ! Second FIP for default setting
        INTEGER     NDEF                    ! Number of MFIPs that are defaults

C.........  Variables for unique speeds list

        REAL        SPDLIST( MXFACS )       ! Unique EMS95 speeds

        INTEGER     SCNT                    ! Actual number of unique speeds
  
C.........  Unit numbers 
        INTEGER     LDEV         ! For IO/API starting unit
        INTEGER     MDEV         ! For MVCYID input file
        INTEGER     PDEV         ! For MPLIST output file
        INTEGER     SDEV         ! For OFNETSPD input file 
        INTEGER     XDEV         ! For MPREF output file

C.........  Other local variables
        REAL        SPEED        ! Temporary speeds
    
        INTEGER     D, I, IOS, IREC, J, K, R  ! Counters and pointers

        INTEGER     CY1, CY2     ! Temporary county IDs
        INTEGER     FIP          ! Temporary FIP ID
        INTEGER     LSID         ! Previous state ID
        INTEGER     LSPID        ! Previous speed index
        INTEGER     LMFIP        ! Previous MFIP
        INTEGER     LPSI         ! Previous PSI 
        INTEGER     MFIP         ! Temporary FIP ID for matching
        INTEGER     PSI          ! Parameter scheme index
        INTEGER     ST, CY, ROAD ! Temporary variables for reading
        INTEGER     ST1, ST2     ! Temporary state IDs
        INTEGER     RCL          ! Temporary road class
        INTEGER     SID          ! State ID
        INTEGER     SPID         ! Temporary speed ID or output index code

        LOGICAL EFLAG                 ! Error flag
        DATA    EFLAG  / .FALSE. /
        LOGICAL SFLAG                 ! For adjusting speeds

        CHARACTER*2   RR
        CHARACTER*256 MESG
       
C***********************************************************************
C   begin body of program MKMPLIST
        
        LDEV = INIT3()
        
        WRITE( *,92000 ) 
     &  ' ',
     &  'Program MKMPLIST ,',
     &  ' ',
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").',
     &  'Optional checking that all species are reported for each ',
     &  'source may be turned on via "setenv RAW_SRC_CHECK Y".',
     &  ' ',
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry. Default responses are indicated',
     &  'in brackets [LIKE THIS].',
     &  ' '
        
        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            WRITE( *,92000 ) ' ', 'Ending program MKMPLIST', ' '
            CALL EXIT( 0 )
        END IF

C.........  Get name of EMS-95 mvcyid.mv file

        MDEV = PROMPTFFILE( 
     &           'Enter logical name for ASCII MVCYID.MV file',
     &           .TRUE., .TRUE., 'MVCYID', 'MKMPLIST' )

C.........  Get name of EMS-95 ofnetspd.mv file
        SDEV = PROMPTFFILE( 
     &           'Enter logical name for ASCII OFNETSPD.MV file',
     &           .TRUE., .TRUE., 'OFNETSPD', 'MKMPLIST' )

        SFLAG = GETYN( 'Adjust speeds to even 2 MPH interval?', .TRUE. )

C.........  Read MVCYID.MV file for which FIPs have same MOBILE5 files

        IREC = 0     !  Counter for records
        I    = 0     !  Counter for valid records

   55   CONTINUE       !  Head of the ADEV-read loop

            READ( MDEV, *, END=66, IOSTAT=IOS ) ST1, CY1, ST2, CY2

            IREC = IREC + 1

C.............  Abort if error reading file
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG, 94010 ) 'Error ', IOS, 
     &               'reading MVCYID.MV file at line', IREC
                CALL M3EXIT( 'MKMPLIST', 0, 0, MESG, 2 )

            ENDIF

            I = I + 1   

            IF( I .LE. MXFIP ) THEN

                INDEXA( I ) = I 

                IFIP1A( I ) = ST1*1000 + CY1
                IFIP2A( I ) = ST2*1000 + CY2

            ENDIF    ! If overflow or not

        GO TO 55     ! To head of mvcyid.mv read loop

   66   CONTINUE
        NFIP = I 

        IF ( NFIP .GT. MXFIP ) THEN

            WRITE( *,92010 )
     &          'Record count                  ', NFIP,
     &          'MXFIP dimensioned record-count', MXFIP 
            CALL M3EXIT ( 'MKMPLIST', 0, 0, 
     &                    'Max record-count exceeded in ' //
     &                    'MVCYID.MV file', 2 )
        END IF

C.........  Sort MVCYID data
        CALL SORTI2( NFIP, INDEXA, IFIP1A, IFIP2A )

C.........  Created sorted MVCYID list

        DO 69 I = 1, NFIP

            IFIP1( I ) = IFIP1A( INDEXA( I ) )
            IFIP2( I ) = IFIP2A( INDEXA( I ) )

   69   CONTINUE

C.........  Read in OFNETSPD.MV file

        IREC = 0     !  Counter for records
        I    = 0     !  Counter for valid records

   77   CONTINUE       !  Head of the ADEV-read loop

            READ( SDEV, 93010, END=88, IOSTAT=IOS )  ! In EMS95 format
     &            ST, CY, ROAD, SPEED

            IREC = IREC + 1

            IF( IOS .NE. 0 ) THEN
                WRITE( MESG, 94010 ) 'Error ', IOS, 
     &               'reading OFNETSPD.MV file at line', IREC
                CALL M3EXIT( 'MKMPLIST', 0, 0, MESG, 2 )

            ELSE
 
                WRITE( RR, '( I2.2 )' ) ROAD
                J = INDEX1( RR, NRCLAS, MARFAC3 )
 
                IF( J .GT. 0 ) THEN

                    RCL = MROADS3( J )
                    FIP = ST*1000 + CY

                ELSE      ! Error, road class not found
 
                    WRITE( MESG, 94010 )
     &                 'Skipping record', IREC,
     &                 ': Bad area/facility combination ' // RR
                    CALL M3MESG( MESG )
 
                    GO TO 77  ! To head of SDEV read loop
 
                ENDIF

            ENDIF

            I = I + 1   

            IF( I .LE. MXFIP * NRCLAS ) THEN

                J = FIND1( FIP, NFIP, IFIP1 )

                IF( J .GT. 0 ) THEN

                    SPDIDX( I ) = I

                    IFIPA ( I ) = FIP
                    IRCLSA( I ) = RCL
                    SPDA  ( I ) = SPEED
                    MFIPA ( I ) = IFIP2( J ) 

                ELSE

                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                     'FIP from speed file:', FIP,
     &                     'is not found in mvcyid data'
                    CALL M3MESG( MESG )

                ENDIF

            ENDIF     ! If overflow or not

        GO TO 77      ! To head of SDEV-read loop

   88   CONTINUE

        NSPD = I

C.........  Check for overflow and drop sources if necessary
        IF ( NSPD .GT. MXFIP * NRCLAS ) THEN

            EFLAG = .TRUE.
            WRITE( *,92010 )
     &          'Record count                ', NSPD,
     &          'MXFIP*NRCLAS dimensioned record-count', MXFIP*NRCLAS
            CALL M3MESG ( 'Max record-count exceeded in ' //
     &                    'MVCYID.MV file' )
        END IF

        IF( EFLAG ) THEN
            CALL M3EXIT( 'MKMPLIST', 0, 0, 
     &                   'Error in speed records', 2 )

        ENDIF

C.........  Round speeds to EMS-95 equivalent if requested
        IF( SFLAG ) THEN

            DO 111 I = 1, NSPD
                SPDA( I ) = NINT( SPDA( I ) / 2. ) * 2.0
  111       CONTINUE

        ENDIF

C.........  Make list of default MFIPs (the ones that are used the most
C.........  in each state).  For a new state, the MFIP with the most 
C.........  FIPs pointing to it is the default.

        NDEF  = 0
        ICNTX = 0 
        ICNTY = 0
        FIPX  = IMISS3
        FIPY  = IMISS3
        LSID  = IMISS3

        DO 122 I = 1, NFIP

            SID  = IFIP1( I ) / 1000

            IF( SID .EQ. LSID ) THEN              ! Same state, so count...

                IF( IFIP2( I ) .EQ. FIPX ) THEN   ! Add 1 for current max

                    ICNTX = ICNTX + 1

                ELSE                              ! Not the max FIP

                    IF( IFIP2( I ) .NE. FIPY ) THEN ! New FIP

                        FIPY  = IFIP2( I )
                        ICNTY = 1

                    ELSE

                        ICNTY = ICNTY + 1           ! Increase counter for new 
                        IF( ICNTY .GT. ICNTX ) THEN ! Replace if new > old
                            FIPX  = FIPY
                            ICNTX = ICNTY
                        ENDIF

                    ENDIF
                ENDIF

            ELSE  ! New state

                IF( I .GT. 1 ) THEN
                    NDEF  = NDEF  + 1

                    IF( NDEF .LE. MXSTA ) THEN
                        DEFLST( NDEF ) = FIPX   ! Store for previous state
                    ENDIF  ! If overflow or not

                ENDIF

                ICNTX = 1               ! initialize for this state
                ICNTY = 1               ! initialize for this state
                FIPX  = IFIP2( I )      ! initialize for this state
                FIPY  = IFIP2( I )      ! initialize for this state
                LSID  = SID

            ENDIF

  122   CONTINUE

C.........  Store for last state in list

        NDEF = NDEF + 1

        IF( NDEF .LE. MXSTA ) THEN
            DEFLST( NDEF ) = FIPX
        ENDIF

        IF( NDEF .GT. MXSTA ) THEN
            WRITE( *,92010 )
     &          'State defaults record count   ', NDEF,
     &          'MXSTA dimensioned record-count', MXSTA
            CALL M3EXIT ( 'MKMPLIST', 0, 0, 
     &                    'Max record-count exceeded', 2 )

        ENDIF

C.........  Initialize DFSTRCL to zero
        DO 128 D = 1, MXSTA
            DO 127 R = 1,NRCLAS

                DFSTRCL( R,D ) = IMISS3

  127       CONTINUE
  128   CONTINUE
 

C.........  Make list of unique speeds (speeds are REAL and not sorted)
C.........  and assign speed list index ID for each speed record.
C.........  Yes, this is a very inefficient algorithm.
        SPDLIST( 1 ) = SPDA( 1 )
        SCNT = 1

        DO 133 I = 1, NSPD

            J = 0 
  129       CONTINUE        

                J = J + 1

                IF( SPDA( I ) .EQ. SPDLIST( J ) ) THEN

                    ISPIDA( I ) = J
                    GO TO 133

                ELSEIF( J .EQ. SCNT ) THEN
 
                    SCNT = SCNT + 1

                    IF( SCNT .LE. MXFACS ) THEN
                        SPDLIST( SCNT ) = SPDA( I )
                        ISPIDA ( I )    = SCNT

                        GO TO 133

                    ENDIF  ! If overflow or not

                ENDIF      ! If new speed or not

            GO TO 129   ! Loop on list of speeds

  133   CONTINUE        ! Loop on all ofnetspd records

        IF ( SCNT .GT. MXFACS ) THEN

            WRITE( *,92010 )
     &          'Unique speeds record count  ', SCNT,
     &          'MXFACS dimensioned record-count', MXFACS
            CALL M3EXIT ( 'MKMPLIST', 0, 0, 
     &                    'Max record-count exceeded', 2 )
        END IF
   
C.........  Sort OFNETSPD records by MFIP, FIP, Speed index, and RCLAS
        CALL SORTI4( NSPD, SPDIDX, MFIPA, ISPIDA, IFIPA, IRCLSA )

C.........  Create list of all individual combinations.  The rules
C.........  are that a new combination exists if the FIP changes (either
C.........  from a default to a default or a non-default to a non-default), 
C.........  or a new speed index with a default or non-default FIP. Keep in
C.........  mind that MFIP is the FIP code that an actual FIP points to
C.........  to get its emission factor parameters.

        J     = 0
        LMFIP = IMISS3
        LSPID = IMISS3
        PSI   = 0

        DO 144 I = 1, NSPD

            MFIP = MFIPA ( SPDIDX( I ) )
            SPID = ISPIDA( SPDIDX( I ) )
            FIP  = IFIPA ( SPDIDX( I ) )
            RCL  = IRCLSA( SPDIDX( I ) )

C.............  Find MFIP in list of defaults
            D = FIND1( MFIP, NDEF, DEFLST )

C.............  Find road class in compile-time list
            R = FIND1( RCL, NRCLAS, MROADS3 )

C.............  If MFIP is a default
            IF( D .GT. 0 ) THEN

C.................  MFIP is not equal to last time ( new MFIP )
                IF( MFIP .NE. LMFIP ) THEN

                    PSI = PSI + 1                      ! Increase PSI
                    FIP = ( MFIP/1000 ) * 1000         ! Output FIP is default

                    IF( DFSTRCL( R,D ) .EQ. IMISS3 )
     &                                      DFSTRCL( R,D ) = SPID

C.................  Same MFIP as last iteration
                ELSE

C.....................  Check if road class value has not yet been stored
                    IF( DFSTRCL( R,D ) .EQ. IMISS3 ) THEN

                        IF( SPID .NE. LSPID ) THEN
                            PSI = PSI + 1              ! Increase PSI
                            FIP = ( MFIP/1000 ) * 1000 ! Output FIP is default

                        ELSE
                            FIP = ( MFIP/1000 ) * 1000 ! Output FIP is default

                        ENDIF

C.........................  Set Default state road-class speed to PSID value 
                        DFSTRCL( R,D ) = SPID

C.....................  Speed index is not the same as what has been stored
C.....................  This means that this record is not a default after all.
                    ELSEIF( DFSTRCL( R,D ) .GT. 0    .AND. 
     &                      DFSTRCL( R,D ) .NE. SPID       ) THEN

                        IF( SPID .NE. LSPID ) THEN
                            PSI = PSI + 1              ! Increase PSI
                        ENDIF

                    ELSE

                        GO TO 144    ! Not a new MPLIST entry -> to end of loop

                    ENDIF
                ENDIF                ! If new MFIP or not

C.............  If MFIP is not a default
            ELSE

C.................  MFIP is not equal to last time or new speed index
                IF( MFIP .NE. LMFIP .OR.
     &              SPID .NE. LSPID      ) THEN

                    PSI = PSI + 1

                ENDIF

            ENDIF

            J = J + 1

            IF( J .LE. MXPLST ) THEN

                INDEXB( J ) = J
                INDEXC( J ) = J
                IFIP  ( J ) = FIP
                IRCLAS( J ) = RCL
                IPSI  ( J ) = PSI
                ISPID ( J ) = SPID

            ENDIF  ! If overflow or not

            LMFIP = MFIP
            LSPID = SPID

  144   CONTINUE       ! To head of loop on speed entries

        NPLIST = J

        IF( NPLIST .GT. MXPLST ) THEN

            WRITE( *,92010 )
     &          'Record count                   ', NPLIST,
     &          'MXPLST dimensioned record-count', MXPLST
            CALL M3EXIT ( 'MKMPLIST', 0, 0,
     &                    'Max record-count exceeded', 2 )
        END IF

C.........  Prompt for MPLIST output file

        PDEV = PROMPTFFILE(
     &           'Enter logical name for MPLIST output file',
     &           .FALSE., .TRUE., 'MPLIST', 'MKMPLIST' )

C.........  Prompt for MPREF output file

        XDEV = PROMPTFFILE(
     &           'Enter logical name for MPREF output file',
     &           .FALSE., .TRUE., 'MPREF', 'MKMPLIST' )

C.........  Sort on proper parameters for outputting to MPLIST

        CALL SORTI2( NPLIST, INDEXB, IFIP, IRCLAS )

C.........  Write out MPLIST file

        DO 155 I = 1, NPLIST

            J = INDEXB( I )

            IF( IPSI( J ) .LT. 1000 ) THEN

                WRITE( PDEV,93020 )
     &                 IFIP( J ), IRCLAS( J ), 0, ( IPSI( J ), K=1,24 )

            ELSE

                WRITE( PDEV,93021 )
     &                 IFIP( J ), IRCLAS( J ), 0, ( IPSI( J ), K=1,24 )

            ENDIF

  155   CONTINUE

C.........  Sort on PSI numbers

        CALL SORTI1( NPLIST, INDEXC, IPSI )

C.........  Write out MPLIST file

        LPSI = IMISS3
        DO 166 I = 1, NPLIST

            J = INDEXC( I )

            IF( IPSI( J ) .NE. LPSI ) THEN

                SPID = NINT( 10. * SPDLIST( ISPID( J ) ) )

                WRITE( XDEV,93030 ) IPSI( J )
                WRITE( XDEV,93040 ) SPID

                LPSI = IPSI( J )

            ENDIF

  166   CONTINUE

  999   CONTINUE          !  exit program

        CALL M3EXIT( 'MKMPLIST', 0, 0,
     &               'Normal completion of Program MKMPLIST', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, :, I10 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( I2, I3, I5, F4.1 ) ! If change F4.1, speed profile naming bad

93020   FORMAT( I5.5, 2X, I3.3, 1X, I3.3, 24( 1X, I3.3 ) )

93021   FORMAT( I5.5, 2X, I3.3, 1X, I3.3, 24( 1X, I4.4 ) )

93030   FORMAT( I3.3, 1X, '0' )

93040   FORMAT( '0 000 0 001 0 00 ', I4.4,
     &          ' 001 1 0 0 0 0     0     0 0 0     0 0 0 0' )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 5( A, :, I7, :, 1X ) )

94020   FORMAT( A, 1X, I5, 1X, A, I7.7, I3.3, 1X, A, I5 )

94040   FORMAT( A, I2.2 )

        END
