
        SUBROUTINE RDPROCEMSPT( FDEV, NCNTY, FIPTOCSRC )

C***********************************************************************
C  subroutine body starts at line 232
C
C  DESCRIPTION:
C      This subroutine reads the EMS-95 point process file.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Copied from rdemspt.f by C. Seppanen (2/03)
C
C****************************************************************************
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

C...........   MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC
        
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: CSOURC, CSCC

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL         CHKINT
        CHARACTER*2     CRLF
        INTEGER         FINDCFIRST
        INTEGER         FIND1
        INTEGER         STR2INT
        
        EXTERNAL        CHKINT, CRLF, FINDCFIRST, FIND1, STR2INT

C...........   SUBROUTINE ARGUMENTS
        INTEGER,          INTENT (IN) :: FDEV   ! unit no. of inv file
        INTEGER,          INTENT (IN) :: NCNTY  ! no. counties in inventory
        INTEGER,          INTENT (IN) :: FIPTOCSRC( NCNTY+1,2 )  ! index into CSOURC

C...........   Other local variables
        INTEGER              I, K1                  ! counters and indices
        
        INTEGER              ACTCNT                 ! actual count of processed records
        INTEGER              IOS                    ! I/O status
        INTEGER              IREC                   ! no. records in file
        INTEGER              IFIP                   ! integer FIPS code
        INTEGER              STIDX, ENDIDX          ! start and end idx into CSOURC
        INTEGER              KEYLEN                 ! length of source key
        
        CHARACTER(LEN=FIPLEN3) CFIP      ! fip code
        CHARACTER(LEN=PLTLEN3) FCID      ! facility ID
        CHARACTER(LEN=CHRLEN3) SKID      ! stack ID
        CHARACTER(LEN=CHRLEN3) DVID      ! device ID
        CHARACTER(LEN=CHRLEN3) PRID      ! process ID
        
        CHARACTER(LEN=SCCLEN3) TSCC                 ! SCC code
        CHARACTER(LEN=ALLLEN3) SRCKEY               ! source key
        CHARACTER(LEN=300)     LINE                 ! line from file
        CHARACTER(LEN=300)     MESG                 ! message buffer

        CHARACTER(LEN=16) :: PROGNAME = 'RDPROCEMSPT' ! program name

C***********************************************************************
C   begin body of subroutine RDPROCEMSPT

        ACTCNT = 0
        IREC = 0

C.........  Loop through lines in file
        DO
    
            READ( FDEV, 93000, IOSTAT=IOS ) LINE
            IREC = IREC + 1
            
C.............  Check for I/O errors
            IF( IOS > 0 ) THEN
                WRITE( MESG, 94010 ) 'I/O error', IOS,
     &              'reading process file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF
    
C.............  Check for end of file
            IF( IOS < 0 ) EXIT
    
C.............  Read FIPS code from file
            CFIP( 1:1 ) = '0'
            CFIP( 2:3 ) = ADJUSTR( LINE( 1:2 ) )
            CFIP( 4:6 ) = ADJUSTR( LINE( 3:5 ) )
                        
C.............  Replace blanks with zeros        
            DO I = 1,FIPLEN3
                IF( CFIP( I:I ) == ' ' ) CFIP( I:I ) = '0'
            END DO
            
C.............  Find county in FIPTOCSRC array
            IFIP = STR2INT( CFIP( 1:6 ) )
            K1 = FIND1( IFIP, NCNTY+1, FIPTOCSRC( :,1 ) )

C.............  Make sure county is in the inventory
            IF( K1 <= 0 ) CYCLE
            
C.............  Set starting and ending indices            
            STIDX = FIPTOCSRC( K1,2 )
            ENDIDX = FIPTOCSRC( K1+1,2 ) - 1
            
C.............  Build source key
            FCID = ADJUSTL( LINE(  6:20 ) )
            SKID = ADJUSTL( LINE( 21:32 ) )
            DVID = ADJUSTL( LINE( 33:44 ) )
            PRID = ADJUSTL( LINE( 45:56 ) )
            
            CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID,
     &                    CHRBLNK3, CHRBLNK3, CHRBLNK3,
     &                    SRCKEY )
     
            KEYLEN = LEN_TRIM( SRCKEY )
            
C.............  Find source key in CSOURC array
            K1 = FINDCFIRST( SRCKEY, ENDIDX-STIDX+1, 
     &                       CSOURC( STIDX:ENDIDX )( 1:KEYLEN ) )
    
C.............  If key not found, go to next line
            IF( K1 <= 0 ) CYCLE

C.............  Shift index to account for only searching part of CSOURC array
            K1 = K1 + STIDX - 1
                
C.............  Read and check SCC
            TSCC = LINE( 57:64 )
            
            IF( .NOT. CHKINT( TSCC ) ) THEN
                WRITE( MESG,94010 ) 
     &                 'Badly formatted SCC in process ' //
     &                 'file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
            IF( TSCC == ' ' ) THEN
                WRITE( MESG,94010 ) 
     &                 'Missing SCC in process file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF
    
            IF( STR2INT( TSCC ) <= 9999999 ) THEN  ! SCC must be 8 digits
                WRITE( MESG,94010 ) 
     &                 'Invalid SCC "' // TSCC // 
     &                 '" in process file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF
    
C.............  Increment count of actual records            
            ACTCNT = ACTCNT + 1
            
C.............  Loop through matching sources and set values
            DO
            
C.................  Check for end of array or end of matching sources
                IF( K1 > NSRC ) EXIT
                IF( CSOURC( K1 )( 1:KEYLEN ) /= SRCKEY ) EXIT
                CALL PADZERO( TSCC )
                CSCC( K1 ) = TSCC
                
                K1 = K1 + 1
            END DO
    
        END DO

C.........  Write count of total records processed            
        WRITE( MESG,94010 )
     &         'PROCESS FILE processed' // CRLF() // BLANK10 //
     &         'Actual PROCESS record-count ', ACTCNT
        CALL M3MSG2( MESG )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94125   FORMAT( I5 )

94300   FORMAT( A, I2.2, A, I2.2, A )

        END SUBROUTINE RDPROCEMSPT
