
        SUBROUTINE RDDEVEMSPT( FDEV, NCNTY, FIPTOCSRC )

C***********************************************************************
C  subroutine body starts at line 232
C
C  DESCRIPTION:
C      This subroutine reads the EMS-95 point device file.
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

C...........   MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC
        
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: CSOURC, ISIC, IDIU, IWEK

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL         CHKINT
        CHARACTER(2)    CRLF
        INTEGER         FINDCFIRST
        INTEGER         FIND1
        INTEGER         STR2INT
        
        EXTERNAL        CHKINT, CRLF, FINDCFIRST, FIND1, STR2INT

C...........   SUBROUTINE ARGUMENTS
        INTEGER,          INTENT (IN) :: FDEV    ! unit no. of inv file
        INTEGER,          INTENT (IN) :: NCNTY   ! no. counties in inventory
        INTEGER,          INTENT (IN) :: FIPTOCSRC( NCNTY+1,2 )  ! index into CSOURC

C...........   Other local variables
        INTEGER              I, K1                  ! counters and indices
        
        INTEGER              ACTCNT                 ! actual count of processed records
        INTEGER              IOS                    ! I/O status
        INTEGER              IREC                   ! no. records in file
        INTEGER              IFIP                   ! integer FIPS code
        INTEGER              STIDX, ENDIDX          ! start and end idx into CSOURC
        INTEGER              KEYLEN                 ! length of source key
        INTEGER              SIC                    ! SIC code
        INTEGER              TDIU                   ! hourly profile code
        INTEGER              TWEK                   ! weekly profile code

        CHARACTER(FIPLEN3) CFIP      ! fip code
        CHARACTER(PLTLEN3) FCID      ! facility ID
        CHARACTER(CHRLEN3) SKID      ! stack ID
        CHARACTER(CHRLEN3) DVID      ! device ID
        
        CHARACTER(ALLLEN3) SRCKEY                 ! source key
        CHARACTER(300)     LINE                   ! line from file
        CHARACTER(300)     MESG                   ! message buffer

        CHARACTER(16) :: PROGNAME = 'RDDEVEMSPT' ! program name

C***********************************************************************
C   begin body of subroutine RDDEVEMSPT

        ACTCNT = 0
        IREC = 0

C.........  Loop through lines in file
        DO
        
            READ( FDEV, 93000, IOSTAT=IOS ) LINE
            IREC = IREC + 1
            
C.............  Check for I/O errors
            IF( IOS > 0 ) THEN
                WRITE( MESG, 94010 ) 'I/O error', IOS,
     &                 'reading device file at line', IREC
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
            
            CALL BLDCSRC( CFIP, FCID, SKID, DVID, CHRBLNK3,
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
        
C.............  Read and check SIC
            IF( .NOT. CHKINT( LINE( 45:48 ) ) ) THEN
                WRITE( MESG,94010 )
     &                 'Badly formatted SIC in device file at line',
     &                 IREC, CRLF() // BLANK10 // 'Setting to 0000'
                CALL M3MESG( MESG )
                SIC = 0
            ELSE
                SIC = STR2INT( LINE( 45:48 ) )
                
                IF( SIC == 0 ) THEN
                    WRITE( MESG,94010 )
     &                     'Default SIC "0000" in device file at line', 
     &                     IREC
                    CALL M3MESG( MESG )
                ELSE IF( SIC < 111 ) THEN   ! valid SIC codes from 0111 to 9999
                    WRITE( MESG,94010 )
     &                     'Invalid SIC "' // LINE( 45:48 ) //
     &                     '" in device file at line', IREC,
     &                     CRLF() // BLANK10 // 'Setting to 0000'
                    CALL M3MESG( MESG )
                    SIC = 0
                END IF
            END IF
        
C.............  Read and check temporal profile numbers
            TDIU = STR2INT( LINE( 121:122 ) )
            TWEK = STR2INT( LINE( 123:124 ) )
            
            IF( TDIU < 0 ) TDIU = 0   ! treat missing as default
            IF( TWEK < 0 ) TWEK = 0
            
            IF( TDIU == 0 ) THEN
                WRITE( MESG,94010 )
     &                 'Default hourly profile', TDIU,
     &                 'in device file at line', IREC
                CALL M3MESG( MESG )
            END IF
            
            IF( TWEK == 0 ) THEN
                WRITE( MESG,94010 )
     &                 'Default weekly profile', TWEK,
     &                 'in device file at line', IREC
                CALL M3MESG( MESG )
            END IF
            
C.............  Increment count of actual records            
            ACTCNT = ACTCNT + 1
            
C.............  Loop through matching sources and set values
            DO

C.................  Check for end of array or end of matching sources
                IF( K1 > NSRC ) EXIT
                IF( CSOURC( K1 )( 1:KEYLEN ) /= SRCKEY ) EXIT
                
                ISIC( K1 ) = SIC
                IDIU( K1 ) = TDIU
                IWEK( K1 ) = TWEK
                
                K1 = K1 + 1
            END DO
        
        END DO
    
C.........  Write count of total records processed    
        WRITE( MESG,94010 )
     &         'DEVICE FILE processed' // CRLF() // BLANK10 //
     &         'Actual DEVICE record-count ', ACTCNT
        CALL M3MSG2( MESG )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94125   FORMAT( I5 )

94300   FORMAT( A, I2.2, A, I2.2, A )

        END SUBROUTINE RDDEVEMSPT
