
        SUBROUTINE RDDEVEMSPT( FDEV )

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
        USE MODSOURC, ONLY: CSOURC, ISIC, IDIU, IWEK

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL         CHKINT
        CHARACTER*2     CRLF
        INTEGER         FINDCFIRST
        INTEGER         STR2INT
        
        EXTERNAL        CHKINT, CRLF, FINDCFIRST, STR2INT

C...........   SUBROUTINE ARGUMENTS
        INTEGER,          INTENT (IN) :: FDEV   ! unit no. of inv file

C...........   Other local variables
        INTEGER              I, K1                  ! counters and indices
        
        INTEGER              ACTCNT                 ! actual count of processed records
        INTEGER              IOS                    ! I/O status
        INTEGER              IREC                   ! no. records in file
        INTEGER              KEYLEN                 ! length of source key
        INTEGER              SIC                    ! SIC code
        INTEGER              TDIU                   ! hourly profile code
        INTEGER              TWEK                   ! weekly profile code

        CHARACTER(LEN=FIPLEN3) CFIP      ! fip code
        CHARACTER(LEN=PLTLEN3) FCID      ! facility ID
        CHARACTER(LEN=CHRLEN3) SKID      ! stack ID
        CHARACTER(LEN=CHRLEN3) DVID      ! device ID
        
        CHARACTER(LEN=ALLLEN3) SRCKEY                 ! source key
        CHARACTER(LEN=300)     LINE                   ! line from file
        CHARACTER(LEN=300)     MESG                   ! message buffer

        CHARACTER(LEN=16) :: PROGNAME = 'RDDEVEMSPT' ! program name

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
        
C.............  Create source key
            CFIP( 1:1 ) = '0'
            CFIP( 2:3 ) = LINE( 1:2 )
            CFIP( 4:6 ) = LINE( 3:5 )
                        
C.............  Replace blanks with zeros        
            DO I = 1,FIPLEN3
                IF( CFIP( I:I ) == ' ' ) CFIP( I:I ) = '0'
            END DO
            
            FCID = ADJUSTL( LINE(  6:20 ) )
            SKID = ADJUSTL( LINE( 21:32 ) )
            DVID = ADJUSTL( LINE( 33:44 ) )
            
            CALL BLDCSRC( CFIP, FCID, SKID, DVID, CHRBLNK3,
     &                    CHRBLNK3, CHRBLNK3, CHRBLNK3,
     &                    SRCKEY )
            
            KEYLEN = LEN_TRIM( SRCKEY )
            
C.............  Find source key in CSOURC array
            K1 = FINDCFIRST( SRCKEY, NSRC, CSOURC( : )( 1:KEYLEN ) )
        
C.............  If key not found, go to next line
            IF( K1 <= 0 ) CYCLE
        
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
