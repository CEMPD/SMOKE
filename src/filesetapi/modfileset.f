
       MODULE MODFILESET

!***********************************************************************
!  Module body starts at line 38
!
!  DESCRIPTION:
!     This module contains the public and private variables and arrays  
!     needed to use the FileSetAPI.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 6/02 by C. Seppanen
!     09/2025 by HT UNC-IE:  Use M3UTILIO; Make use of EMSTRG3.EXT for variable size control 
!
!***************************************************************************
!
! Project Title: FileSetAPI
! File: @(#)$Id$
!
! COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
! All Rights Reserved
!
! Carolina Environmental Program
! University of North Carolina at Chapel Hill
! 137 E. Franklin St., CB# 6116
! Chapel Hill, NC 27599-6116
!
! smoke@unc.edu
!
! Pathname: $Source$
! Last updated: $Date$ 
!
!*************************************************************************
       USE M3UTILIO

       IMPLICIT NONE
       SAVE

!.........  Include files
c      INCLUDE 'PARMS3.EXT'  ! I/O API parameters
c      INCLUDE 'FDESC3.EXT'  ! I/O API file description data structures

!.........  File set information        
       INTEGER              :: NVARSET             ! total number of variables in the file set
       INTEGER              :: NFILESET            ! total number of files in the file set
       INTEGER, ALLOCATABLE :: VARS_PER_FILE( : )  ! number of variables per file

!.........  Arrays for storing variable information (dim: NVARSET)
       INTEGER,       ALLOCATABLE :: VTYPESET( : )  ! variable types
       CHARACTER(16), ALLOCATABLE :: VNAMESET( : )  ! variable names
       CHARACTER(16), ALLOCATABLE :: VUNITSET( : )  ! variable units
       CHARACTER(80), ALLOCATABLE :: VDESCSET( : )  ! variable descriptions
c      CHARACTER(IOVLEN3), ALLOCATABLE :: VNAMESET( : )  ! variable names
c      CHARACTER(IOULEN3), ALLOCATABLE :: VUNITSET( : )  ! variable units
c      CHARACTER(IODLEN3), ALLOCATABLE :: VDESCSET( : )  ! variable descriptions

!.........  Internal wrapper data
       TYPE :: CHAR_PTR_ARRAY
           LOGICAL                :: RDONLY            ! read-only status
           CHARACTER(16), POINTER :: LNAMES( : )       ! logical file names
           CHARACTER(16), POINTER :: VARS( :,: )       ! variable names
       END TYPE
        
       INTEGER                :: NOPENSETS = 0         ! total number of open file sets
       CHARACTER(16)          :: RNAMES( MXFILE3 )     ! logical file names for open file sets
       TYPE( CHAR_PTR_ARRAY ) :: FILE_INFO( MXFILE3 )  ! file information for open file sets
       
       END MODULE MODFILESET
       
