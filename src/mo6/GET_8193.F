      subroutine get_8193(INPREC,INERR,RC) 
C 
C GET_8193 calls the various subroutines involved in reading user input 
C records for 1981-93 emissions. 
C 
C   Changes : (Last change first) 
C 
C    4 Aug 1999 @Dyntel-nh  commented out calls to pr_sr_p and pr_avg_p and 
C               added call to init_sr_avg 
C   16 Jul 1999 @Dyntel-bag 2-000 allocate statement after call to count_user 
C   28 Jun 1999 @Dyntel-NH added comments. 
C   17 Jun 1999 @Dyntel-NH. New subroutine. 
C 
C 
C CALLED BY DSPTCH 
C 
C CALLS User_file, count_user,get_user,pr_sr_p, pr_avg_p 
C 
C  Input on call or from calls: 
C 
C    argument list: INPREC,INERR,RC 
C    common blocks: 
C 
C Argument list 
C    INPREC  - The last input record read from the main input file. 
C    INERR - A running tally of the number of errors. 
C    RC    - The return code (zero if there are no errors). 
C          Return codes less than zero mean end-of-file conditions. 
C          Return codes greater than zero mean various error conditions. 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ----------------------------------------------------- 
C NUM_NH    I    a count of the normal starts records 
C num_a     I    a count of the average running records 
C 
      USE userdat 
      implicit none 
      CHARACTER*(*) inprec 
      INTEGER inerr,rc 
C 
      call user_file(INPREC,INERR,RC) !Check existence of user file 
C 
      call init_sr_avg 
C 
      call get_user() ! Read input file and store user records 
c      call pr_sr_p() ! Print stored running/start records 
c      call pr_avg_p()! Print stored average records. 
C 
      end subroutine 
 
