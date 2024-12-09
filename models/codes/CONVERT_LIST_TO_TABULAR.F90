PROGRAM CONVERT_LIST_TO_TABULAR
USE SUBS
IMPLICIT NONE

! Reads in a models for B in list form, and outputs a tabular version
! (e.g. CHAOS models to Maglib format).

The

INTEGER ::  LMAX_B, IOS, NP, i, HARMONIC, J, M, L, SINCOS
CHARACTER(300) :: FILENAME_B,  JUNK, FILENAME_B_OUT
CHARACTER(1) :: STRING_HARM_TYPE2

REAL( KIND = 8), ALLOCATABLE, DIMENSION(:) ::  GAUSS


PRINT*, 'ENTER MAX DEGREE FOR B'
READ*, LMAX_B

PRINT*, 'ENTER FILENAME FOR MODEL'
READ*, FILENAME_B

PRINT*, 'ENTER OUTPUT FILENAME FOR MODEL'
READ*, FILENAME_B_OUT


! Read in Gauss coefficients for magnetic field in list format.

 OPEN(11, FILE = FILENAME_B, STATUS = 'OLD', FORM = 'FORMATTED', & 
                             IOSTAT = IOS, ACTION = 'READ')
    IF( IOS .NE. 0) THEN
    PRINT*, 'ERROR IN OPENING FILE ', FILENAME_B
    STOP
    ENDIF
   

     NP = LMAX_B*(LMAX_B +2)
     ALLOCATE( GAUSS(1:NP) )
     GAUSS(:) = 0.0_8
     READ(11, *) JUNK
     DO J = 1, NP
     READ(11,*)  GAUSS(J)
     ENDDO
  
     CLOSE(11)

! only one set of harmonic coefficients are required - the SV coefficients which have the highest truncation. Those for B[obs] and u are simply just the first section of the coefficients.      
      
      ALLOCATE(HARMONICS(1: LMAX_B * (LMAX_B + 2) ) )
      HARMONIC = 1       
      DO L = 1, LMAX_B
      DO M = 0, L
      DO SINCOS = COSINE_HARMONIC, SINE_HARMONIC  !cos is 1, sin is 2.
      IF( M .eq. 0 .AND. SINCOS .eq. 2) CYCLE
      HARMONICS(HARMONIC)%M = M
      HARMONICS(HARMONIC)%SINCOS = SINCOS
      HARMONICS(HARMONIC)%L = L
      HARMONIC= HARMONIC+1
      ENDDO
      ENDDO
      ENDDO

          
       OPEN(22, FILE = FILENAME_B_OUT, STATUS = 'REPLACE', FORM = 'FORMATTED')
       WRITE(22,'(A)') TRIM(JUNK)
       DO J=1, LMAX_B * (LMAX_B + 2)
       IF( HARMONICS(J)%M .NE. 0 .AND. HARMONICS(J)%SINCOS .EQ. COSINE_HARMONIC)  THEN
       WRITE(22,'(I3,x,I3,2(x,ES15.5,x))') HARMONICS(J)%L, HARMONICS(J)%M, GAUSS(J),  GAUSS(J+1) 
       ENDIF

       IF( HARMONICS(J)%M .EQ. 0 .AND. HARMONICS(J)%SINCOS .EQ. COSINE_HARMONIC)  THEN
       WRITE(22,'(I3,x,I3,2(x,ES15.5,x))') HARMONICS(J)%L, HARMONICS(J)%M, GAUSS(J), 0.0
       ENDIF
       ENDDO

       CLOSE(22)


         
STOP
   
 
END PROGRAM CONVERT_LIST_TO_TABULAR