F90 = ifort # gfortran -ffree-line-length-0
CONVERT_LIST_TO_TABULAR :  CONVERT_LIST_TO_TABULAR.F90 subs.F90
	$(F90) -c   subs.F90
	$(F90) -o CONVERT_LIST_TO_TABULAR CONVERT_LIST_TO_TABULAR.F90 subs.o  -lfftw3 
    
