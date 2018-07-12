######  Fortran Complier  ##################################
#
#
######  Generic fortran compilation flag  ##################
#
# FFLAGS= -O2 -framework accelerate -fcheck=bounds
# FFLAGS= -O2 -Wl,-rpath,/Users/jeremy/anaconda3/lib -framework accelerate -fcheck=bounds
# FFLAGS= -O3 -Wl,-rpath,/Users/jeremy/anaconda3/lib -framework accelerate -ffast-math -fopenmp -funroll-loops
# FFLAGS= -O0 -Wl,-rpath,/Users/jeremy/anaconda3/lib -framework accelerate -fopenmp
# FFLAGS= -O3 -Wl,-rpath,/Users/jeremy/anaconda3/lib -framework accelerate
# FFLAGS= -O2 -llapack -ffast-math
# FFLAGS= -O -framework accelerate -fcheck=bounds -fopenmp
# FFLAGS= -O3 -framework accelerate -ffast-math -fopenmp
# FFLAGS= -O3 -framework accelerate -fopenmp
# FFLAGS= -O3 -framework accelerate -fcheck=bounds
FFLAGS= -O0 -framework accelerate -fcheck=bounds -fopenmp -Wall
# FFLAGS= -O3 -framework accelerate 
#
############################################################
#
objects = nrtype.o molprops.o boys.o \
         one_electron_terms.o s_sp_l_terms.o s_p_d_terms.o s_d_l_terms.o \
         hf.o read_infile.o read_basfile.o spread_basis.o \
         time_checker.o  nuclear_repulsion.o   \
         cmpt1e.o cmpt2e.o index_functions.o precalc_kabs.o \
         index_sort.o  \
         incore_hartreefock.o u_incore_hartreefock.o direct_hartreefock.o \
         transform.o ccsd.o mp2.o \
         wrapper_dsygv.o fock_build_incore.o \
         fock_build_initial_diag.o fock_build_non_diag.o fock_build_diag.o evaluate_maxdm.o \
         dealloc_integs.o
hf.x : $(objects)
	gfortran ${FFLAGS} -o hf.x $(objects)

$(objects): %.o : %.f90
	gfortran ${FFLAGS} -c $< -o $@
#
############################################################
#



