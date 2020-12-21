#=============================================================================#
#== The following is a tutorial on how you can use the debug.h5 debug       ==#
#== database. this demo uses Julia's HDF5 interface; but any language       ==# 
#== that supports HDF5 can do this as well, albeit with different syntax.   ==# 
#=============================================================================#

#== first, we need to load the HDF5 interface/module ==#
using HDF5

#=============================================================================#
#== The layout of information in each debug.h5 can be divided into two      ==#
#== parts - the pre-SCF information, and information for each SCF           ==#
#== iteration. Pre-SCF iteration information is contained in the            ==#
#== /RHF/Iteration-None H5 group, and information for SCF iteration ITER    ==#
#== is contained in the /RHF/Iteration-ITER H5 group.                       ==#
#=============================================================================#

#=============================================================================#
#== The /RHF/Iteration-None group contains the following datasets:          ==# 
#==   - /S: the overlap matrix                                              ==# 
#==   - /X: the orthogonalization/transformation matrix                     ==# 
#==   - /Guess: the initial guess matrix                                    ==# 
#==   - /E_nuc: the nuclear repulsion energy                                ==# 
#==   - /H: the one-electron hamiltonian matrix                             ==#
#=============================================================================#

#=============================================================================#
#==  Let's access this debug.h5's overlap matrix and print it...            ==# 
#=============================================================================#
println("debug.h5 overlap matrix:")
S_debug = h5read("debug.h5", "/RHF/Iteration-None/S")
display(S_debug); println(); println()

#=============================================================================#
#== Note that all matrices are stored as 1D arrays, and the full matrix     ==#
#== representation is stored.                                               ==#   
#=============================================================================#

#=============================================================================#
#== The /RHF/Iteration-ITER group contains the following groups and         ==#
#== datasets:                                                               ==# 
#==   - /C: the MO coefficients                                             ==#
#==   - /D: the updated density matrix                                      ==#
#==   - /E: a subgroup containing energy information, which includes:       ==#
#==     - /E/EHF1: the D*F subcomponent of the electronic energy            ==# 
#==     - /E/EHF2: the D*H subcomponent of the electronic energy            ==# 
#==     - /E/EHF: the total electronic energy                               ==# 
#==   - /F: a subgroup containing Fock matrix information, which includes:  ==#
#==     - /F/Skeleton: the "skeleton" Fock matrix                           ==# 
#==     - /F/Total: the "total" Fock matrix, including Hcore                ==# 
#==   - /F_evec: the eigenvectors of the total Fock matrix                  ==# 
#=============================================================================#

#=============================================================================#
#== Let's access this debug.h5's first-iteration MO coefficients matrix     ==#
#== and print it...                                                         ==# 
#=============================================================================#
println("debug.h5 first-iteration MO coefficients matrix:")
C_debug = h5read("debug.h5", "/RHF/Iteration-1/C")
display(C_debug); println(); println()

#=============================================================================#
#== Now let's access this debug.h5's first-iteration total Fock matrix      ==#
#== and print it...                                                         ==# 
#=============================================================================#
println("debug.h5 first-iteration total Fock matrix:")
F_total_debug = h5read("debug.h5", "/RHF/Iteration-1/F/Total")
display(F_total_debug); println(); println()

#=============================================================================#
#== And that's it! If you have any questions, let me know.                  ==#
#=============================================================================#
println("READ THE COMMENTS OF THIS DEBUG.JL FILE FOR A FULL TUTORIAL!")
