!----- PARAMETERS: MODIFY TO CHANGE SIMULATION SETUP ----------------
module parameters
    use bctypes
    use implicit_types
    use gridtypes
    implicit none

    ! Grid points
    integer :: Nx = 192
    integer :: Ny = 128
    integer :: Nz = 128

    ! ! Domain size
    ! real :: Lx = 1.0
    ! real :: Ly = 1.0
    ! real :: Lz = 1.0
    
    real :: Lx = 2.0 * 3.141592653589793
    real :: Ly = 1.0 * 3.141592653589793
    real :: Lz = 2.0

    ! Time-stepping: Nt0, Nt can be optionally specified as command-line arguments
    integer :: Nt0 = 0 ! 0 for initial condition, /=0 for restart
    integer :: Nt = 60000 ! No. of timesteps
    real :: dt = 1.0e-4 ! Timestep
    
    logical :: cflmode = .true.
    real :: want_cfl = 1.0


    ! Flow parameters
    real :: nu = 1.0 / 180.0
    real :: mean_dpdx = -1.0

    ! Grid stretching
    ! UNIFORM , TANH, COSINE, ERF, NATURAL
    integer :: gridtype = TANH
    real :: str_coeff = 2.2

    logical :: implicitXYmode = .false.
    integer :: implicit_type = HELMHOLTZ ! ADI or HELMHOLTZ, ignored if implicitXYmode is false

    ! Add-ons
    logical :: scalarmode = .false.
    real :: prandtl = 1.0
    real :: beta_gx = 0.0 ; real :: beta_gy = 0.0 ; real :: beta_gz = 1.0 ! buoyancy forcing


    real :: time


    ! Wall boundary-conditions, DIRICHLET or NEUMANN(not implemented yet)

    integer :: bctype_utop = DIRICHLET
    integer :: bctype_ubot = DIRICHLET
    integer :: bctype_vtop = DIRICHLET
    integer :: bctype_vbot = DIRICHLET
    integer :: bctype_wtop = DIRICHLET
    integer :: bctype_wbot = DIRICHLET

    real :: bcval_utop = 0.0
    real :: bcval_ubot = 0.0
    real :: bcval_vtop = 0.0
    real :: bcval_vbot = 0.0
    real :: bcval_wtop = 0.0
    real :: bcval_wbot = 0.0
    
    integer :: bctype_Ttop = DIRICHLET
    integer :: bctype_Tbot = DIRICHLET

    real :: bcval_Ttop = 0.0
    real :: bcval_Tbot = 1.0

    ! Data writing
    integer :: tframe2d = 1000
    integer :: tstat = 1000

end module parameters

!------------------ END PARAMETER SETUP MODIFICATION ----------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------



! Computational grid
module grid
    implicit none

    real :: dx, dy
    real, allocatable, dimension(:) :: dz
    real, allocatable, dimension(:) :: xc,xm
    real, allocatable, dimension(:) :: yc,ym
    real, allocatable, dimension(:) :: zc,zm

end module grid


! Field variables
module velfields
    implicit none
    real, allocatable, dimension(:,:,:) :: u,v,w,p
    real, allocatable, dimension(:,:,:) :: pseudo_p

end module velfields

module scalarfields
    implicit none
    real, allocatable, dimension(:,:,:) :: temp
end module scalarfields

! Working memory for velocity fields
module velMemory
    implicit none

    ! Work arrays for RHSs
    real, allocatable, dimension(:,:,:) :: rhs_u, rhs_v,rhs_w
    real, allocatable, dimension(:,:,:) :: expl_u, expl_u_m1, expl_v, expl_v_m1, expl_w, expl_w_m1
end module velMemory

module scalarMemory
    implicit none
    ! Work arrays for RHSs
    real, allocatable, dimension(:,:,:) :: rhs_temp
    real, allocatable, dimension(:,:,:) :: expl_c, expl_c_m1

end module scalarMemory

! Working memory for FFT / Poisson solver
module fftMemory
    use fftw3
    implicit none
    ! Work arrays for Helmholtz/FFT solver
    real, allocatable, dimension(:) :: lmb_x_on_dx2
    real, allocatable, dimension(:) :: lmb_y_on_dy2

    type(C_PTR) :: fftw_plan_fwd
    type(C_PTR) :: fftw_plan_bwd
    type(C_PTR) :: ptr1, ptr2, ptr3, ptr4

    complex(C_DOUBLE_COMPLEX), pointer :: rhs_hat(:,:,:), soln_hh_hat(:,:,:)
    real(C_DOUBLE), pointer :: soln_hh(:,:,:), rhs_poisson(:,:,:) ! Poisson equation

end module fftMemory

! Working memory for implicit solve
module implicit
    implicit none
    
    ! TDM coefficients
    real, allocatable, dimension(:) :: ami(:), aci(:), api(:)
    real, allocatable, dimension(:) :: amj(:), acj(:), apj(:)
    real, allocatable, dimension(:) :: amk(:), ack(:), apk(:)

    ! RHS/solution buffers for TDM
    real, allocatable, dimension(:) :: tdm_rhsX1(:), tdm_rhsX2(:) ! Two buffers for two solves in Sherman--Morrison
    real, allocatable, dimension(:) :: tdm_rhsY1(:), tdm_rhsY2(:) ! Two buffers for two solves in Sherman--Morrison
    real, allocatable, dimension(:) :: tdm_rhsZ_r(:), tdm_rhsZ_c ! Real & complex buffer


    ! Provisional solution array (delta u etc to pass in x->y->z)
    real, allocatable, dimension(:,:,:) :: impl_delta

end module implicit

module stats
    implicit none

    ! Mean profiles
    real, allocatable, dimension(:) :: uavg(:), vavg(:), wavg(:)
    real, allocatable, dimension(:) :: cavg(:)

    ! Second-order statistics
    real, allocatable, dimension(:) :: uu_avg(:), vv_avg(:), ww_avg(:)
    real, allocatable, dimension(:) :: uw_avg(:)
    real, allocatable, dimension(:) :: cc_avg(:), cw_avg(:)

    ! Integrated statistics
    real :: uvolavg, tvolavg ! bulk velocities, temperatures
    real :: nududz_bot, nududz_top ! Wall shear stress
    real :: kapdTdz_bot, kapdTdz_top ! Wall heat fluxes


    ! SPECTRA
    



    

end module stats
