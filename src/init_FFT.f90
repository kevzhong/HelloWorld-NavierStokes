subroutine init_fft
    use parameters, only: Nx, Ny, Nz
    use grid, only : dx, dy
    use fftw3
    use fftMemory
    use velMemory
    implicit none

    integer :: i
    real :: PI = 2.d0*dasin(1.d0) 


    !---------------------------- Modified wavenumbers --------------------------------
    allocate( lmb_x_on_dx2(Nx/2+1) )
    allocate( lmb_y_on_dy2(1:Ny) )

    do i = 1,Nx/2+1
        lmb_x_on_dx2(i) = ( 2.0 * cos(2.0*PI*(real(i)-1.0) / Nx ) - 2.0 ) / dx**2
    enddo

    do i = 1,Ny
        lmb_y_on_dy2(i) = ( 2.0 * cos(2.0*PI*(real(i)-1.0) / Ny ) - 2.0 ) / dy**2
    enddo

    !--------------- FFT plans for Helmholtz solver using C-type allocations --------

    ! Solution array real(Nx,Ny,Nz)
    ptr1 = fftw_alloc_real(int(Nx * Ny * Nz, C_SIZE_T))
    call c_f_pointer(ptr1, soln_hh, [Nx,Ny,Nz])

    ptr2 = fftw_alloc_complex(int((Nx/2+1) * Ny * Nz, C_SIZE_T))
    call c_f_pointer(ptr2, soln_hh_hat, [Nx/2+1,Ny,Nz])

    ptr3 = fftw_alloc_complex(int((Nx/2+1) * Ny * Nz, C_SIZE_T))
    call c_f_pointer(ptr3, rhs_hat, [Nx/2+1,Ny,Nz])

    ptr4 = fftw_alloc_real(int(Nx * Ny * Nz, C_SIZE_T))
    call c_f_pointer(ptr4, rhs_poisson, [Nx,Ny,Nz])

    ! ! 2D real to complex plan:  double(Nx,1) ---> complex(Nx/2+1, 1)
    ! fftw_plan_fwd = fftw_plan_dft_r2c_2d(Nx,Ny, rhs_poisson(:,:,1), rhs_hat(:,:,1), FFTW_ESTIMATE)

    ! ! 2D complex to real transform: complex(Nx/2+1, 1) ----> double(Nx,1)
    ! fftw_plan_bwd = fftw_plan_dft_c2r_2d(Nx,Ny, rhs_hat(:,:,1), soln_hh(:,:,1), FFTW_ESTIMATE)

    fftw_plan_fwd = fftw_plan_dft_r2c_2d( Ny, Nx, rhs_poisson(:,:,1), rhs_hat(:,:,1), FFTW_ESTIMATE)

    fftw_plan_bwd = fftw_plan_dft_c2r_2d( Ny, Nx, rhs_hat(:,:,1), soln_hh(:,:,1), FFTW_ESTIMATE)


end subroutine init_fft


subroutine dealloc_fft
    use fftMemory
    implicit none

    if(allocated(lmb_x_on_dx2)) deallocate(lmb_x_on_dx2)
    if(allocated(lmb_y_on_dy2)) deallocate(lmb_y_on_dy2)

    ! Destroy the FFTW plan
    call fftw_destroy_plan(fftw_plan_fwd)
    call fftw_destroy_plan(fftw_plan_bwd)
    call fftw_free(ptr1)
    call fftw_free(ptr2)
    call fftw_free(ptr3)
    call fftw_free(ptr4)

end subroutine dealloc_fft
