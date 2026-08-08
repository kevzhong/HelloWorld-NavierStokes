subroutine decide_dt
    use parameters
    use velfields
    use scalarfields
    use grid
    implicit none
    real :: cfl_buffer
    integer :: i,j,k


    cfl_buffer = -huge(0.0)

    !$omp parallel do default(none) collapse(3) &
    !$omp private(i, j, k) &
    !$omp shared(Nx, Ny, Nz, u, v, w, dx, dy, dz, nu, want_cfl,dt) &
    !$omp reduction(max:cfl_buffer)
    do k = 1,Nz
        do j = 1,Ny
            do i = 1,Nx
                ! Advection
                cfl_buffer = max( cfl_buffer, abs(u(i,j,k))* dt / dx    + &
                                              abs(v(i,j,k))* dt / dy    + &
                                              abs(w(i,j,k))* dt / dz(k)    )
                                                    

                ! Viscous
                !cfl_buffer = max(cfl_buffer, 6.0 * nu * dt / min(dx**2,dy**2, dz(k)**2 ) )
            enddo
        enddo
    enddo
    !$omp end parallel do

    if (cfl_buffer .gt. 0.0 ) then
        dt = dt * want_cfl / cfl_buffer
    endif


end subroutine decide_dt

! subroutine decide_dt
!     use parameters
!     use velfields
!     use scalarfields
!     use grid
!     implicit none
!     real :: my_dt, small
!     integer :: i,j,k


!     my_dt = huge(0.0)
!     small = tiny(0.0)

!     !$omp parallel default(none) &
!     !$omp private(i, j, k) &
!     !$omp shared(Nx, Ny, Nz, u, v, w, dx, dy, dz, nu, want_cfl,small) &
!     !$omp reduction(min:my_dt)
!     do k = 1,Nz
!         do j = 1,Ny
!             do i = 1,Nx
!                 ! Advection
!                 my_dt = min( my_dt, want_cfl * (dx    / abs(u(i, j,k ) + small) + &
!                                                 dy    / abs(v(i, j,k ) + small) + &
!                                                 dz(k) / abs(w(i, j, k) + small) ) )

!                 ! Viscous
!                 my_dt = min(my_dt, want_cfl * 0.5 / (3.0*nu) * min(dx**2, dy**2, dz(k)**2 ) )
!             enddo
!         enddo
!     enddo
!     !$omp end parallel

!     dt = my_dt

! end subroutine decide_dt