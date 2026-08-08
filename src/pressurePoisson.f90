subroutine pressurePoisson
    use parameters
    use velfields, only: pseudo_p
    use fftMemory
    use ghost
    use bctypes
    implicit none

    ! See section 3.4.1 of Kajishima & Taira (2016) for treatment of Poisson eqn. on staggered grids
    ! with details on non-uniform grids too

    call build_rhsPoisson
    !call solve_pressurePoisson    
    call solve_helmholtz(pseudo_p(1:Nx,1:Ny,1:Nz), rhs_poisson, 0.0, -1.0, PRESSUREBC, PRESSUREBC )
    call update_ghost_pressure(pseudo_p)
    call projectionUpdate


end subroutine pressurePoisson

subroutine build_rhsPoisson
    use velfields, only: u,v,w
    use grid
    use rk3
    use parameters
    use fftMemory, only: rhs_poisson
    implicit none
    integer :: i,j,k
    real :: divu


    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,j,k,divu) &
    !$omp shared(u,v,w,rhs_poisson,dx,dy,dz,Nx,Ny,Nz,aldt)
    do k = 1,Nz
        do j = 1,Ny
            do i = 1,Nx
                divu = ( u(i+1,j,k) - u(i,j,k) ) / dx + &
                       ( v(i,j+1,k) - v(i,j,k) ) / dy + &
                       ( w(i,j,k+1) - w(i,j,k) ) / dz(k)

                rhs_poisson(i,j,k) = divu / aldt
            enddo
        enddo
    enddo
    !$omp end parallel do
end subroutine build_rhsPoisson

subroutine projectionUpdate
    use velfields
    use rk3
    use grid
    use ghost
    use parameters
    implicit none
    integer :: i, j, k
    real :: half_nualdt
    real :: dzmh, dzph

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,j,k,dzmh) &
    !$omp shared(u,v,w,p,pseudo_p,dx,dy,dz,Nx,Ny,Nz,aldt)
    do k = 1,Nz
        dzmh = 0.5*( dz(k-1) + dz(k  ) )

        do j = 1,Ny
            do i = 1,Nx

                u(i,j,k) = u(i,j,k) - aldt * ( pseudo_p(i,j,k) - pseudo_p(i-1, j  , k  ) ) / dx
                v(i,j,k) = v(i,j,k) - aldt * ( pseudo_p(i,j,k) - pseudo_p(i  , j-1, k  ) ) / dy
                w(i,j,k) = w(i,j,k) - aldt * ( pseudo_p(i,j,k) - pseudo_p(i  , j  , k-1) ) / dzmh

                p(i,j,k) = p(i,j,k) + pseudo_p(i,j,k)

            enddo
        enddo
    enddo
    !$omp end parallel do


    ! Additional Laplacian contribution from implicit treatment
    half_nualdt = 0.5 * nu * aldt

    if (implicitXYmode) then
    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,j,k,dzmh,dzph) &
    !$omp shared(p,pseudo_p,dx,dy,dz,Nx,Ny,Nz,dt,half_nualdt)
    do k = 1,Nz
        dzmh = 0.5*( dz(k-1) + dz(k  ) )
        dzph = 0.5*( dz(k  ) + dz(k+1) )

        do j = 1,Ny
            do i = 1,Nx

            p(i,j,k) = p(i,j,k) - half_nualdt * ( ( pseudo_p(i-1,j,k) - 2.0 * pseudo_p(i,j,k) + pseudo_p(i+1,j,k) ) / dx**2 &
                                                + ( pseudo_p(i,j-1,k) - 2.0 * pseudo_p(i,j,k) + pseudo_p(i,j+1,k) ) / dy**2 &
                                                + (   &
                                                   ( pseudo_p(i,j,k+1) - pseudo_p(i,j,k  ) ) / dzph  - &
                                                   ( pseudo_p(i,j,k  ) - pseudo_p(i,j,k-1) ) / dzmh      ) / dz(k)  &
                                                )

            enddo
        enddo
    enddo
    !$omp end parallel do

else ! Only add z-implicit correction
    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,j,k,dzmh,dzph) &
    !$omp shared(p,pseudo_p,dx,dy,dz,Nx,Ny,Nz,dt,half_nualdt)
    do k = 1,Nz
        dzmh = 0.5*( dz(k-1) + dz(k  ) )
        dzph = 0.5*( dz(k  ) + dz(k+1) )
        do j = 1,Ny
            do i = 1,Nx
            p(i,j,k) = p(i,j,k) - half_nualdt * ( ( &
                                                   ( pseudo_p(i,j,k+1) - pseudo_p(i,j,k  ) ) / dzph  - &
                                                   ( pseudo_p(i,j,k  ) - pseudo_p(i,j,k-1) ) / dzmh      ) / dz(k)  )

            enddo
        enddo
    enddo
    !$omp end parallel do
endif

    call update_ghost_wallsU(u,bctype_ubot,bctype_utop,bcval_ubot,bcval_utop)
    call update_ghost_wallsU(v,bctype_vbot,bctype_vtop,bcval_vbot,bcval_vtop)
    call update_ghost_wallsW(w,bctype_wbot,bctype_wtop,bcval_wbot,bcval_wtop)
    call update_ghost_pressure(p)

end subroutine projectionUpdate