! Setup the compuational grid
! Cell centres: xm, cell edges: xc

subroutine generate_grid
    use parameters, only: Lx, Lz, Nx, Nz, gridtype, str_coeff
    use grid
    use ghost, only: halosize
    use gridtypes
    implicit none
    integer :: i

    ! Grid spacing
    dx = Lx / dble(Nx)


    ! Allocate memory for spatial locations

    allocate( xc( -halosize+1:Nx+halosize  )  )
    allocate( xm( -halosize+1:Nx+halosize  )  )

    allocate( zc( -halosize+1:Nz+halosize  )  )
    allocate( zm( -halosize+1:Nz+halosize  )  )
    allocate( dz( -halosize+1:Nz+halosize  )  )

    ! x always uniform
    do i = -halosize+1,Nx+halosize
        xc(i) = dble(i - 1) * dx
        xm(i) = ( dble(i-1) + 0.5 ) * dx
    enddo

    !------------ Wall-normal grid ------------------------------
    select case  (gridtype)

        case (UNIFORM)
            call generate_uniform_grid(dz,zc,zm,   Nz,halosize,Lz)

        case (TANH)
            call generate_tanh_grid(dz,zc,zm,   Nz,halosize,Lz,str_coeff)

    end select


    call write_grid

end subroutine generate_grid

subroutine generate_uniform_grid(dz,zc,zm,   Nz,halosize,Lz)
    implicit none
    integer :: i
    integer, intent(in) :: Nz, halosize
    real, intent(in) :: Lz
    real, dimension(1-halosize:Nz+halosize), intent(inout) :: dz, zc, zm

    do i = -halosize+1,Nz+halosize
        dz(i) = Lz / dble(Nz)
        zc(i) = dble(i - 1) * dz(i)
        zm(i) = ( dble(i-1) + 0.5 ) * dz(i)
    enddo

end subroutine generate_uniform_grid

subroutine generate_tanh_grid(dz,zc,zm,   Nz,halosize,Lz,str_coeff)
    implicit none
    integer :: i
    integer, intent(in) :: Nz, halosize
    real, intent(in) :: Lz, str_coeff
    real, dimension(1-halosize:Nz+halosize), intent(inout) :: dz, zc, zm
    real, allocatable :: z(:)
    real :: ii
    
    allocate( z(1:(Nz+1)) )

    do i = 0,Nz
        ii = dble(i) / dble(Nz)
        z(i+1) = 0.5 * Lz * (1.0 + tanh(str_coeff*(2.0*ii - 1.0 ) ) / tanh(str_coeff) )
    enddo



    do i = 1,Nz
        dz(i) = z(i+1) - z(i)
        zc(i) = z(i)
        zm(i) = 0.5 * ( z(i) + z(i+1)  )
    enddo

    ! Ghost-cell symmetric padding
    dz(1-halosize) = dz(1) ; dz(Nz + halosize) = dz(Nz)
    zc(1-halosize) = zc(1) - dz(1) ; zc(Nz + halosize) = zc(Nz) + dz(Nz)
    zm(1-halosize) = zm(1) - dz(1) ; zm(Nz + halosize) = zm(Nz) + dz(Nz)

    deallocate(z)

end subroutine generate_tanh_grid


subroutine dealloc_grid
    use grid

    if(allocated(xc)) deallocate(xc)
    if(allocated(xm)) deallocate(xm)

    if(allocated(zc)) deallocate(zc)
    if(allocated(zm)) deallocate(zm)
    if(allocated(dz)) deallocate(dz)

end subroutine dealloc_grid
