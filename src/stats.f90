! Calculate various integrated statistics and append to ASCII

subroutine bulk_stats
    use stats
    use parameters
    use velfields
    use scalarfields
    use grid
    implicit none
    integer :: i,j,k


    uvolavg = 0.0
    !$omp parallel do collapse(3) default(none) &
    !$omp private(i, j, k) &
    !$omp shared(Nx, Ny, Nz, u, dx,dy,dz) &
    !$omp reduction(+:uvolavg)
    do k = 1,Nz
        do j = 1,Ny
            do i = 1,Nx
                uvolavg = uvolavg + u(i,j,k) * dz(k)
            enddo
        enddo
    enddo
    !$omp end parallel do

    nududz_bot = 0.0
    nududz_top = 0.0

    !$omp parallel do collapse(2) default(none) &
    !$omp private(i, j) &
    !$omp shared(Nx, Ny, Nz, u, dz,nu) &
    !$omp reduction(+:nududz_bot, nududz_top)
    do j = 1,Ny
        do i = 1,Nx
            nududz_bot = nududz_bot + nu * ( u(i,j,1) - u(i,j,0) ) /   dz(1)  
            nududz_top = nududz_top + nu * ( u(i,j,Nz+1) - u(i,j,Nz) ) /  dz(Nz)  
        enddo
    enddo
    !$omp end parallel do

    nududz_bot = nududz_bot / ( dble(Nx * Ny) )
    nududz_top = nududz_top / ( dble(Nx * Ny) )



    uvolavg = uvolavg / (dble(Nx * Ny ) * Lz)
    write(*,*) "Bulk velocity = ", uvolavg, " Wall stress (bot,top) = ", nududz_bot, nududz_top


    ! Append to ASCII file
    open(unit=10, file='outputdir/stats.txt', status='unknown', action='write', position='append')
    write(10,*) time, uvolavg, nududz_bot, nududz_top
    close(10)


end subroutine bulk_stats

subroutine dump_profiles(it)
    use stats
    use velfields
    use parameters
    use grid
    use scalarfields
    implicit none
    integer :: i,j,k
    integer, intent(in) :: it
    character(len=50) :: filename 

    ! Dump files at iteration it
    ! Dump profiles to ASCII files


    uavg(:) = 0.0
    vavg(:) = 0.0
    wavg(:) = 0.0

    uu_avg(:) = 0.0
    vv_avg(:) = 0.0
    ww_avg(:) = 0.0
    uw_avg(:) = 0.0

    !$omp parallel do collapse(3) default(none) &
    !$omp private(i, j, k) &
    !$omp shared(Nx, Ny, Nz, u,v,w) &
    !$omp reduction(+:uavg, vavg, wavg, uu_avg, vv_avg, ww_avg, uw_avg)
    do k = 1,Nz
        do j = 1,Ny
            do i = 1,Nx
                ! Means
                uavg(k) = uavg(k) + u(i,j,k)
                vavg(k) = vavg(k) + v(i,j,k) 
                wavg(k) = wavg(k) + w(i,j,k)

                ! Second-order statistics
                uu_avg(k) = uu_avg(k) + u(i,j,k)**2
                vv_avg(k) = vv_avg(k) + v(i,j,k)**2
                ww_avg(k) = ww_avg(k) + w(i,j,k)**2
                uw_avg(k) = uw_avg(k) + u(i,j,k)*w(i,j,k) ! neglecting stagger for simplicity
            enddo
        enddo
    enddo
    !$omp end parallel do

    ! Normalize
    uavg(:) = uavg(:) / (dble(Nx*Ny))
    vavg(:) = vavg(:) / (dble(Nx*Ny))
    wavg(:) = wavg(:) / (dble(Nx*Ny))

    uu_avg(:) = uu_avg(:) / (dble(Nx*Ny))
    vv_avg(:) = vv_avg(:) / (dble(Nx*Ny))
    ww_avg(:) = ww_avg(:) / (dble(Nx*Ny))
    uw_avg(:) = uw_avg(:) / (dble(Nx*Ny))


    ! Write mean profiles to files
    write(filename, '("outputdir/uavg_it", I0, ".txt")') it
    open(unit=11, file=filename, status='unknown', action='write', position='append')
    do k = 1,Nz
        write(11,*) zm(k), uavg(k), vavg(k), wavg(k)
    enddo
    close(11)

    ! Second-order statistics
    write(filename, '("outputdir/uu_avg_it", I0, ".txt")') it
    open(unit=12, file=filename, status='unknown', action='write', position='append')
    do k = 1,Nz
        write(12,*) zm(k), uu_avg(k), vv_avg(k), ww_avg(k), uw_avg(k)
    enddo
    close(12)



    ! Scalar
    if (scalarmode) then

    cavg(:) = 0.0
    cc_avg(:) = 0.0
    cw_avg(:) = 0.0

    !$omp parallel do collapse(3) default(none) &
    !$omp private(i, j, k) &
    !$omp shared(Nx, Ny, Nz, temp, w) &
    !$omp reduction(+:cavg, cc_avg, cw_avg)
    do k = 1,Nz
        do j = 1,Ny
            do i = 1,Nx
                cavg(k) = cavg(k) + temp(i,j,k)

                cc_avg(k) = cc_avg(k) + temp(i,j,k)**2
                cw_avg(k) = cw_avg(k) + temp(i,j,k)*w(i,j,k) ! neglecting stagger for simplicity

            enddo
        enddo
    enddo
    !$omp end parallel do

    cavg(:) = cavg(:) / (dble(Nx*Ny))
    cc_avg(:) = cc_avg(:) / (dble(Nx*Ny))
    cw_avg(:) = cw_avg(:) / (dble(Nx*Ny))


    ! Write mean profiles to files
    write(filename, '("outputdir/cavg_it", I0, ".txt")') it
    open(unit=11, file=filename, status='unknown', action='write', position='append')
    do k = 1,Nz
        write(11,*) zm(k), cavg(k)
    enddo
    close(11)

    ! Second-order statistics
    write(filename, '("outputdir/cc_avg_it", I0, ".txt")') it
    open(unit=12, file=filename, status='unknown', action='write', position='append')
    do k = 1,Nz
        write(12,*) zm(k), cc_avg(k), cw_avg(k)
    enddo
    close(12)

    endif


end subroutine dump_profiles