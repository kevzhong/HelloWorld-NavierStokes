!----OMP-accelerated 2D Navier--Stokes finite difference solver on a staggered grid -----
! USAGE: Go to params.f90 and modify parameters module as desired

program main
    use parameters
    use rk3
    use fftw3
    use omp_lib
    use grid
    use velfields
    use scalarfields
    use velMemory
    implicit none
    integer :: i,nrk3
    !real :: time
    character(100) :: arg

    ! For timing
    real :: start_time, end_time

    !$omp parallel
    !$omp master
    write(*,*) "Executing solver with ", OMP_GET_NUM_THREADS() , " threads!"
    !$omp end master
    !$omp end parallel

    if (command_argument_count().eq.2) then
        call get_command_argument(1,arg)
        read(arg,'(i10)') nt0
        call get_command_argument(2,arg)
        read(arg,'(i10)') Nt
    endif


    ! Variable initialisation
    call allocFields
    call generate_grid
    call init_fft

    if ( Nt0 .ne. 0) then
        write(*,*) "Reading restart file at timestep ", nt0
        call readRestart
    else
        write(*,*) "Starting from prescribed initial condition..."
        call initialCondition
    endif

    time = 0.0d0

    !--------- Begin time-marching ------------------
    do i = nt0+1,nt
        call cpu_time(start_time)

        if (cflmode) call decide_dt
        
        do nrk3 = 1,3
            call next_rk3(nrk3)
            call updateFields
            call pressurePoisson ! Build div(ustar), solve Poisson, projection update to n+1
            ! call check
        enddo

        call cpu_time(end_time)

        ! call postpro
        call bulk_stats

        time = time + dt

        if ( mod(i,tframe2d) .eq. 0 ) then

            ! 3D
            !if (scalarmode .eqv. .true.) call write3DField(temp(1:Nx,1:Ny,1:Nz),Nx,Ny,Nz,'c',i)

            ! ! 2D slice
            call write2DField(u(1:Nx,Ny/2,1:Nz),Nx,Nz,'uxz',i)
            ! call write2DField(w(1:Nx,Ny/2,1:Nz),Nx,Nz,'w',i)
            if (scalarmode ) call write2DField(temp(1:Nx,Ny/2,1:Nz),Nx,Nz,'cxz',i)
            
        endif

        if ( mod(i,tstat) .eq. 0 ) call dump_profiles(i) ! Write statistics to ASCII files

        write(*,*) "Timestep ", i," Wall-time per step = ", end_time - start_time, "dt = ", dt

    enddo
    !--------- End time-marching -------------------

    ! Write restart files
    call write3DField(u(1:Nx,1:Ny,1:Nz),Nx,Ny,Nz,'u',nt)
    call write3DField(v(1:Nx,1:Ny,1:Nz),Nx,Ny,Nz,'v',nt)
    call write3DField(w(1:Nx,1:Ny,1:Nz),Nx,Ny,Nz,'w',nt)
    if (scalarmode) call write3DField(temp(1:Nx,1:Ny,1:Nz),Nx,Ny,Nz,'c',nt)

    ! Memory deallocation
    call deallocFields
    call dealloc_grid
    call dealloc_fft

    write(*,*) "Finished time-stepping!, exiting program"

  end program main