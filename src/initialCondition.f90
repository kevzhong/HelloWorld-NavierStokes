subroutine initialCondition
    use velfields
    use grid
    use parameters
    use ghost
    use scalarfields
    implicit none
    integer :: i, j,k
    real(8) :: PI = 2.d0*dasin(1.d0) 
    real :: eps
    real :: Ub, Av, Aw

    integer, parameter :: nroll  = 3   ! Number of roller pairs across Ly
    integer, parameter :: nxmode = 1   ! Streamwise modulation mode

    real, parameter :: roll_mod = 0.25

    real :: h
    real :: alpha, beta
    real :: phase_roll, phase_streak

    real :: eta_u, eta_v, eta_w
    real :: G_u, G_v, G_w
    real :: F_w
    real :: H_v, dHdy_w
    real :: xfac
    real :: u_lam, u_streak

    Ub = 20.0
    Av = 0.03
    Aw = 0.08


h = 0.5*Lz

alpha = 2.0*PI*real(nxmode)/Lx
beta  = 2.0*PI*real(nroll )/Ly

phase_roll   = PI/3.0
phase_streak = PI/5.0

!$omp parallel do collapse(3)                                      &
!$omp default(none)                                                &
!$omp private(i,j,k,eta_u,eta_v,eta_w,G_u,G_v,G_w,F_w,            &
!$omp         H_v,dHdy_w,xfac,u_lam,u_streak)                     &
!$omp shared(u,v,w,xm,yc,ym,zc,zm,Nx,Ny,Nz,Lx,Ly,Lz,PI,          &
!$omp        Ub,Av,Aw,h,alpha,beta,phase_roll,phase_streak)
do k = 1,Nz
    do j = 1,Ny
        do i = 1,Nx

            ! ----------------------------------------------------------
            ! Normalised wall-normal coordinates
            ! ----------------------------------------------------------

            ! u and v are located at wall-normal cell centres zm(k)
            eta_u = (zm(k) - h)/h
            eta_v = (zm(k) - h)/h

            ! w is located at wall-normal faces zc(k)
            eta_w = (zc(k) - h)/h

            G_u = 1.0 - eta_u*eta_u
            G_v = 1.0 - eta_v*eta_v
            G_w = 1.0 - eta_w*eta_w

            ! Guard against tiny negative values caused by roundoff at
            ! the physical walls.
            G_u = max(G_u,0.0)
            G_v = max(G_v,0.0)
            G_w = max(G_w,0.0)

            F_w = G_w*G_w

            ! ----------------------------------------------------------
            ! Weak streamwise modulation of the rollers
            !
            ! This breaks the exact streamwise-independent symmetry.
            ! Multiplying both v and w by the same x-dependent factor
            ! does not destroy their y-z divergence-free property.
            ! ----------------------------------------------------------

            xfac = 1.0 + roll_mod*cos(alpha*xm(i))

            ! ----------------------------------------------------------
            ! Laminar plane-Poiseuille profile
            !
            ! This profile has:
            !   u(0)    = 0
            !   u(Lz)   = 0
            !   u(Lz/2) = 1.5*Ub
            !   volume average = Ub
            ! ----------------------------------------------------------

            u_lam = 6.0*Ub*(zm(k)/Lz)*(1.0 - zm(k)/Lz)

            ! ----------------------------------------------------------
            ! Streamwise streak perturbation
            !
            ! This has zero spanwise mean and therefore does not alter
            ! the bulk velocity.
            ! ----------------------------------------------------------

            u_streak = Aw*Ub*G_u *                                &
                       ( sin(beta*ym(j))                            &
                       + 0.5*sin(2.0*beta*ym(j) + phase_streak) )

            u(i,j,k) = u_lam + u_streak

            ! ----------------------------------------------------------
            ! Roller streamfunction:
            !
            ! Psi = Av*Ub*h*(1-eta^2)^2*H(y)*xfac
            !
            ! v = d(Psi)/dz
            ! w = -d(Psi)/dy
            ! ----------------------------------------------------------

            H_v = cos(beta*yc(j))                                  &
                + 0.5*cos(2.0*beta*yc(j) + phase_roll)

            ! d[(1-eta^2)^2]/deta = -4*eta*(1-eta^2)
            !
            ! Since deta/dz = 1/h, the factor h in Psi cancels.
            v(i,j,k) = -4.0*Av*Ub*eta_v*G_v*H_v*xfac

            ! Derivative of H evaluated at the w location:
            !
            ! dH/dy = -beta*sin(beta*y)
            !         -beta*sin(2*beta*y + phase_roll)
            dHdy_w = -beta*sin(beta*ym(j))                          &
                     -beta*sin(2.0*beta*ym(j) + phase_roll)

            w(i,j,k) = -Av*Ub*h*F_w*dHdy_w*xfac

        enddo
    enddo
enddo
!$omp end parallel do

    call update_ghost_wallsU(u,bctype_ubot,bctype_utop,bcval_ubot,bcval_utop)
    call update_ghost_wallsU(v,bctype_vbot,bctype_vtop,bcval_vbot,bcval_vtop)
    call update_ghost_wallsW(w,bctype_wbot,bctype_wtop,bcval_wbot,bcval_wtop)
    !call update_ghost_pressure(p)

    if (scalarmode .eqv. .true.) then
        call random_seed()
        !$omp parallel do &
        !$omp default(none) &
        !$omp private(eps,i,j,k) &
        !$omp shared(temp,zm,bcval_Tbot,Lz,bcval_Ttop,Nx,Ny,Nz)
        do k = 1,Nz
            do j = 1,Ny
                do i = 1,Nx
                    temp(i,j,k) = (bcval_Ttop - bcval_Tbot) / Lz * zm(k) + bcval_Tbot 

                    ! Super-impose perturbation
                    call random_number(eps)
                    temp(i,j,k) = temp(i,j,k) + 0.2 * (eps - 0.5)
                enddo
            enddo
        enddo
        !$omp end parallel do

    call update_ghost_wallTemp(temp,bctype_Tbot,bctype_Ttop,bcval_Tbot,bcval_Ttop)
    
    endif

end subroutine initialCondition


subroutine readRestart
    use velfields
    use grid
    use parameters
    use ghost
    use scalarfields
    implicit none

    call read3DField(u(1:Nx,1:Ny,1:Nz),Nx,Ny,Nz,'u',nt0)
    call read3DField(v(1:Nx,1:Ny,1:Nz),Nx,Ny,Nz,'v',nt0)
    call read3DField(w(1:Nx,1:Ny,1:Nz),Nx,Ny,Nz,'w',nt0)
    call update_ghost_wallsU(u,bctype_ubot,bctype_utop,bcval_ubot,bcval_utop)
    call update_ghost_wallsU(v,bctype_vbot,bctype_vtop,bcval_vbot,bcval_vtop)
    call update_ghost_wallsW(w,bctype_wbot,bctype_wtop,bcval_wbot,bcval_wtop)

    if (scalarmode .eqv. .true. ) then
        call read3DField(temp(1:Nx,1:Ny,1:Nz),Nx,Ny,Nz,'c',nt0)
        call update_ghost_wallTemp(temp,bctype_Tbot,bctype_Ttop,bcval_Tbot,bcval_Ttop)
    endif

end subroutine readRestart