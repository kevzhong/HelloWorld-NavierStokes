! Accumulates the RHS terms in the linear system for velocity components
! That is, the advective, diffusive, and pressure-gradient terms
! Since the current treatment is fully explicit for all terms using OpenMP
! everything has been accumulated in a single (Nx,Nz) loop so that only a single instance
! of OMP thread creation/destruction is needed for brevity

subroutine rhsVelocity
    implicit none
    call build_rhsVelocity
end subroutine rhsVelocity

! Build all RHS vectors for velocity
subroutine build_rhsVelocity
    use velfields
    use rk3
    use parameters
    use grid
    use velMemory
    use scalarfields
    implicit none
    integer :: i, j,k
    real :: duudx, duvdy, duwdz, nudd2_u, nudd2xy_u, dpdx !u-velocity
    real :: dvudx, dvvdy, dvwdz, nudd2_v, nudd2xy_v, dpdy !v-velocity
    real :: dwudx, dwvdy, dwwdz, nudd2_w, nudd2xy_w, dpdz !w-velocity
    real :: Tu_interp, Tv_interp, Tw_interp
    real :: dzph, dzmh

    ! Staggered grid arrangement, the stagger is by -1/2, index convention:

    !                     v(i,j+1,k)
    !                      ^
    !                      |
    !                      |
    !    __________________|___________________
    !   |                                      |
    !   |                                      |
    !   |                                      |
    !   |                                      |
    !   |  u(i,j,k)                            |
    ! ----->                               ------>  u(i+1,j,k)
    !   |                  O                   |
    !   |              temp(i,j,k)             |
    !   |              or p(i,j,k)             |
    !   |                                      |
    !   |                                      |
    !   |                  ^                   |
    !   |__________________|___________________
    !                      |
    !                      |
    !                     v(i,j,k)

    !$omp parallel do &
    !$omp default(none) &
    !$omp private(i,j,k) &
    !$omp private(dzmh,dzph) &
    !$omp private(duudx, duvdy, duwdz, nudd2_u, nudd2xy_u, dpdx) &
    !$omp private(dvudx, dvvdy, dvwdz, nudd2_v, nudd2xy_v, dpdy) &
    !$omp private(dwudx, dwvdy, dwwdz, nudd2_w, nudd2xy_w, dpdz) &
    !$omp private(Tu_interp, Tv_interp, Tw_interp) &
    !$omp shared(Nx,Ny,Nz,dx,dy,dz,nu,mean_dpdx,scalarmode,implicitXYmode,beta_gx,beta_gy,beta_gz,gamdt,zetdt,aldt) &
    !$omp shared(u,v,w,p,temp,rhs_u,rhs_v,rhs_w,expl_u,expl_v,expl_w,expl_u_m1,expl_v_m1,expl_w_m1)
    do k = 1,Nz
        dzmh = 0.5*( dz(k-1) + dz(k  ) )
        dzph = 0.5*( dz(k  ) + dz(k+1) )
        do j = 1,Ny
            do i = 1,Nx

            !--------------------- U VELOCITY -----------------------!

            ! d uu  |                 1     [                  ]
            ! ----- |             =  ----   |   uu    -   uu   | 
            !  dx   |i+1/2,j,k        dx    [    i+1        i  ]

            duudx =( (u(i+1,j,k)+u(i,j,k)) &
                     *(u(i+1,j,k)+u(i,j,k)) &
                     -(u(i-1,j,k)+u(i,j,k)) &
                     *(u(i-1,j,k)+u(i,j,k)) &
                     ) / (4.0 * dx)

            ! d uv  |                1     [                  ]
            ! ----- |            =  ----   |   uv    -   uv   | 
            !  dy   |i+1/2,j,k       dy    [    i+1        i  ]

            duvdy = ( (v(i,j+1,k)+v(i-1,j+1,k))*(u(i,j+1,k)+u(i,j,k)) &
                        -(v(i,j,k)+v(i-1,j,k))*(u(i,j,k)+u(i,j-1,k)) &
                         ) / (4.0 * dy )  

            ! d uw  |               1     [                  ]
            ! ----- |           =  ----   |   uw    -   uw   | 
            !  dz   |i+1/2, k       dz    [    i+1        i  ] 

            duwdz = ( (w(i,j,k+1)+w(i-1,j,k+1))*(u(i,j,k+1)+u(i,j,k)) &
                    -(w(i,j,k)+w(i-1,j,k))*(u(i,j,k)+u(i,j,k-1)) &
                     ) / (4.0 * dz(k) )      

            ! d2u / dxj2
            nudd2xy_u = nu *  ( ( u(i-1,j,k) -2.0*u(i,j,k) + u(i+1,j,k)  ) / dx**2  + &
                                ( u(i,j-1,k) -2.0*u(i,j,k) + u(i,j+1,k)  ) / dy**2  )

            nudd2_u = nu  * (   ( u(i,j,k+1) - u(i,j,k  ) ) / dzph  - &
                                ( u(i,j,k  ) - u(i,j,k-1) ) / dzmh   ) / dz(k)

            ! dp/dx | i+1/2
            dpdx = ( -p(i-1,j,k) + p(i,j,k) ) / dx

            !--------------------- V VELOCITY -----------------------!

            ! d vv  |                 1     [                  ]
            ! ----- |             =  ----   |    vv   -   vv   | 
            !  dy   |i,j+1/2,k        dy    [    j+1        j  ]

            
            dvvdy = ( (v(i,j+1,k)+v(i,j,k)) &
                    *(v(i,j+1,k)+v(i,j,k)) &
                    -(v(i,j-1,k)+v(i,j,k)) &
                    *(v(i,j-1,k)+v(i,j,k)) &
                    ) / (4.0 * dy)


            ! d vu  |                 1     [                  ]
            ! ----- |             =  ----   |    vu   -   vu   | 
            !  dx   |i,j+1/2,k        dx    [    j+1        j  ]

            dvudx = ( (u(i+1,j,k)+u(i+1,j-1,k))*(v(i+1,j,k)+v(i,j,k)) &
                    -(u(i,j,k)+u(i,j-1,k))*(v(i,j,k)+v(i-1,j,k)) &
                     ) / (4.0 * dx )  

            ! d vw  |                1     [                  ]
            ! ----- |            =  ----   |   vw    -   vw   | 
            !  dz   |i+1/2,j,k       dz    [    j+1        j  ] 

            dvwdz = ( (w(i,j,k+1) + w(i,j-1,k+1)) * &
                       (v(i,j,k+1) + v(i,j,k))       &
                      -(w(i,j,k)   + w(i,j-1,k))     * &
                       (v(i,j,k)   + v(i,j,k-1))     &
                     ) / (4.0 * dz(k))

            ! d2v / dxj2
            nudd2xy_v = nu *  ( ( v(i-1,j,k) -2.0*v(i,j,k) + v(i+1,j,k)  ) / dx**2  + &
                                ( v(i,j-1,k) -2.0*v(i,j,k) + v(i,j+1,k)  ) / dy**2  )

             nudd2_v = nu  * (   ( v(i,j,k+1) - v(i,j,k  ) ) / dzph  - &
                                 ( v(i,j,k  ) - v(i,j,k-1) ) / dzmh   ) / dz(k)

            ! dp/dy | j+1/2
             dpdy = ( -p(i,j-1,k) + p(i,j,k) ) / dy

            !--------------------- W VELOCITY -----------------------!
                     
            ! d ww  |                1     [                  ]
            ! ----- |            =  ----   |   ww    -   ww   | 
            !  dz   |i,j,k+1/2       dz    [    k+1        k  ] 


            !dwwdz = 0.5*( dx*dy*w(i,j,k) + dx*dy*w(i,j,k+1) ) * 0.5 * ( w(i,j,k) + w(i,j,k+1) ) - & ! mww, k+1/2
            !        0.5*( dx*dy*w(i,j,k) + dx*dy*w(i,j,k-1) ) * 0.5 * ( w(i,j,k) + w(i,j,k-1) ) ! mww, k-1/2

            !dwwdz = dwwdz / (dx * dy * dzmh )

            dwwdz = ( 0.25*(w(i,j,k) + w(i,j,k+1))**2 - &
                      0.25*(w(i,j,k) + w(i,j,k-1))**2 ) / dzmh

            ! d wu  |                1     [                  ]
            ! ----- |            =  ----   |   wu    -   wu   | 
            !  dx   |i,j,k+1/2       dx    [    k+1        k  ]


            dwudx = 0.5*( dz(k)*u(i+1,j,k) + dz(k-1)*u(i+1,j,k-1) ) * 0.5 * ( w(i,j,k) + w(i+1,j,k) ) - & ! muw, i+1/2
                    0.5*( dz(k)*u(i  ,j,k) + dz(k-1)*u(i  ,j,k-1) ) * 0.5 * ( w(i,j,k) + w(i-1,j,k) )     ! muw, i-1/2

            dwudx = dwudx / (dx * dzmh)

            ! d wv  |                1     [                  ]
            ! ----- |           =   ----   |   wv    -   wv   | 
            !  dy   |i,j,k+1/2       dy    [    k+1        k  ]

            dwvdy = 0.5*( dz(k)*v(i,j+1,k) + dz(k-1)*v(i,j+1,k-1) ) * 0.5 * ( w(i,j,k) + w(i,j+1,k) ) - & ! mvw, j+1/2
                    0.5*( dz(k)*v(i,j  ,k) + dz(k-1)*v(i,j  ,k-1) ) * 0.5 * ( w(i,j,k) + w(i,j-1,k) )     ! mvw, j-1/2
                   
            dwvdy = dwvdy / (dy * dzmh)

            ! d2w / dxj2
            nudd2xy_w = nu *  ( ( w(i-1,j,k) -2.0*w(i,j,k) + w(i+1,j,k)  ) / dx**2  + &
                                ( w(i,j-1,k) -2.0*w(i,j,k) + w(i,j+1,k)  ) / dy**2  )

            nudd2_w = nu  * (   ( w(i,j,k+1) - w(i,j,k  ) ) / dz(k  )  - &
                                ( w(i,j,k  ) - w(i,j,k-1) ) / dz(k-1)   ) / dzmh

            ! dp/dz | k+1/2
            dpdz = ( -p(i,j,k-1) + p(i,j,k) ) / dzmh


            !--------------------- BUILD RHS -----------------------!

            expl_u(i,j,k) = -( duudx + duvdy + duwdz)
            expl_v(i,j,k) = -( dvudx + dvvdy + dvwdz)
            expl_w(i,j,k) = -( dwudx + dwvdy + dwwdz)

            if (implicitXYmode .eqv. .true.) then
                nudd2_u = nudd2_u + nudd2xy_u
                nudd2_v = nudd2_v + nudd2xy_v
                nudd2_w = nudd2_w + nudd2xy_w
            else
                expl_u(i,j,k) = expl_u(i,j,k) + nudd2xy_u
                expl_v(i,j,k) = expl_v(i,j,k) + nudd2xy_v
                expl_w(i,j,k) = expl_w(i,j,k) + nudd2xy_w
            endif

            if (scalarmode .eqv. .true.) then

                Tu_interp = 0.5 * ( temp(i-1,j,k) + temp(i,j,k) )
                Tv_interp = 0.5 * ( temp(i,j-1,k) + temp(i,j,k) )
                Tw_interp = 0.5 * ( temp(i,j,k-1) + temp(i,j,k) )

                expl_u(i,j,k) = expl_u(i,j,k) + ( beta_gx * Tu_interp )
                expl_v(i,j,k) = expl_v(i,j,k) + ( beta_gy * Tv_interp )
                expl_w(i,j,k) = expl_w(i,j,k) + ( beta_gz * Tw_interp )

            endif

            rhs_u(i,j,k) = gamdt*expl_u(i,j,k) + zetdt*expl_u_m1(i,j,k)  & ! Fully-explicit terms (advection + buoyancy)
                          -aldt * dpdx  & ! Pressure gradient
                          +aldt * nudd2_u & ! Diffusive term: semi-implicit part
                          -aldt * mean_dpdx ! Mean pressure gradient body forcing

            rhs_v(i,j,k) = gamdt*expl_v(i,j,k) + zetdt*expl_v_m1(i,j,k) &
                          -aldt * dpdy  &
                          +aldt * nudd2_v ! Diffusive term: semi-implicit part

            rhs_w(i,j,k) = gamdt*expl_w(i,j,k) + zetdt*expl_w_m1(i,j,k)  & ! Fully-explicit terms (advection + buoyancy)
                          -aldt * dpdz  & ! Pressure gradient
                          +aldt * nudd2_w ! Diffusive term: semi-implicit part

            enddo
    enddo
    enddo
    !$omp end parallel do


end subroutine build_rhsVelocity