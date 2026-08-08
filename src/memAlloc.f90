! Allocate memory arrays: fields and working memory arrays
subroutine allocFields
    use parameters
    use stats
    use ghost
    use velfields
    use velMemory
    use scalarfields
    use scalarMemory
    use implicit
    use implicit_types
    implicit none


    ! Fields
    allocate( u( -halosize+1:Nx+halosize , -halosize+1:Ny+halosize , -halosize+1:Nz+halosize )  )
    allocate( v( -halosize+1:Nx+halosize , -halosize+1:Ny+halosize , -halosize+1:Nz+halosize )  )
    allocate( w( -halosize+1:Nx+halosize , -halosize+1:Ny+halosize , -halosize+1:Nz+halosize )  )
    allocate( p( -halosize+1:Nx+halosize , -halosize+1:Ny+halosize , -halosize+1:Nz+halosize )  )

    allocate( pseudo_p( -halosize+1:Nx+halosize , -halosize+1:Ny+halosize , -halosize+1:Nz+halosize )  )

    ! Working memory: RHS vectors
    allocate( rhs_u( 1:Nx , 1:Ny , 1:Nz )  )
    allocate( rhs_v( 1:Nx , 1:Ny , 1:Nz )  )
    allocate( rhs_w( 1:Nx , 1:Ny , 1:Nz )  )

    ! RK3
    allocate( expl_u( 1:Nx , 1:Ny , 1:Nz ) ) ; allocate( expl_u_m1( 1:Nx , 1:Ny , 1:Nz ) ) 
    allocate( expl_v( 1:Nx , 1:Ny , 1:Nz ) ) ; allocate( expl_v_m1( 1:Nx , 1:Ny , 1:Nz ) ) 
    allocate( expl_w( 1:Nx , 1:Ny , 1:Nz ) ) ; allocate( expl_w_m1( 1:Nx , 1:Ny , 1:Nz ) ) 

    ! Implicit solver working memory
    allocate(ami( 1:Nx )) ; allocate(aci( 1:Nx )) ; allocate(api( 1:Nx ))
    allocate(amj( 1:Ny )) ; allocate(acj( 1:Ny )) ; allocate(apj( 1:Ny ))
    allocate(amk( 1:Nz )) ; allocate(ack( 1:Nz )) ; allocate(apk( 1:Nz ))
    allocate(tdm_rhsX1( 1:Nx )) ; allocate(tdm_rhsX2( 1:Nx )) 
    allocate(tdm_rhsY1( 1:Ny )) ; allocate(tdm_rhsY2( 1:Ny )) 
    allocate(tdm_rhsZ_r( 1:Nz )) ; allocate(tdm_rhsZ_c( 1:Nz ))

    if ( (implicitXYmode .eqv. .true. ) .and. (implicit_type .eq. ADI) ) then
        allocate(impl_delta( 1:Nx, 1:Ny, 1:Nz ))
    endif

    if (scalarmode .eqv. .true.) then
        allocate( temp( -halosize+1:Nx+halosize , -halosize+1:Ny+halosize , -halosize+1:Nz+halosize )  )
        allocate( rhs_temp( 1:Nx , 1:Ny , 1:Nz )  )
        allocate( expl_c( 1:Nx , 1:Ny , 1:Nz ) ) ; allocate( expl_c_m1( 1:Nx , 1:Ny , 1:Nz ) ) 
    endif

    ! Means
    allocate( uavg(1:Nz) )
    allocate( vavg(1:Nz) )
    allocate( wavg(1:Nz) )

    allocate( uu_avg(1:Nz) )
    allocate( vv_avg(1:Nz) )
    allocate( ww_avg(1:Nz) )
    allocate( uw_avg(1:Nz) )

    

    if (scalarmode ) then
        allocate( cavg(1:Nz) )
        allocate( cc_avg(1:Nz) )
        allocate( cw_avg(1:Nz) )
    endif

end subroutine allocFields


subroutine deallocFields
    use stats
    use velfields
    use velMemory
    use scalarfields
    use scalarMemory
    use implicit
    implicit none

    if(allocated(u)) deallocate(u)
    if(allocated(v)) deallocate(v)
    if(allocated(w)) deallocate(w)
    if(allocated(p)) deallocate(p)
    if(allocated(pseudo_p)) deallocate(pseudo_p)

    if(allocated(rhs_u)) deallocate(rhs_u)
    if(allocated(rhs_v)) deallocate(rhs_v)
    if(allocated(rhs_w)) deallocate(rhs_w)
    if(allocated(expl_u)) deallocate(expl_u) ; if(allocated(expl_u_m1)) deallocate(expl_u_m1) 
    if(allocated(expl_v)) deallocate(expl_v) ; if(allocated(expl_v_m1)) deallocate(expl_v_m1) 
    if(allocated(expl_w)) deallocate(expl_w) ; if(allocated(expl_w_m1)) deallocate(expl_w_m1) 

    if(allocated(ami)) deallocate(ami) ; if(allocated(aci)) deallocate(aci) ; if(allocated(api)) deallocate(api)
    if(allocated(amj)) deallocate(amj) ; if(allocated(acj)) deallocate(acj) ; if(allocated(apj)) deallocate(apj)
    if(allocated(amk)) deallocate(amk) ; if(allocated(ack)) deallocate(ack) ; if(allocated(apk)) deallocate(apk)
    if(allocated(tdm_rhsX1)) deallocate(tdm_rhsX1) ; if(allocated(tdm_rhsX2)) deallocate(tdm_rhsX2) 
    if(allocated(tdm_rhsY1)) deallocate(tdm_rhsY1) ; if(allocated(tdm_rhsY2)) deallocate(tdm_rhsY2) 
    if(allocated(tdm_rhsZ_r)) deallocate(tdm_rhsZ_r) ;     if(allocated(tdm_rhsZ_c)) deallocate(tdm_rhsZ_c)
    if(allocated(impl_delta)) deallocate(impl_delta)

    if(allocated(temp)) deallocate(temp)
    if(allocated(rhs_temp)) deallocate(rhs_temp)
    if(allocated(expl_c)) deallocate(expl_c) ; if(allocated(expl_c_m1)) deallocate(expl_c_m1) 

    if(allocated(uavg)) deallocate(uavg) ; if(allocated(vavg)) deallocate(vavg) ; if(allocated(wavg)) deallocate(wavg) 

    if(allocated(uu_avg)) deallocate(uu_avg) ; if(allocated(vv_avg)) deallocate(vv_avg) ; if(allocated(ww_avg)) deallocate(ww_avg)
    if(allocated(uw_avg)) deallocate(uw_avg)

    if(allocated(cavg)) deallocate(cavg)
    if(allocated(cc_avg)) deallocate(cc_avg) ; if(allocated(cw_avg)) deallocate(cw_avg)

end subroutine deallocFields
