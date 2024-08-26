	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>initialisation for the shallow water model
    module initialisation
    use numerics_type

    private
    public :: allocate_and_set, add_noise
    contains
	! ALLOCATE AND SET ARRAYS
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>allocate arrays on each PE, and initialise them
	!>@param[inout] ip for this grid
	!>@param[inout] jp for this grid
	!>@param[inout] ntim - number of time-steps
	!>@param[inout] f - rotation rate of planet
	!>@param[inout] re - radius of planet
	!>@param[inout] g - gravity of planet
	!>@param[inout] rho - density
	!>@param[inout] dphi - step in longitude
	!>@param[inout] dtheta - step in latitude
	!>@param[inout] dphin - step in longitude - staggered
	!>@param[inout] dthetan - step in latitude
	!>@param[inout] f_cor - coriolis parameter
	!>@param[inout] h - depth of fluid
	!>@param[inout] hs - height of surface above a reference
	!>@param[inout] u - east-west wind
	!>@param[inout] v - north-south wind
	!>@param[inout] height - total height of fluid
	!>@param[inout] dt - time-step
	!>@param[inout] dx - east-west step
	!>@param[inout] dy - north-south step
	!>@param[inout] x - east-west distance
	!>@param[inout] y - north-south distance
	!>@param[inout] phi - longitude
	!>@param[inout] theta - latitude
	!>@param[inout] phin - staggered longitude
	!>@param[inout] thetan - staggerd latitude
	!>@param[inout] recqdp - for efficiency
	!>@param[inout] recqdp_s - for efficiency
	!>@param[inout] recqdq_s - for efficiency
	!>@param[inout] redq_s - for efficiency
	!>@param[inout] redq - for efficiency
	!>@param[inout] recq - for efficiency
	!>@param[inout] rect - for efficiency
	!>@param[inout] rect_s - for efficiency
	!>@param[inout] cq_s - for efficiency
	!>@param[inout] cq - for efficiency
	!>@param[inout] dp1 - for efficiency
	!>@param[inout] dq - for efficiency
	!>@param[inout] recqdq - for efficiency
	!>@param[inout] u_nudge - wind to nudge to
	!>@param[in] o_halo - number of halos required
	!>@param[inout] ipstart - start of i indexing for this PE
	!>@param[inout] jpstart - start of j indexing for this PE
	!>@param[inout] coords - coordinates of cartesian topology
	!>@param[in] inputfile - netcdf file of saturn winds
	!>@param[in] add_random_height_noise - add noise to get going
	!>@param[in] initially_geostrophic - set to geostrophic balance
	!>@param[in] initial_winds - flag: saturn, or jet?
	!>@param[in] ideal jet parameters: u_jet, jet_noise, theta_jet, h_jet
	!>@param[in] ideal polar vortex parameters: u_vortex, h_vortex
	!>@param[in] perturbation parameters: perturb_strength
	!>@param[in] ip - global ip (all pes)
	!>@param[in] jp - global jp (all pes)
	!>@param[in] wind_factor - factor to multiply mean wind by
	!>@param[in] wind_shift - amount to shift mean wind by
	!>@param[in] wind_reduce - amount to reduce mean wind by
	!>@param[in] runtime - length of run
	!>@param[in] dt_nm - time-step
	!>@param[in] grav - gravity read in from namelist
	!>@param[in] rho_nm - density of fluid read in from namelist
	!>@param[in] re_nm - radius of planet read in from namelist
	!>@param[in] rotation_period_hours - rotation period of planet
	!>@param[in] scale_height - scale height of atmosphere
	!>@param[inout] slat - most southerly latitude
	!>@param[inout] nlat - most northerly latitude
	!>@param[in] slat_thresh - most southerly latitude
	!>@param[in] nlat_thresh - most northerly latitude
	!>@param[in] dims - dimensions of cartesian topology
	!>@param[in] id - id of this PE
	!>@param[in] comm2d - communicator for cartesian topology
	!>@param[in] polar_vortex: whether or not to include a polar vortex
	subroutine allocate_and_set(ipp,jpp,ntim, f, &
				re, g, rho, dphi, dtheta, dphin,dthetan, &
				f_cor,h,hs, u, v, &
				height, dt, dx, dy, x, y, &
				phi, theta, phin, thetan, &
				recqdp, recqdp_s, recqdq_s, redq_s,redq, &
				recq, rect, rect_s, cq_s, cq, dp1, dq, &
				recqdq, &
				u_nudge, o_halo, ipstart, jpstart, coords, &
				inputfile, add_random_height_noise, &
				initially_geostrophic, initial_winds, &
				u_jet, u_vortex, jet_noise, theta_jet, h_jet, &
				h_vortex, perturb_strength, ip, jp, &
				wind_factor, wind_shift, wind_reduce, runtime, &
				dt_nm, grav, rho_nm, re_nm, &
				rotation_period_hours, scale_height, slat, nlat, &
				slat_thresh, nlat_thresh, dims, id, comm2d, polar_vortex)
				
		use numerics_type
		use mpi
		use netcdf
		use numerics, only : find_pos, poly_int
		use random, only : random_normal
		use mpi_module
		
		implicit none
		! grid variables to be set
		integer(i4b), intent(inout) :: ipp, jpp, ntim, o_halo
    	real(wp), intent(inout) :: f, re, g, rho, dt
    	real(wp), intent(inout), allocatable, dimension(:,:) :: &
    									f_cor, h, hs, u, v, height, &
    											dx, dy, x, y, &
    				recqdp, recqdp_s, recqdq_s, redq_s, redq, &
    				recq, rect, rect_s, cq_s, cq, dp1, dq, recqdq
    	real(wp), intent(inout), allocatable, dimension(:) :: &
    									phi, theta, phin, thetan,u_nudge, &
    									dphi, dtheta, dphin, dthetan

		! namelist variables used to set the grid
		character (len=*), intent(in) :: inputfile
		logical, intent(in) :: add_random_height_noise, initially_geostrophic, polar_vortex
		integer(i4b), intent(in) :: initial_winds, ip, jp
		real(wp), intent(in) :: wind_factor, wind_shift, wind_reduce, runtime, &
							dt_nm, grav, rho_nm, re_nm, &
							rotation_period_hours, scale_height, &
							slat_thresh, nlat_thresh, &
							u_jet, u_vortex, jet_noise, theta_jet, h_jet, &
							h_vortex, perturb_strength
		real(wp), intent(inout) :: slat, nlat
		integer(i4b), dimension(2), intent(in) :: dims
		integer(i4b), dimension(2), intent(inout) :: coords
		integer(i4b), intent(inout) :: ipstart, jpstart
		integer(i4b), intent(in) :: id, comm2d
		! locals
		integer(i4b) :: iloc, error, AllocateStatus, ncid, varid1,varid2, dimid, nlats, &
						i, j
		real(wp), dimension(:), allocatable :: latitude, wind
		real(wp) :: var, dummy, delta_omega, slat_thresh2, nlat_thresh2
		! for random number:
		real(wp) :: r
		real(wp), dimension(10,10) :: rs
		integer(i4b) :: k, nbottom, ntop, tag1
		integer(i4b), allocatable, dimension(:) :: seed
		

		if(id>=dims(1)*dims(2)) return

		! start of initialisation code
		! 1. basic equalities:
		re=re_nm
		g=grav
		! rho=rho_nm
		dt=dt_nm

		! 2. scalar formulae:
		ntim=ceiling(runtime/dt)
		f=2._wp*PI / (rotation_period_hours*3600._wp)*2707788_wp
		
		! 3. arrays:


		! FIND THE NUMBER OF GRID POINTS ON EACH PE
		call MPI_CART_COORDS(comm2d, id, 2, coords, error)
		! print *,'Coords of ',id,' are ',coords

		! number of grid points in all but last:
		ipp = floor(real(ip,wp)/real(dims(1),wp)) 
		ipstart = ipp*(coords(1))     
		if(coords(1) == (dims(1)-1)) then
			ipp=ip-(dims(1)-1)*ipp ! number of grid points in last
		endif
		! number of grid points in all but last:
		jpp = floor(real(jp,wp)/real(dims(2),wp))      
		jpstart = jpp*(coords(2))     
		if(coords(2) == (dims(2)-1)) then
			jpp=jp-(dims(2)-1)*jpp ! number of grid points in last
		endif
		! print *,ip,jp,ipp,jpp,ipstart, jpstart,coords


		! ALLOCATE ARRAYS
		allocate( f_cor(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( h(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( hs(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( u(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( v(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( height(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( x(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( y(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dx(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dy(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		
		allocate( recqdp(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( recqdp_s(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( recqdq_s(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( redq_s(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( redq(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( recq(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( rect(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( rect_s(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( cq_s(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( cq(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dp1(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dq(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( recqdq(1-o_halo:ipp+o_halo,1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	

		allocate( phi(1-o_halo:ipp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dphi(1-o_halo:ipp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( phin(1-o_halo:ipp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dphin(1-o_halo:ipp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( theta(1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dtheta(1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( thetan(1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( dthetan(1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( u_nudge(1-o_halo:jpp+o_halo), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"


		! READ NETCDF FILE
		! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
		! the file.
		call check( nf90_open(inputfile, NF90_NOWRITE, ncid) )

		call check( nf90_inq_dimid(ncid, "nlats", dimid) )
		call check( nf90_inquire_dimension(ncid, dimid, len = nlats) )


		! ALLOCATE ARRAYS
		allocate( latitude(1:nlats), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( wind(1:nlats), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"


		! Get the varid of the latitude variable, based on its name.
		call check( nf90_inq_varid(ncid, "latitude", varid1) )
		! Get the varid of the wind variable, based on its name.
		call check( nf90_inq_varid(ncid, "wind", varid2) )
		
		call check( nf90_get_var(ncid, varid1, latitude, start = [1] ) )
		call check( nf90_get_var(ncid, varid2, wind, start = [1] ) )
		
		! Close the file, freeing all resources.
		call check( nf90_close(ncid) )
		


		! interpolate to grid, do finite diffs, pass halos, etc
		dphi=(2._wp*PI) / real(ip-1,wp) ! lon
		dphin=(2._wp*PI) / real(ip-1,wp) ! lon

		! set up longitude array:
		phi=dphi(1)*(/(i,i=ipstart+1-o_halo-1,ipstart+ipp+o_halo-1)/)
		phin=phi+dphi/2._wp

		nlat=min(nlat,nlat_thresh)
		slat=max(slat,slat_thresh)
			
		! latitude array:	
		dtheta=(nlat-slat) &
				/ real(jp-1,wp) * PI/180._wp  ! lat
		
		! ! deal with singularity at the poles:
		! if((coords(2) == (dims(2)-1))  .and. (nlat .gt. 0._wp)) then
		! 	dtheta(jpp) = 2._wp*(90._wp-nlat)*PI/180._wp
		! endif
		! if((coords(2) == 0)  .and. (slat .lt. 0._wp)) then
		! 	dtheta(0) = 2._wp*(90._wp+slat)*PI/180._wp
		! endif

		
		! SET UP LATITUDE ARRAY NEED MPI TO ADD THEM UP:
		call MPI_CART_SHIFT( comm2d, 1, 1, nbottom, ntop, error)
		if(nbottom .ne. -1) then
			tag1=001
			call MPI_Recv(theta(1-o_halo), 1, MPI_REAL8, nbottom, &
					tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE,error)
		endif

		if(coords(2) == 0) then
			theta(0) = slat*PI/180._wp-dtheta(0)
		endif

		do i=1,jpp+o_halo
			theta(i)=theta(i-1)+dtheta(i-1)
		enddo
		
		if(ntop .ne. -1) then
			tag1=001
			call MPI_Send(theta(jpp), 1, MPI_REAL8, ntop, &
				tag1, MPI_COMM_WORLD,error)
		endif
		
		
		thetan=theta+dtheta/2._wp
		do i=1-o_halo,jpp
			dthetan(i)=thetan(i+1)-thetan(i)
		enddo
		
		! if this is the top then set dtheta:
		if(ntop == -1 ) then
			dthetan(jpp+1)=2._wp*(90._wp*PI/180._wp -dthetan(jpp))
		endif 
		if(ntop /= -1 ) then
			call MPI_Recv(dthetan(jpp+o_halo), 1, MPI_REAL8, ntop, &
					tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE,error)
		endif 
		if(nbottom /= -1 ) then
			call MPI_Send(dthetan(1), 1, MPI_REAL8, nbottom, &
					tag1, MPI_COMM_WORLD,error)
		endif 	
		! dthetan=(nlat-slat) / real(jp-1,wp) * PI/180._wp  ! lat


		select case (initial_winds)    
			case (1)  ! saturn winds
			! interpolate winds to u_nudge
			delta_omega=2._wp*PI/(3600._wp)*(1._wp/10.656_wp-1._wp/rotation_period_hours)
			do i=1-o_halo,jpp+o_halo
				iloc=find_pos(latitude(1:nlats),asin(sin(theta(i)))*180._wp/PI)
				iloc=min(nlats-1,iloc)
				iloc=max(1,iloc)
				! linear interp theta

				call poly_int(latitude(iloc:iloc+1), wind(iloc:iloc+1), &
					min(max( asin(sin(theta(i)))*180._wp/PI, &
					latitude(nlats)),latitude(1)), var,dummy)
				
				if((theta(i)*180._wp/pi>90._wp) .or. &
					(theta(i)*180._wp/pi<-90._wp) ) var=-var

				var=var+delta_omega*re*cos(theta(i))
			
				u_nudge(i)=var
			enddo

			case (2) ! ideal jet:
			do i=1-o_halo,jpp+o_halo
				if (polar_vortex) then
					u_nudge(i)=u_jet* exp(-perturb_strength*(theta(i)-theta_jet*pi/180._wp)**2._wp &
									/ (h_jet*pi/180._wp)**2._wp ) + &
									u_vortex*exp(-perturb_strength*(theta(i)-(theta_jet+12.6_wp)*pi/180._wp)**2._wp &
									/ (h_vortex*pi/180._wp)**2._wp )
				else
					u_nudge(i)=u_jet*exp(-perturb_strength*(theta(i)-theta_jet*pi/180._wp)**2._wp &
					/ (h_jet*pi/180._wp)**2._wp )
				endif
			enddo

			case default
				print *,'error initial_winds',initial_winds
				stop
		end select
		
		! calculate Coriolis param:
		do i=1-o_halo,ipp+o_halo
			f_cor(i,:)=2._wp*f*sin(theta)
		enddo
			
		! calculate x, y, dx, dy:
		do i=1-o_halo,ipp+o_halo
			x(i,:)=phi(i)*cos(asin(sin(theta)))*re			
		enddo
		do i=1-o_halo,ipp+o_halo
			y(i,:)=(theta)*re			
		enddo
		
		do i=1-o_halo,ipp+o_halo
			dx(i,:)=dphi(i)*cos(asin(sin(theta)))*re			
		enddo
		do i=1-o_halo,ipp+o_halo
			dy(i,:)=(dtheta)*re			
		enddo
		
		
		! set u and v:
		v(:,:)=0._wp
		do j=1-o_halo,jpp+o_halo
			u(:,j)=u_nudge(j)
		enddo
		
		! set surface to zero.
		hs(:,:)=0._wp



		! CALCULATE THE HEIGHT FIELD FROM WINDS - MPI NEEDED TO SPAN SUB-DOMAINS
		call MPI_CART_COORDS(comm2d, id, 2, coords, error)
		call MPI_CART_SHIFT( comm2d, 1, 1, nbottom, ntop, error)	
		
		! if the y coordinate is not the most northerly
		! then receive from domain north of here:
		if ( (coords(2)+1) /= dims(2) ) then
			tag1=2010
			call MPI_Recv(height(:,jpp+1:jpp+o_halo),& ! the data packet to receive into
				ipp+2*o_halo, & ! size of the data packet
				MPI_REAL8, &
				ntop, & ! receive from above (north)
				tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE,error)
		
		endif		
		! if most northerly grid point (note set the jpp+1 point)
		if( (coords(2)+1) == dims(2) ) then
			height(:,jpp+1)=scale_height
		endif

		do j=jpp,0,-1
			height(:,j)=height(:,j+1)+ &	
					0.25_wp*(f_cor(:,j+1)+f_cor(:,j))* &
					(u(:,j+1)+u(:,j))* &
					re/g*(dtheta(j))
		enddo				
		! if the y coordinate is not most southerly
		if ( (coords(2)+1) /= 1 ) then
			tag1=2010
			call MPI_Send(height(:,1:1),& ! the data packet to send
				ipp+2*o_halo, & ! size of the data packet
				MPI_REAL8, &
				nbottom, & ! send to below (south)
				tag1, MPI_COMM_WORLD,error)
		
		endif
		
		! Note: Currently commented out due to adding a circular band of zeros on the v
		! velocity. Not certain what causes this - the code within the method is essentially
		! the same as the method and it works fine.
 		! ! Call the subroutine to add noise to the height field
		! call add_noise(height, initial_winds, ip, jp, ipstart, jpstart, ipp, jpp, theta, jet_noise, polar_vortex, theta_jet, h_jet)		

		! CALCULATE AND ADD NOISE	
		call random_seed(size=k)
		allocate(seed(1:k))
		seed(:)=2
		call random_seed(put=seed)
		select case (initial_winds)
			case (1) ! saturn winds
			do j=1,jp
				do i=1,ip
					r=random_normal() ! from the Netlib

					if((i > ipstart) .and. (i <=ipstart+ipp) &
						.and. (j > jpstart) .and. (j <= jpstart+jpp) ) then
					
						if ((theta(j-jpstart)*180._wp/PI) > 75._wp &
							.and. (theta(j-jpstart)*180._wp/PI) < 80._wp) then
						
							height(i-ipstart,j-jpstart) = &
								height(i-ipstart,j-jpstart) + &
								r*1000_wp*0.6e5_wp/height(i-ipstart,j-jpstart)
						endif
					endif

				enddo
			enddo

			case (2) ! ideal jet
			do j=1,jp
				do i=1,ip
					r=random_normal() ! from the Netlib

					if((i > ipstart) .and. (i <=ipstart+ipp) &
						.and. (j > jpstart) .and. (j <= jpstart+jpp) ) then
					
						if ((theta(j-jpstart)*180._wp/PI) > (theta_jet-h_jet*3._wp) &
						.and. (theta(j-jpstart)*180._wp/PI) <(theta_jet+h_jet*3._wp)) then
						
							height(i-ipstart,j-jpstart) = &
								height(i-ipstart,j-jpstart) + &
								 r*jet_noise/height(i-ipstart,j-jpstart) 
						endif
						if (polar_vortex .and. (theta(j-jpstart)*180._wp/PI) > (theta_jet+12.6-h_vortex*2._wp) &
						.and. (theta(j-jpstart)*180._wp/PI) <(theta_jet+12.6+h_vortex*2._wp)) then
						
							height(i-ipstart,j-jpstart) = &
								height(i-ipstart,j-jpstart) + &
								 r*jet_noise/height(i-ipstart,j-jpstart)
						endif
					endif

				enddo
			enddo

			case default
				print *,'error initial_winds',initial_winds
				stop
		end select
		
		deallocate(seed)

		! set halos in Coriolis array
		call exchange_halos(comm2d, id, ipp, jpp, o_halo, f_cor)
		! set halos in height array
		call exchange_halos(comm2d, id, ipp, jpp, o_halo, height)

	
		! CALCULATE U AND V FROM HEIGHT FIELD (WITH NOISE) - NO MESSAGE PASSING NEEDED
		! IF HALOS IN HEIGHT ARE SET CORRECTLY	
		if(initially_geostrophic) then
			do j=1,jpp
				u(:,j)=-g*(height(:,j+1)-height(:,j-1)) / &
						(re*(dtheta(j)+dtheta(j-1))*f_cor(:,j))
			enddo
			do i=1,ipp
				v(i,:)=g*(height(i+1,:)-height(i-1,:)) / &
						(re*(dphi(i)+dphi(i-1))*cos(theta(:))*f_cor(i,:))
			enddo		
		endif	


		! set halos in u array
		call exchange_halos(comm2d, id, ipp, jpp, o_halo, u)
		! set halos in v array
		call exchange_halos(comm2d, id, ipp, jpp, o_halo, v)


		! set h array
		h(:,:)=height(:,:)-hs(:,:)
		
		
		! set some variables for efficiency
		do j=1-o_halo,jpp+o_halo
			recqdp(:,j)=re*cos(theta(j))*dphi(:)
			recqdp_s(:,j)=re*cos(theta(j))*dphin(:)
			recqdq_s(:,j)=re*cos(thetan(j))*dthetan(j)
			redq_s(:,j)=re*dthetan(j)
			redq(:,j)=re*dtheta(j)

			recq(:,j)=re*cos(theta(j))
			recqdq(:,j)=re*cos(theta(j))*dtheta(j)
			cq_s(:,j)=cos(thetan(j))
			cq(:,j)=cos(theta(j))
			dp1(:,j)=dphi(:)
			dq(:,j)=dtheta(j)

			rect(:,j)=tan(theta(j))/re
			rect_s(:,j)=tan(thetan(j))/re
		enddo		

		deallocate(wind)
		deallocate(latitude)
		
	end subroutine allocate_and_set
	

	! HELPER ROUTINE
	subroutine check(status)
	use netcdf
	use numerics_type
	integer(i4b), intent ( in) :: status

	if(status /= nf90_noerr) then
		print *, trim(nf90_strerror(status))
		stop "Stopped"
	end if
	end subroutine check

	! ADD NOISE ROUTINE
	subroutine add_noise(height, initial_winds, ip, jp, ipstart, jpstart, ipp, jpp, theta, jet_noise, polar_vortex, theta_jet, h_jet)
		use random, only : random_normal
		implicit none
		
		! Define types and variables
		integer, intent(in) :: initial_winds, ip, jp, ipstart, jpstart, ipp, jpp
		real(wp), dimension(ipp,jpp), intent(inout) :: height
		real(wp), intent(in) :: theta(:), jet_noise, theta_jet, h_jet
		logical, intent(in) :: polar_vortex

		! locals
		integer(i4b) :: i, j
	
		! for random number:
		real(wp) :: r
		integer(i4b) :: k
		integer(i4b), allocatable, dimension(:) :: seed
		
		! Initialize the random number generator with a fixed seed
		call random_seed(size=k)
		allocate(seed(1:k))
		seed(:)=2
		call random_seed(put=seed)
	
		! Calculate and add noise based on the initial_winds parameter
		select case (initial_winds)
			case (1) ! Saturn winds
				do j = 1, jp
					do i = 1, ip
						r=random_normal() ! from the Netlib
						if((i > ipstart) .and. (i <=ipstart+ipp) &
							.and. (j > jpstart) .and. (j <= jpstart+jpp) ) then
							if ((theta(j-jpstart)*180._wp/PI) > 75._wp &
								.and. (theta(j-jpstart)*180._wp/PI) < 80._wp) then
								height(i-ipstart,j-jpstart) = &
									height(i-ipstart,j-jpstart) + &
									r*1000_wp*0.6e5_wp/height(i-ipstart,j-jpstart)
							endif
						endif
					enddo
				enddo
				
			case (2) ! Ideal jet
				do j=1,jp
					do i=1,ip
						r=random_normal() ! from the Netlib

						if((i > ipstart) .and. (i <=ipstart+ipp) &
							.and. (j > jpstart) .and. (j <= jpstart+jpp) ) then
						
							if ((theta(j-jpstart)*180._wp/PI) > (theta_jet-h_jet*3._wp) &
							.and. (theta(j-jpstart)*180._wp/PI) <(theta_jet+h_jet*3._wp)) then
							
								height(i-ipstart,j-jpstart) = &
									height(i-ipstart,j-jpstart) + &
									 r*jet_noise/height(i-ipstart,j-jpstart) 
							endif
							if (polar_vortex .and. (theta(j-jpstart)*180._wp/PI) > (theta_jet+12.6-h_jet*2._wp) &
							.and. (theta(j-jpstart)*180._wp/PI) <(theta_jet+12.6+h_jet*2._wp)) then
							
								height(i-ipstart,j-jpstart) = &
									height(i-ipstart,j-jpstart) + &
									 r*jet_noise/height(i-ipstart,j-jpstart)
							endif
						endif

					enddo
				enddo
				
			case default
				print *, 'error initial_winds', initial_winds
				stop
		end select
		deallocate(seed)
	end subroutine add_noise	
	

	end module initialisation
	