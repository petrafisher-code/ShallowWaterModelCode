	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>advection routines
    module advection
    use numerics_type
    private
    public :: lax_wendroff_conservative, lax_wendroff_ll, dissipation, smagorinsky
    contains
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>advects a scalar field on a ll grid 
	!>@param[in] ip: number of east-west points
	!>@param[in] jp: ditto for north-south
	!>@param[in] o_halo: halos required for advection scheme
	!>@param[in] dt:  timestep
	!>@param[in] g: gravity
	!>@param[inout] u,v,h: prognostic variables
	!>@param[in] hs: surface height
	!>@param[in] re: radius of planet
	!>@param[in] theta: latitude
	!>@param[in] thetan: latitude - staggered
	!>@param[in] dtheta: latitude step
	!>@param[in] dthetan: latitude step - staggered
	!>@param[in] phi: phi
	!>@param[in] phin: phin
	!>@param[in] dphi: dphi
	!>@param[in] dphin: dphin
	!>@param[in] f_cor: Coriolis parameter
	!>@param[in] recqdq - for efficiency
	!>@param[in] recqdp - for efficiency
	!>@param[in] recqdp_s - for efficiency
	!>@param[in] recqdq_s - for efficiency
	!>@param[in] redq_s - for efficiency
	!>@param[in] redq - for efficiency
	!>@param[in] rect - for efficiency
	!>@param[in] rect_s - for efficiency
	!>@param[in] cq - for efficiency
	!>@param[in] cq_s - for efficiency
	!>solves the 1-d advection equation:
	!>\f$ \frac{\partial \psi}{\partial t} + \frac{\partial u \psi}{\partial x} = 0 \f$
    subroutine lax_wendroff_conservative(ip,jp,o_halo,dt,g,u,v,h,hs,re,&
		theta,thetan,dtheta,dthetan, phi, phin, dphi, dphin, f_cor, &
		recqdq, recqdp, recqdp_s, recqdq_s, redq_s, redq, rect, rect_s, cq, cq_s)

		use numerics_type
		implicit none
		integer(i4b), intent(in) :: ip,jp,o_halo
		real(wp), intent(in) :: dt, g, re
		real(wp), intent(in), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: &
																		hs, f_cor, &
					recqdq, recqdp, recqdp_s, recqdq_s, redq_s, redq, rect, &
					rect_s, cq, cq_s
		real(wp), dimension(1-o_halo:jp+o_halo), intent(in) :: &
											theta,thetan, dtheta, dthetan
		real(wp), dimension(1-o_halo:ip+o_halo), intent(in) :: &											
											phi, phin, dphi, dphin
		real(wp), intent(inout), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: &
																				h, u, v
																				
		! local variables:
		real(wp), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: & 
					dy1, v1, h1, vh, uh, vh1, Ux, Uy, Vx, Vy, Vy2, uv, u2
		real(wp), dimension(1:ip,1:jp) :: &
										uh_new, vh_new, h_new
										
		real(wp), dimension(0:ip,1:jp) :: h_mid_xt, uh_mid_xt, Ux_mid_xt, Vx_mid_xt
		real(wp), dimension(0:ip,1:jp) :: vh_mid_xt
		real(wp), dimension(1:ip,0:jp) :: h_mid_yt, uh_mid_yt, vh_mid_yt,  &
												Uy_mid_yt, Vy_mid_yt, Vy_mid_yt2
		integer(i4b) :: j, i
			
		
		v1=v*cq ! cq=cos(theta)
		h1=h*cq
		
		uh=u *h
		vh=v*h
		vh1=v1*h
		
		! continuity equation (calculate mid-point values at 0.5*dt):
		h_mid_xt = 0.5_wp*(h(1:ip+1,1:jp)+h(0:ip,1:jp)) &
		-(0.5_wp*dt/(recqdp_s(0:ip,1:jp))) &
		*(uh(1:ip+1,1:jp)-uh(0:ip,1:jp))
		
		h_mid_yt = 0.5_wp*(h(1:ip,1:jp+1)+h(1:ip,0:jp)) &
		-(0.5_wp*dt/(recqdq_s(1:ip,0:jp))) &
		*(vh1(1:ip,1:jp+1)-vh1(1:ip,0:jp))

		! v-phi, or u momentum equation (calculate mid-point values at 0.5*dt):
		Ux = uh*u + g*h**2*(0.5_wp)
		Uy = uh*v1
		uh_mid_xt(0:ip,:) = 0.5_wp*(uh(1:ip+1,1:jp)+uh(0:ip,1:jp)) &
		-(0.5_wp*dt/(recqdp_s(0:ip,1:jp)))* &
			(Ux(1:ip+1,1:jp)-Ux(0:ip,1:jp)) &
		+0.125_wp*dt*(f_cor(1:ip+1,1:jp)+f_cor(0:ip,1:jp))* &
			(vh(1:ip+1,1:jp)+vh(0:ip,1:jp))
		
		uh_mid_yt = 0.5_wp*(uh(1:ip,1:jp+1)+uh(1:ip,0:jp)) &
		-(0.5_wp*dt/(recqdq_s(1:ip,0:jp)))* &
			(Uy(1:ip,1:jp+1)-Uy(1:ip,0:jp)) &
		+0.125_wp*dt*(f_cor(1:ip,1:jp+1)+f_cor(1:ip,0:jp))* &
			(vh(1:ip,1:jp+1)+vh(1:ip,0:jp))

		! v-theta, or v momentum equation (calculate mid-point values at 0.5*dt):
		Vx = uh*v
		Vy = vh1*v
		Vy2 = 0.5_wp*g*h**2
		vh_mid_xt(0:ip,1:jp) = 0.5_wp*(vh(1:ip+1,1:jp)+vh(0:ip,1:jp)) &
		-(0.5_wp*dt/(recqdp_s(0:ip,1:jp)))*(Vx(1:ip+1,1:jp)-Vx(0:ip,1:jp)) &
		-0.125_wp*dt*(f_cor(1:ip+1,1:jp)+f_cor(0:ip,1:jp))*(uh(1:ip+1,1:jp)+uh(0:ip,1:jp))

		vh_mid_yt(1:ip,0:jp) = 0.5_wp*(vh(1:ip,1:jp+1)+vh(1:ip,0:jp)) &
		-(0.5_wp*dt/(recqdq_s(1:ip,0:jp)))*(Vy(1:ip,1:jp+1)-Vy(1:ip,0:jp)) &
		-(0.5_wp*dt/(redq_s(1:ip,0:jp)))*(Vy2(1:ip,1:jp+1)-Vy2(1:ip,0:jp)) &
		-0.125_wp*dt*(f_cor(1:ip,1:jp+1)+f_cor(1:ip,0:jp))*(uh(1:ip,1:jp+1)+uh(1:ip,0:jp))

		! Now use the mid-point values to predict the values at the next timestep
		! continuity:
		h_new = h(1:ip,1:jp) &
		- (dt/(recqdp(0:ip-1,1:jp)))*(uh_mid_xt(1:ip,1:jp)-uh_mid_xt(0:ip-1,1:jp)) &
		- (dt/(recqdq(1:ip,0:jp-1))) * &
		(vh_mid_yt(1:ip,1:jp)*cq_s(1:ip,1:jp)-vh_mid_yt(1:ip,0:jp-1)*cq_s(1:ip,0:jp-1))

		! u-momentum equation:
		Ux_mid_xt = uh_mid_xt*uh_mid_xt/h_mid_xt + 0.5_wp*g*h_mid_xt**2
		Uy_mid_yt = uh_mid_yt*vh_mid_yt/h_mid_yt*cq_s(1:ip,0:jp)
		uh_new = uh(1:ip,1:jp) &
		- (dt/(recqdp(0:ip-1,1:jp)))*  (Ux_mid_xt(1:ip,1:jp)-Ux_mid_xt(0:ip-1,1:jp)) &
		- (dt/(recqdq(1:ip,0:jp-1)))*(Uy_mid_yt(1:ip,1:jp)-Uy_mid_yt(1:ip,0:jp-1))

		! v-momentum equation:
		Vx_mid_xt = uh_mid_xt*vh_mid_xt/h_mid_xt
		Vy_mid_yt = vh_mid_yt*vh_mid_yt/h_mid_yt*cq_s(1:ip,0:jp)
		Vy_mid_yt2 = 0.5_wp*g*h_mid_yt**2
		vh_new = vh(1:ip,1:jp) &
		- (dt/(recqdp(0:ip-1,1:jp)))*(Vx_mid_xt(1:ip,1:jp)-Vx_mid_xt(0:ip-1,1:jp)) &
		- (dt/(recqdq(1:ip,0:jp-1)))*(Vy_mid_yt(1:ip,1:jp)-Vy_mid_yt(1:ip,0:jp-1)) &
		- (dt/(redq(1:ip,0:jp-1) ))* &
		(Vy_mid_yt2(1:ip,1:jp)-Vy_mid_yt2(1:ip,0:jp-1))

		! add on Coriolis and contribution of orography to pressure gradient:
		uv = u*v
		u2 = u*u
		uh_new=uh_new + dt*.5_wp*(f_cor(1:ip,1:jp)*v(1:ip,1:jp) & 
			- uv(1:ip,1:jp)*rect(1:ip,1:jp))*(h(1:ip,1:jp)+h_new)

		vh_new=vh_new - dt*.5_wp*(f_cor(1:ip,1:jp)*u(1:ip,1:jp) & 
			+ u2(1:ip,1:jp)*rect(1:ip,1:jp))*(h(1:ip,1:jp)+h_new)

		! re-calculate u and v.
		u(1:ip,1:jp) = uh_new/(h_new)
		v(1:ip,1:jp) = vh_new/h_new
		h(1:ip,1:jp) = h_new

	end subroutine lax_wendroff_conservative
	
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>advects a scalar field on a ll grid 
	!>@param[in] ip: number of east-west points
	!>@param[in] jp: ditto for north-south
	!>@param[in] o_halo: halos required for advection scheme
	!>@param[in] dt:  timestep
	!>@param[in] g: gravity
	!>@param[inout] u,v,h: prognostic variables
	!>@param[in] hs: surface height
	!>@param[in] re: radius of planet
	!>@param[in] theta: latitude
	!>@param[in] thetan: latitude - staggered
	!>@param[in] dtheta: latitude step
	!>@param[in] dthetan: latitude step - staggered
	!>@param[in] phi: phi
	!>@param[in] phin: phin
	!>@param[in] dphi: dphi
	!>@param[in] dphin: dphin
	!>@param[in] f_cor: Coriolis parameter
	!>@param[in] recqdq - for efficiency
	!>@param[in] recqdp - for efficiency
	!>@param[in] recqdp_s - for efficiency
	!>@param[in] recqdq_s - for efficiency
	!>@param[in] redq_s - for efficiency
	!>@param[in] redq - for efficiency
	!>@param[in] rect - for efficiency
	!>@param[in] rect_s - for efficiency
	!>@param[in] cq - for efficiency
	!>@param[in] cq_s - for efficiency
	!>solves the 1-d advection equation:
	!>\f$ \frac{\partial \psi}{\partial t} + \frac{\partial u \psi}{\partial x} = 0 \f$
    subroutine lax_wendroff_ll(ip,jp,o_halo,dt,g,u,v,h,hs,re,&
    		theta,thetan,dtheta,dthetan, phi, phin, dphi, dphin, f_cor, &
    		recqdq, recqdp, recqdp_s, recqdq_s, redq_s, redq, rect, rect_s, cq, cq_s)

		use numerics_type
		implicit none
		integer(i4b), intent(in) :: ip,jp,o_halo
		real(wp), intent(in) :: dt, g, re
		real(wp), intent(in), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: &
																		hs, f_cor, &
    				recqdq, recqdp, recqdp_s, recqdq_s, redq_s, redq, rect, &
    				rect_s, cq, cq_s
		real(wp), dimension(1-o_halo:jp+o_halo), intent(in) :: &
											theta,thetan, dtheta, dthetan
		real(wp), dimension(1-o_halo:ip+o_halo), intent(in) :: &											
											phi, phin, dphi, dphin
		real(wp), intent(inout), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: &
																				h, u, v
																				
		! local variables:
		real(wp), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: & 
					dy1, v1, h1, vh, uh, vh1, Ux, Uy, Vx, Vy, Vy2
		real(wp), dimension(1:ip,1:jp) :: &
										uh_new, vh_new, h_new
										
		real(wp), dimension(0:ip,1:jp) :: h_mid_xt, uh_mid_xt, Ux_mid_xt, Vx_mid_xt
		real(wp), dimension(0:ip,1:jp) :: vh_mid_xt
		real(wp), dimension(1:ip,0:jp) :: h_mid_yt, uh_mid_yt, vh_mid_yt,  &
												Uy_mid_yt, Vy_mid_yt, Vy_mid_yt2
		integer(i4b) :: j, i
			
		
		v1=v*cq ! cq=cos(theta)
		h1=h*cq
		
		uh=u *h
		vh=v*h
		vh1=v1*h
		
		! continuity equation (calculate mid-point values at 0.5*dt):
		h_mid_xt = 0.5_wp*(h(1:ip+1,1:jp)+h(0:ip,1:jp)) &
		  -(0.5_wp*dt/(recqdp_s(0:ip,1:jp))) &
		  *(uh(1:ip+1,1:jp)-uh(0:ip,1:jp))
		  
		h_mid_yt = 0.5_wp*(h(1:ip,1:jp+1)+h(1:ip,0:jp)) &
		  -(0.5_wp*dt/(recqdq_s(1:ip,0:jp))) &
		  *(vh1(1:ip,1:jp+1)-vh1(1:ip,0:jp))

		! v-phi, or u momentum equation (calculate mid-point values at 0.5*dt):
		Ux = uh*u + g*h**2*(0.5_wp)
		Uy = uh*v1
		uh_mid_xt(0:ip,:) = 0.5_wp*(uh(1:ip+1,1:jp)+uh(0:ip,1:jp)) &
		  -(0.5_wp*dt/(recqdp_s(0:ip,1:jp)))* &
		  	(Ux(1:ip+1,1:jp)-Ux(0:ip,1:jp)) &
		  +0.125_wp*dt*(f_cor(1:ip+1,1:jp)+f_cor(0:ip,1:jp))* &
		  	(vh(1:ip+1,1:jp)+vh(0:ip,1:jp))
		
		uh_mid_yt = 0.5_wp*(uh(1:ip,1:jp+1)+uh(1:ip,0:jp)) &
		  -(0.5_wp*dt/(recqdq_s(1:ip,0:jp)))* &
		  	(Uy(1:ip,1:jp+1)-Uy(1:ip,0:jp)) &
		  +0.125_wp*dt*(f_cor(1:ip,1:jp+1)+f_cor(1:ip,0:jp))* &
		  	(vh(1:ip,1:jp+1)+vh(1:ip,0:jp))

		! v-theta, or v momentum equation (calculate mid-point values at 0.5*dt):
		Vx = uh*v
		Vy = vh1*v
		Vy2 = 0.5_wp*g*h**2
		vh_mid_xt(0:ip,1:jp) = 0.5_wp*(vh(1:ip+1,1:jp)+vh(0:ip,1:jp)) &
		  -(0.5_wp*dt/(recqdp_s(0:ip,1:jp)))*(Vx(1:ip+1,1:jp)-Vx(0:ip,1:jp)) &
		  -0.125_wp*dt*(f_cor(1:ip+1,1:jp)+f_cor(0:ip,1:jp))*(uh(1:ip+1,1:jp)+uh(0:ip,1:jp))

		vh_mid_yt(1:ip,0:jp) = 0.5_wp*(vh(1:ip,1:jp+1)+vh(1:ip,0:jp)) &
		  -(0.5_wp*dt/(recqdq_s(1:ip,0:jp)))*(Vy(1:ip,1:jp+1)-Vy(1:ip,0:jp)) &
		  -(0.5_wp*dt/(redq_s(1:ip,0:jp)))*(Vy2(1:ip,1:jp+1)-Vy2(1:ip,0:jp)) &
		  -0.125_wp*dt*(f_cor(1:ip,1:jp+1)+f_cor(1:ip,0:jp))*(uh(1:ip,1:jp+1)+uh(1:ip,0:jp))

		! Now use the mid-point values to predict the values at the next timestep
		! continuity:
		h_new = h(1:ip,1:jp) &
		  - (dt/(recqdp(0:ip-1,1:jp)))*(uh_mid_xt(1:ip,1:jp)-uh_mid_xt(0:ip-1,1:jp)) &
		  - (dt/(recqdq(1:ip,0:jp-1))) * &
		  (vh_mid_yt(1:ip,1:jp)*cq_s(1:ip,1:jp)-vh_mid_yt(1:ip,0:jp-1)*cq_s(1:ip,0:jp-1))

		! u-momentum equation:
		Ux_mid_xt = uh_mid_xt*uh_mid_xt/h_mid_xt + 0.5_wp*g*h_mid_xt**2
		Uy_mid_yt = uh_mid_yt*vh_mid_yt/h_mid_yt*cq_s(1:ip,0:jp)
		uh_new = uh(1:ip,1:jp) &
		  - (dt/(recqdp(0:ip-1,1:jp)))*  (Ux_mid_xt(1:ip,1:jp)-Ux_mid_xt(0:ip-1,1:jp)) &
		  - (dt/(recqdq(1:ip,0:jp-1)))*(Uy_mid_yt(1:ip,1:jp)-Uy_mid_yt(1:ip,0:jp-1))

		! v-momentum equation:
		Vx_mid_xt = uh_mid_xt*vh_mid_xt/h_mid_xt
		Vy_mid_yt = vh_mid_yt*vh_mid_yt/h_mid_yt*cq_s(1:ip,0:jp)
		Vy_mid_yt2 = 0.5_wp*g*h_mid_yt**2
		vh_new = vh(1:ip,1:jp) &
		  - (dt/(recqdp(0:ip-1,1:jp)))*(Vx_mid_xt(1:ip,1:jp)-Vx_mid_xt(0:ip-1,1:jp)) &
		  - (dt/(recqdq(1:ip,0:jp-1)))*(Vy_mid_yt(1:ip,1:jp)-Vy_mid_yt(1:ip,0:jp-1)) &
		  - (dt/(redq(1:ip,0:jp-1) ))* &
		  (Vy_mid_yt2(1:ip,1:jp)-Vy_mid_yt2(1:ip,0:jp-1))

		! add on Coriolis and contribution of orography to pressure gradient:
		uh_new=uh_new  +dt*.5_wp*(f_cor(1:ip,1:jp)*v(1:ip,1:jp) - &
			g*(hs(2:ip+1,1:jp)-hs(0:ip-1,1:jp))/(recqdp(1:ip,1:jp)+recqdp(0:ip-1,1:jp)))* &
			(h(1:ip,1:jp)+h_new)

		vh_new=vh_new  -dt*.5_wp*(f_cor(1:ip,1:jp)*u(1:ip,1:jp) + &
			g*(hs(1:ip,2:jp+1)-hs(1:ip,0:jp-1))/(redq(1:ip,1:jp)+redq(1:ip,0:jp-1)))* &
			(h(1:ip,1:jp)+h_new)

		! re-calculate u and v.
		u(1:ip,1:jp) = uh_new/(h_new)
		v(1:ip,1:jp) = vh_new/h_new
		h(1:ip,1:jp) = h_new

	end subroutine lax_wendroff_ll

	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates del2 of prognostic variable
	!>@param[in] ip: number of east-west points
	!>@param[in] jp: ditto for north-south
	!>@param[in] o_halo: halos required for advection scheme
	!>@param[in] dt:  timestep
	!>@param[in] f: prognostic variable
	!>@param[inout] delsq: delsq of f
	!>@param[in] re: radius of planet
	!>@param[in] theta: latitude
	!>@param[in] thetan: latitude - staggered
	!>@param[in] dtheta: latitude step
	!>@param[in] dthetan: latitude step - staggered
	!>@param[in] phi: phi
	!>@param[in] phin: phin
	!>@param[in] dphi: dphi
	!>@param[in] dphin: dphin
	!>@param[in] recq: for efficiency
	!>@param[in] cq_s: for efficiency
	!>@param[in] dp1: for efficiency
	!>@param[in] dq: for efficiency
	!>calculates del**2:
	!>\f$ visterm = \frac{1}{re^2\cos\theta}
	!>  \frac{\partial }{\partial \theta} 
	!>\left(\cos\theta\frac{\partial f}{\partial\theta} \right) + 
	!> \frac{1}{re^2\cos^2\theta}\frac{\partial^2 f}{\partial \phi ^2}\f$
    subroutine dissipation(ip,jp,o_halo,dt,f,delsq,re,&
    		theta,thetan,dtheta,dthetan, phi, phin, dphi, dphin, &
    		recq, cq_s, dp1, dq)

		use numerics_type
		implicit none
		integer(i4b), intent(in) :: ip,jp,o_halo
		real(wp), intent(in) :: dt, re
		real(wp), dimension(1-o_halo:jp+o_halo), intent(in) :: &
											theta,thetan, dtheta, dthetan
		real(wp), dimension(1-o_halo:ip+o_halo), intent(in) :: &											
											phi, phin, dphi, dphin
		real(wp), intent(in), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: &
																				f, &
															recq, cq_s, dp1, dq
		real(wp), intent(inout), dimension(1:ip,1:jp) :: delsq
		! local variables:
		integer(i4b) :: j, i
			
		! calculate del^2 using 2nd order difference 
		! (central difference of forward and backward):
		delsq(1:ip,1:jp)  =1._wp/(re*recq(1:ip,1:jp))* &
			( cq_s(1:ip,1:jp)*(f(1:ip,2:jp+1)-f(1:ip,1:jp))/dq(1:ip,1:jp) - &
			  cq_s(1:ip,0:jp-1)*(f(1:ip,1:jp)-f(1:ip,0:jp-1))/dq(1:ip,0:jp-1) ) / &
			  dq(1:ip,1:jp)
			  
		delsq(1:ip,1:jp)  = delsq(1:ip,1:jp) + &
			1._wp/(recq(1:ip,1:jp)**2._wp)* &
			( (f(2:ip+1,1:jp)-f(1:ip,1:jp))/dp1(1:ip,1:jp) - &
			  (f(1:ip,1:jp)-f(0:ip-1,1:jp))/dp1(1:ip,0:jp-1) ) / &
			  dp1(1:ip,1:jp)

	end subroutine dissipation

	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates smagorinsky-lilly viscosity
	!>@param[in] ip: number of east-west points
	!>@param[in] jp: ditto for north-south
	!>@param[in] o_halo: halos required for advection scheme
	!>@param[in] cvis:  coefficient for viscosity
	!>@param[in] u,v: u and v winds
	!>@param[inout] vis: viscosity
	!>@param[in] re: radius of planet
	!>@param[in] recq: for efficiency
	!>@param[in] dp1: for efficiency
	!>@param[in] dq: for efficiency
	!>calculates smagorinsky-lilly viscosity:
	!>\f$ visco = C_s^2\Delta x\Delta y|S|\f$
    subroutine smagorinsky(ip,jp,o_halo,cvis,u,v,vis,re,&
    		recq, dp1, dq)

		use numerics_type
		implicit none
		integer(i4b), intent(in) :: ip,jp,o_halo
		real(wp), intent(in) :: cvis, re
		real(wp), intent(in), dimension(1-o_halo:ip+o_halo,1-o_halo:jp+o_halo) :: &
																		u,v, &
															recq, dp1, dq
		real(wp), intent(inout), dimension(1:ip,1:jp) :: vis
		! local variables:
		integer(i4b) :: j, i	
		
		! calculate viscosity using centred differences:
		vis(1:ip,1:jp) = cvis**2._wp*re*recq(1:ip,1:jp)*dp1(1:ip,1:jp)*dq(1:ip,1:jp)* &
			sqrt( ( (u(2:ip+1,1:jp)-u(0:ip-1,1:jp))/ &
			(recq(1:ip,1:jp)*(dp1(0:ip-1,1:jp)+dp1(1:ip,1:jp))) )**2._wp + &
			( (v(1:ip,2:jp+1)-v(1:ip,0:jp-1))/ &
			(re*(dq(1:ip,1:jp)+dq(1:ip,0:jp-1))) )**2._wp + &
			0.5_wp*(  &
			(u(1:ip,2:jp+1)-u(1:ip,0:jp-1))/ &
			(re*(dq(1:ip,1:jp)+dq(1:ip,0:jp-1)))+ &
			(v(2:ip+1,1:jp)-v(0:ip-1,1:jp))/ &
			(recq(1:ip,1:jp)*(dp1(1:ip,1:jp)+dp1(0:ip-1,1:jp))) &
			)**2._wp )

	end subroutine smagorinsky

	end module advection