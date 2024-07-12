	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>mpi routines for shallow water model
    module mpi_module
    use numerics_type
    use mpi
    implicit none
    
#if VAR_TYPE==0
		integer(i4b), parameter :: MPIREAL=MPI_REAL4
#endif
#if VAR_TYPE==1
		integer(i4b), parameter :: MPIREAL=MPI_REAL8
#endif
    private
    public :: mpi_define, block_ring, exchange_halos
    
	contains
	
	! DEFINE SOME TYPES
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>define some types to be used in the model
	!>@param[inout] MPI_INTEGER9_DEF - type to be defined as integer kind=9
	subroutine mpi_define(MPI_INTEGER9_DEF)
		implicit none
		integer(i4b), intent(inout) :: MPI_INTEGER9_DEF
		
		integer(i4b) :: error
		
		call MPI_TYPE_CREATE_F90_INTEGER (9, MPI_INTEGER9_DEF, error)
	end subroutine mpi_define
	

	! EXCHANGE HALOS FOR A VARIABLE USING CARTESIAN TOPOLOGY
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>define some types to be used in the model
	!>@param[in] comm2, id, ipp, jpp, o_halo
	!>@param[inout] array: the array to exchange_halos on
	subroutine exchange_halos(comm2d, id, ipp, jpp, o_halo, array)
		implicit none
		integer(i4b), intent(in) :: comm2d, id, ipp, jpp, o_halo
		real(wp), intent(inout), &
			 dimension(1-o_halo:o_halo+ipp,1-o_halo:o_halo+jpp) :: array
		
		integer(i4b) :: error, nbrleft, nbrright, nbrbottom, nbrtop, tag1, &
						request
		integer(i4b), dimension(MPI_STATUS_SIZE) :: status
		
		! Find the processors neighbours
		call MPI_CART_SHIFT( comm2d, 0, 1, nbrleft, nbrright, error)	
		call MPI_CART_SHIFT( comm2d, 1, 1, nbrbottom, nbrtop, error)
		
		! print *,id,nbrleft, nbrright, nbrbottom, nbrtop

		
		! now receive data from top, bottom, left, and right
		if (nbrleft /= id) then
			tag1=010
			! send from left (specify destination):
			call MPI_Issend(array(ipp+1-o_halo:ipp,1:jpp), jpp, MPIREAL, nbrright, &
				tag1, MPI_COMM_WORLD, request,error)
			! receive from left (specify source):
			call MPI_Recv(array(1-o_halo:0,1:jpp), jpp, MPIREAL, nbrleft, &
				tag1, MPI_COMM_WORLD, status,error)
			call MPI_Wait(request, status, error)
		else
			array(1-o_halo:0,1:jpp)=array(ipp+1-o_halo:ipp,1:jpp)
		endif
		
		
		if (nbrright /= id) then
			tag1=010
			! send from right (specify destination):
			call MPI_Issend(array(1:o_halo,1:jpp), jpp, MPIREAL, nbrleft, &
				tag1, MPI_COMM_WORLD, request,error)
			! receive from right (specify source):
			call MPI_Recv(array(ipp+1:ipp+o_halo,1:jpp), jpp, MPIREAL, nbrright, &
				tag1, MPI_COMM_WORLD, status,error)
			call MPI_Wait(request, status, error)
		else
			array(ipp+1:ipp+o_halo,1:jpp)=array(1:o_halo,1:jpp)
		endif
		
		
		if ((nbrtop /= id)) then
			tag1=110
			! send from top (specify destination):
			call MPI_Issend(array(1:ipp,jpp+1-o_halo:jpp), ipp, MPIREAL, nbrtop, &
				tag1, MPI_COMM_WORLD, request,error)
			! receive from top (specify source):
			call MPI_Recv(array(1:ipp,1-o_halo:0), ipp, MPIREAL, nbrbottom, &
				tag1, MPI_COMM_WORLD, status,error)
			call MPI_Wait(request, status, error)
		else
			! array(1:ipp,1-o_halo:0)=array(1:ipp,jpp-o_halo:jpp)
		endif


		if ((nbrbottom /= id)) then
			tag1=110
			! send from bottom (specify destination):
			call MPI_Issend(array(1:ipp,1:o_halo), ipp, MPIREAL, nbrbottom, &
				tag1, MPI_COMM_WORLD, request,error)
			! receive from bottom (specify source):
			call MPI_Recv(array(1:ipp,jpp+1:jpp+o_halo), ipp, MPIREAL, nbrtop, &
				tag1, MPI_COMM_WORLD, status,error)
			call MPI_Wait(request, status, error)
		else
			! array(1:ipp,jpp+1:jpp+o_halo)=array(1:ipp,1:o_halo)
		endif	
		
					
	end subroutine exchange_halos
	
	
	! BLOCK VIA RING
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>block by sending / receiving a message around the ring
	!>@param[in] ring_comm - comm of the cart topology
	!>@param[in] id - id of this process
	!>@param[in] world_process - id of world process
	!>@param[in] rank - rank of mpi job
	subroutine block_ring(ring_comm,id,world_process,rank)
		implicit none
		integer(i4b), intent(in) :: ring_comm, id, world_process, rank
		integer(i4b) :: error, tag1=2010
		character (len=3) :: mesg='Yo!'

        ! ESSENTIALLY BLOCKS UNTIL ALL PROCESSORS CATCH UP
        call MPI_Barrier(ring_comm, error)
	end subroutine block_ring
	
	
	end module mpi_module
	
	