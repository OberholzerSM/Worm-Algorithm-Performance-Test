!v2: Rewrite the entire program.
module Worm7Vertex_v2
      
    implicit none
    public :: generateWormData
    private
    
    integer :: N0,N1                        !Lattice dimensions.
    integer :: i,j,k						!Loop indices.
    real	:: M,p0							!Mass and proposal probability.
    integer, allocatable :: lattice(:,:,:)	!Lattice of vertices with dimensions N0xN1x4
    
    contains

      !Create n_data many closed loop configurations. Will be called from the main C++ project.
    subroutine generateWormData(N0C,N1C,MC,n_data) bind(C,name="generate7VertexWormDataFortran_v2") 
    
        use iso_c_binding, only: c_int, c_float                 !C Types.
        !DEC$ ATTRIBUTES DLLEXPORT :: generateWormData
        integer(kind=c_int), intent(in) :: N0C,N1C,n_data       !Input Variables from C of type int.
        real(kind=c_float),  intent(in) :: MC                   !Input Variables from C of type double.
        
        !Convert from C-Type to Fortran-Type
        N0 = INT(N0C)
        N1 = INT(N1C)
        M  = REAL(MC)
        !p0 = 1.0 / REAL(N0*N1)
        p0 = 0.5
        
        !Allocate lattice.
        allocate(lattice(N0,N1,4))
        lattice(:,:,:) = 0
        
        !Begin Worm
		do i = 2,n_data
			call Worm    
        end do
        
        !Deallocate lattice.
        deallocate(lattice)
        
    end subroutine generateWormData

    subroutine Worm
    
		logical :: check_end, check_end_prop										!check_end: Check if the main loop should end. check_end_prop: To propose the main loop to end.
		integer :: direction														!The direction the worm will move towards next.
        integer, dimension(2) :: start,x,y											!Start location, current wormhead location, next proposed wormhead location.
        integer, dimension(4) :: oldvertex_x, oldvertex_y, newvertex_x, newvertex_y	!Vertices.
        real :: q1,q2,pA,u															!Proposal and Acceptance probabilities.
        real :: weightx, weighty, weightx_new, weighty_new							!Vertex-Weights.
        
        !Suggest a new starting location for the worm.
        start(1) = randomint(1,N0)
        start(2) = randomint(1,N1)
        
        !Determine the acceptance probability.
        oldvertex_x = lattice(start(1),start(2),:)
		weightx = weight(oldvertex_x,.false.)
		weightx_new = weight(oldvertex_x,.true.)
		q1 = 1.0 / REAL(N0*N1)
        q2 = p0
		pA = (weightx_new/weightx) * (q2/q1)
        
        !Test if the new start will be accepted.
        call random_number(u)
        if(u < pA) then
            check_end = .false.
            x = start
        else
            check_end = .true.
        end if
        
        !If the start is accepted, begin the main worm loop.
        WormLoop: do while(check_end == .false.)
            
            !Set the current vertex at location x.
            oldvertex_x = lattice(x(1),x(2),:)
            
            !If a loop has been completed, propose to remove the worm.
            check_end_prop = .false.
            if(x(1) == start(1) .and. x(2) == start(2)) then
				call random_number(u)
                if(u < p0) then
					check_end_prop = .true.
                    y(:) = 0
                    newvertex_x = oldvertex_x
                    oldvertex_y = oldvertex_x
                    newvertex_y = oldvertex_x
                end if
            end if
            
            !Otherwise, propose a new location for the worm to move to.
            if(check_end_prop == .false.) then
                
                !Propose a new direction to move towards.
                direction = randomint(1,4)
                
                select case(direction)
                  
                !If moving right.
                case(1)
                    !Define the new location y.
                    if(x(1) < N0) then
                        y(1) = x(1) + 1
                    else
                        y(1) = 1    !Boundry Condition
                    end if
                    y(2) = x(2)
                    oldvertex_y = lattice(y(1),y(2),:)
                    
                    !Define the new vertices.
                    !If no bond exists, create a new one.
                    if (oldvertex_x(direction)==0) then
				        newvertex_x = oldvertex_x + (/ 1,0,0,0 /)
				        newvertex_y = oldvertex_y + (/ 0,0,1,0 /)
			        !If a bond already exists, delete it.
                    else 
				        newvertex_x = oldvertex_x + (/ -1,0,0,0 /)
				        newvertex_y = oldvertex_y + (/ 0,0,-1,0 /)
                    end if
                    
                !If moving up.
                case(2)
                    !Define the new location y.
                    if(x(2) < N1) then
                        y(2) = x(2) + 1
                    else
                        y(2) = 1    !Boundry Condition
                    end if
                    y(1) = x(1)
                    oldvertex_y = lattice(y(1),y(2),:)
                    
                    !Define the new vertices.
                    !If no bond exists, create a new one.
                    if (oldvertex_x(direction)==0) then
				        newvertex_x = oldvertex_x + (/ 0,1,0,0 /)
				        newvertex_y = oldvertex_y + (/ 0,0,0,1 /)
			        !If a bond already exists, delete it.
                    else 
				        newvertex_x = oldvertex_x + (/ 0,-1,0,0 /)
				        newvertex_y = oldvertex_y + (/ 0,0,0,-1 /)
                    end if
                
                !If moving left.
                case(3)
                    !Define the new location y.
                    if(x(1) > 1) then
                        y(1) = x(1) - 1
                    else
                        y(1) = N0    !Boundry Condition
                    end if
                    y(2) = x(2)
                    oldvertex_y = lattice(y(1),y(2),:)
                    
                    !Define the new vertices.
                    !If no bond exists, create a new one.
                    if (oldvertex_x(direction)==0) then
				        newvertex_x = oldvertex_x + (/ 0,0,1,0 /)
				        newvertex_y = oldvertex_y + (/ 1,0,0,0 /)
			        !If a bond already exists, delete it.
                    else 
				        newvertex_x = oldvertex_x + (/ 0,0,-1,0 /)
				        newvertex_y = oldvertex_y + (/ -1,0,0,0 /)
                    end if
                        
                !If moving down.
                case(4)
                    !Define the new location y.
                    if(x(2) > 1) then
                        y(2) = x(2) - 1
                    else
                        y(2) = N1    !Boundry Condition
                    end if
                    y(1) = x(1)
                    oldvertex_y = lattice(y(1),y(2),:)
                    
                    !Define the new vertices.
                    !If no bond exists, create a new one.
                    if (oldvertex_x(direction)==0) then
				        newvertex_x = oldvertex_x + (/ 0,0,0,1 /)
				        newvertex_y = oldvertex_y + (/ 0,1,0,0 /)
			        !If a bond already exists, delete it.
                    else 
				        newvertex_x = oldvertex_x + (/ 0,0,0,-1 /)
				        newvertex_y = oldvertex_y + (/ 0,-1,0,0 /)
                    end if
                    
                end select     
            end if
            
            !!Determine the vertex-weights.
            !Determine the vertex-weights if one removes the worm.
            if(check_end_prop == .true.) then
                weightx = weight(oldvertex_x,.true.)
                weightx_new = weight(newvertex_x,.false.)
                weighty = 1.0
                weighty_new = 1.0
            else 
                !Determine the vertex-weights at x if one is at the start-location.
                if (x(1) == start(1) .and. x(2) == start(2)) then
                    weightx = weight(oldvertex_x,.true.)
                    !Do not accept leaving behind a 3-vertex at the start!
                    if( sum(newvertex_x) == 3 ) then
                        weightx_new = 0.0
                    else
                        weightx_new = weight(newvertex_x,.true.) !.true., since the source remains at the start.
                    end if
                !No special conditions: Worm away from x.
                else
                    weightx = weight(oldvertex_x,.true.)
                    weightx_new = weight(newvertex_x,.false.)
                end if
            
                !Determine the vertex-weights at y if one moves to the start-location.
                if (y(1) == start(1) .and. y(2) == start(2)) then
                    weighty = weight(oldvertex_y,.true.)        !.true., since there is already a source at the start.
                    weighty_new = weight(newvertex_y,.true.)
                !No special conditions: Worm moves to y.
                else
                    weighty = weight(oldvertex_y,.false.)
                    weighty_new = weight(newvertex_y,.true.)
                end if
                
            end if
            
            !Determine the proposal and accpetance probability
            q1 = proposal(x,y,oldvertex_x)
            q2 = proposal(y,x,newvertex_y)
            pA = (weightx_new/weightx) * (weighty_new/weighty) * (q2/q1)
 
            !Test if the step will be accepted.
            call random_number(u)
            if(u < pA) then
                !If one removes the worm.
                if(check_end_prop == .true.) then
                    check_end = .true.
                    exit Wormloop
                !If one continues, apply the changes.
                else
                    lattice(x(1),x(2),:) = newvertex_x
                    lattice(y(1),y(2),:) = newvertex_y
                    x = y
                end if
            end if
            
        end do WormLoop
        
    end subroutine Worm
    
	function weight(vertex,sink) !Weight of a Vertex. sink: Source/Sink at the Vertex.

		integer, dimension(4), intent(in) :: vertex
        logical, intent(in) :: sink
        logical :: straightLine
		integer :: nBonds
		real :: weight
		
		nBonds = sum(vertex)
	
		if (sink) then !If there is a sink (or source) at the Vertex.
	
			select case (nBonds)
		
				case(0) !Sink and Source at the Vertex.
				weight = 4.0
			
				case(1)
				weight = 1.0
			
                case(2)
                straightLine = (vertex(1)==1 .and. vertex(3)==1) .or. (vertex(2)==1 .and. vertex(4)==1)
				if ( straightLine ) then !Straight horizontal or vertical line.
					weight = 1.0
				else !Corner
					weight = 0.5
				end if
			
				case(3)
				weight = 0.25**(1.0/3.0)
			
				case default
				weight = 0.0
		
			end select
	
		else !If there is no sink (or source) at the Vertex.
	
			select case (nBonds)
		
				case(0)
				weight=M**2
			
				case(2)
				straightLine = (vertex(1)==1 .and. vertex(3)==1) .or. (vertex(2)==1 .and. vertex(4)==1)
				if ( straightLine ) then !Straight horizontal or vertical line.
					weight = 1.0
				else !Corner
					weight = 0.5
				end if
			
				case default
				weight = 0.0
		
			end select
	
		end if

	end function weight
    
    function proposal(x,y,oldvertex_x) !Returns the odds of proposing the step x -> y.

		integer, dimension(2), intent(in) :: x,y
		integer, dimension(4), intent(in) :: oldvertex_x
		integer :: n_bonds
		real :: proposal
	
		if ( x(1) == 0 ) then !If a loop has been completed.
			proposal = 1.0 / real(N0*N1)
		else
			!Determine the number of bonds at the worm-head.
			n_bonds = sum(oldvertex_x)

			select case(n_bonds)
		
				case(0)!Valid Vertex
				if ( sum(y) == 0 ) then
					proposal = p0
				else
					proposal = (1.0-p0) / 4.0
				end if
		
				case(1) !1-point Vertex, all directions possible.
				proposal = 1.0 / 4.0
		
				case(2)!Valid Vertex
				if ( sum(y) == 0 ) then
					proposal = p0
				else
					proposal = (1.0-p0) / 4.0
				end if
		
				case(3) !3-point Vertex, all directions proposed.
				proposal = 1.0 / 4.0
		
				case default
				proposal = 0.0
		
			end select
		end if

	end function proposal
    
    function randomint(a,b) !Random Integer between and including a and b.

        integer, intent(in) :: a,b
	    integer :: randomint
	    real :: u
	
	    call random_number(u)
	    randomint = a + floor(real(b-a+1)*u)
	
    end function randomint
    
end module Worm7Vertex_v2