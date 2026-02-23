!v5: Same as v4, but the weight and proposal functions were improved.
module Worm7Vertex_v5
    
    implicit none
    public :: generateWormData
    private
    integer :: i !Loop Variable
    
    contains
    
    !Create n_data many closed loop configurations. Will be called from the main C++ project.
    subroutine generateWormData(N0,N1,M,n_data) bind(C,name="generate7VertexWormDataFortran_v5") 
    
        use iso_c_binding, only: c_int, c_float                 !C Types.
        !DEC$ ATTRIBUTES DLLEXPORT :: generateWormData
        integer(kind=c_int), intent(in) :: N0,N1,n_data         !Input Variables from C of type int.
        real(kind=c_float),  intent(in) :: M                    !Input Variables from C of type double.
        integer,dimension(N0,N1,4) :: lattice
        
        !Set lattice to zero.
        lattice(:,:,:) = 0
        
        !Create n_data more closed loop configurations.
        do i=2,n_data
            call Worm(INT(N0),INT(N1),REAL(M),lattice)
        end do
        
    end subroutine generateWormData

	subroutine Worm(N0,N1,M,lattice)

		implicit none
        integer, intent(in) :: N0,N1
        integer,dimension(N0,N1,4), intent(inout) :: lattice
		integer :: check_end,check_end_prop
		integer :: x0,x1,y0,y1,start0,start1,direction
		integer, dimension(4) :: oldvertex_x,oldvertex_y,newvertex_x,newvertex_y
		real :: M,p0,q1,q2,pA,u
		real :: weightx,weighty,weightx_new,weighty_new
	
		!Schlage einen Startort für den neuen Wurm vor
		start0 = randomint(1,N0)
		start1= randomint(1,N1)
	
		!Wahrscheinlichkeit, eine Quelle/Senke wegzunehmen
		!p0 = 1.0 / REAL(N0*N1)
        p0 = 0.5
	
		!Bestimme, ob man die Startkonfiguration akzeptiert
		check_end = 1
		q1 = 1.0/real(N0*N1)
		q2 = p0
		weightx = weight(M,lattice(start0,start1,:),0)
		weightx_new = weight(M,lattice(start0,start1,:),1)
		q1 = 1.0 / REAL(N0*N1)
        q2 = p0
		pA = (weightx_new/weightx) * (q2/q1)
	
		call random_number(u)
		if (pA>u) then !Falls die Startkonfiguration akzeptiert wird
			check_end=0
			x0=start0
			x1=start1
		end if
	
		!Wurmloop: Führe den Wurm weiter, bis man zu einer validen "wurmlosen" Kofniguration zurückwechselt
		Wurmloop: do while(check_end==0)
		
			!Falls man einen Loop gemacht hat, schlage mit p0 Wahrscheinlichkeit vor zurück in die wurmlose config. zu wechseln
			oldvertex_x=lattice(x0,x1,:)
			check_end_prop=0
			if (x0==start0 .and. x1==start1) then
				call random_number(u)
				if (p0>u) then
					check_end_prop=1
					y0=0
					y1=0
					newvertex_x = oldvertex_x
					oldvertex_y = oldvertex_x
					newvertex_y = oldvertex_x
				end if
			end if
		
			!Falls man nicht vorschlägt aufzuhören, schlage vor zum Punkt y weiter zu gehen
			if (check_end_prop==0) then
		
				!Schlage eine neue zufällige Richtung vor
				direction = randomint(1,4)
			
				!Definiere den neuen potentiellen Wurmkopf y sowie die neuen Vertex bei x und y
				select case(direction)
				
				case(1)
				!Definiere den neuen Punkt (y0,y1)
				if (x0 < N0) then
					y0 = x0 + 1
				else !Randbedingung
					y0 = 1
				end if
				y1 = x1
				oldvertex_y = lattice(y0,y1,:)
				!Falls kein Hop in der Richtung existiert, erstelle einen neuen Hop
				if (oldvertex_x(direction)==0) then
					newvertex_x = oldvertex_x + (/ 1,0,0,0 /)
					newvertex_y = oldvertex_y + (/ 0,0,1,0 /)
				else !Falls bereits ein Hop in der Richtung existiert, lösche den Hop
					newvertex_x = oldvertex_x + (/ -1,0,0,0 /)
					newvertex_y = oldvertex_y + (/ 0,0,-1,0 /)
				end if
					
				case(2)
				!Definiere den neuen Punkt (y0,y1)
				y0 = x0
				if (x1 < N1) then
					y1 = x1 + 1
				else !Randbedingung
					y1 = 1
				end if
				oldvertex_y = lattice(y0,y1,:)
				!Falls kein Hop in der Richtung existiert, erstelle einen neuen Hop
				if (oldvertex_x(direction)==0) then
					newvertex_x = oldvertex_x + (/ 0,1,0,0 /)
					newvertex_y = oldvertex_y + (/ 0,0,0,1 /)
				else !Falls bereits ein Hop in der Richtung existiert, lösche den Hop
					newvertex_x = oldvertex_x + (/ 0,-1,0,0 /)
					newvertex_y = oldvertex_y + (/ 0,0,0,-1 /)
				end if
				
				case(3)
				!Definiere den neuen Punkt (y0,y1)
				if (x0 > 1) then
					y0 = x0 - 1
				else !Randbedingung
					y0 = N0
				end if
				y1 = x1
				oldvertex_y = lattice(y0,y1,:)
				!Falls kein Hop in der Richtung existiert, erstelle einen neuen Hop
				if (oldvertex_x(direction)==0) then
					newvertex_x = oldvertex_x + (/ 0,0,1,0 /)
					newvertex_y = oldvertex_y + (/ 1,0,0,0 /)
				else !Falls bereits ein Hop in der Richtung existiert, lösche den Hop
					newvertex_x = oldvertex_x + (/ 0,0,-1,0 /)
					newvertex_y = oldvertex_y + (/ -1,0,0,0 /)
				end if

				case(4)
				!Definiere den neuen Punkt (y0,y1)
				y0 = x0
				if (x1 > 1) then
					y1 = x1 - 1
				else !Randbedingung
					y1 = N1
				end if
				oldvertex_y = lattice(y0,y1,:)
				!Falls kein Hop in der Richtung existiert, erstelle einen neuen Hop
				if (oldvertex_x(direction)==0) then
					newvertex_x = oldvertex_x + (/ 0,0,0,1 /)
					newvertex_y = oldvertex_y + (/ 0,1,0,0 /)
				else !Falls bereits ein Hop in der Richtung existiert, lösche den Hop
					newvertex_x = oldvertex_x + (/ 0,0,0,-1 /)
					newvertex_y = oldvertex_y + (/ 0,-1,0,0 /)
				end if
			
				end select
			end if
		
			!Bestimme die Gewichte bei x
			if (y0==0 .and. y1==0) then !Falls man aufhört, nimm den Wurmkopf weg
				weightx = weight(M,oldvertex_x,1)
				weightx_new = weight(M,newvertex_x,0)
			else if (x0==start0 .and. x1==start1) then !Falls man am Start ist, bleibt der Wurmkopf
				weightx = weight(M,oldvertex_x,1)
				weightx_new = weight(M,newvertex_x,1)
				if ( sum(newvertex_x) ==3) then !Hinterlasse keinen 3-er Vertex beim Start
					weightx_new = 0.0
				end if
			else !Sonst bewegt sich der Wurmkopf weiter zu y
				weightx = weight(M,oldvertex_x,1)
				weightx_new = weight(M,newvertex_x,0)
			end if
		
			!Bestimme die Gewichte bei y
			if (y0==0 .and. y1==0) then !Falls man aufhört, ist y=0.
				weighty = 1.0
				weighty_new = 1.0
			else if (y0==start0 .and. y1==start1) then !Falls man zum Start geht, bleibt der Wurmkopf
				weighty = weight(M,oldvertex_y,1)
				weighty_new = weight(M,newvertex_y,1)
			else !Sonst bewegt sich der Wurmkopf weiter zu y
				weighty = weight(M,oldvertex_y,0)
				weighty_new = weight(M,newvertex_y,1)
			end if
		
			!Bestimme die Proposal-Wahrscheinlichkeiten, von x->y zu gehen, sowie die Akzeptanzwahrscheinlichkeit
			q1 = proposal(N0,N1,x0,x1,y0,y1,p0,oldvertex_x)
			q2 = proposal(N0,N1,y0,y1,x0,x1,p0,newvertex_y)
			pA = (weightx_new/weightx) * (weighty_new/weighty) * (q2/q1)
		
			!Bestimme, ob der Schritt akzeptiert wird
			call random_number(u)
			if ( pA > u ) then
			
				if (y0==0 .and. y1==0) then !Falls man aufhört
					check_end=1
					exit Wurmloop
				else !Falls man mit dem Wurm weitermacht
		
					lattice(x0,x1,:) = newvertex_x
					lattice(y0,y1,:) = newvertex_y
					x0 = y0
					x1 = y1
	
				end if
			end if
		
		end do Wurmloop

	end subroutine Worm

	function weight(M,mu,check_x)!Gewicht von einem Vertex. check_x=1: Quelle/Senke am Vertex.

		implicit none
		integer, intent(in) :: check_x
        integer :: summe
		real, intent(in) :: M
        real :: weight
		integer, dimension(4), intent(in) :: mu
	
		summe = sum(mu)
	
		if (check_x==1) then !Falls man eine Quelle/Senke am Vertex hat
	
			select case (summe)
		
				case(0) !Quelle und Senke am gleichen Ort
				weight = 4.0 !w=4, da man verbundene und unverbundene zusammenzählt
			
				case(1)
				weight = 1.0
			
				case(2)
				if (mu(1)==1 .and. mu(3)==1) then !Gerade horizontale Linie
					weight = 1.0
				else if(mu(2)==1 .and. mu(4)==1) then !Gerade vertikale Linie
					weight = 1.0
				else !Ecke
					weight = 0.5
				end if
			
				case(3)
				weight = 0.25**(1.0/3.0)
			
				case default
				weight = 0.0
		
			end select
	
		else !Falls man keine Quellen/Senken hat
	
			select case (summe)
		
				case(0)
				weight=M**2
			
				case(2)
				if (mu(1)==1 .and. mu(3)==1) then !Gerade horizontale Linie
					weight = 1.0
				else if(mu(2)==1 .and. mu(4)==1) then !Gerade vertikale Linie
					weight = 1.0
				else !Ecke
					weight = 0.5
				end if
			
				case default
				weight = 0.0
		
			end select
	
		end if

	end function weight

	function proposal(N0,N1,x0,x1,y0,y1,p0,oldvertex_x) !Gibt die proposal Wahrscheinlichkeit, von oldvertex_x -> newvertex_x zu gehen.

		implicit none
		integer, intent(in) :: N0,N1,x0,x1,y0,y1
        integer :: n_bonds
		integer, dimension(4), intent(in) :: oldvertex_x
		real, intent(in) :: p0
        real :: proposal
	
		if ( x0 == 0 .and. x1 == 0 ) then
			proposal = 1.0 / real(N0*N1)
		else
			!Bestimme die Anzahl Bonds beim Wurmkopf
			n_bonds = sum(oldvertex_x)
	
			select case(n_bonds)
		
			case(0)!Valider Vertex
			if ( y0 == 0 .and. y1 == 0 ) then
				proposal = p0
			else
				proposal = (1.0-p0) / 4.0
			end if
		
			case(1) !1-point Vertex, alle Richtungen möglich
			proposal = 1.0 / 4.0
		
		
			case(2)!Valider Vertex
			if ( y0 == 0 .and. y1 == 0 ) then
				proposal = p0
			else
				proposal = (1.0-p0) / 4.0
			end if
		
			case(3) !3-point Vertex, nur drei Richtungen möglich (aber alle werden vorgeschlagen)
			proposal = 1.0 / 4.0
		
			case default
			proposal = 0.0
		
			end select
		end if

	end function proposal

	function randomint(a,b) !Zufälliges Integer zwischen a und b (inklusive a und b)

		implicit none
		integer :: randomint,a,b
		real :: u
	
		call random_number(u)
		randomint = a + floor(real(b-a+1)*u)
	
	end function randomint
    
end module Worm7Vertex_v5