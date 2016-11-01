	program secTermThird
	
		implicit none
	
		integer, dimension(0:4) :: params

		double precision :: h , a, b , datan, dsin, x, dble, xArray(0:160), yArray(0:160), s(0:4*160), xArrayEx(0:4*160), max, tmp, prev

		logical :: compute

		integer :: i, j
	
		common /coeff/ h, a , b, xArray, yArray, s, xArrayEx

		params = (/10, 20, 40, 80, 160 /)

		a = 0.0d0

		b = 4.D0*DATAN(1.D0)
	
		open (unit = 1, file = "Z:\\res\\secTermThird.txt")
		
		write(1,*), "Zhiltsov Nikita 23413/1"

		write(1, *), "f(x) = x * sin(x)"

		write(1, *), "a = 0, b = Pi"

		prev = 0
		
		do i=0,4
		
			print*,compute(params(i))
			
			write(1,'(A, I3.0)'), "param (p)= ", params(i)

			max = 0

			do j = 0, 4*params(i)

				tmp = s(j)-(xArrayEx(j)*dsin(xArrayEx(j)))

				if (tmp > max) then

					max = tmp

				end if
			
			end do
			
			write(1,'(A, F14.12 )'),"Max Error = " , max

			if (prev .NE. 0.0d0) then

				write(1, '(A,F5.2)'),"Relation = " ,prev/max
			
			end if

			prev = max
				
			write(1,*), "_______________________________________________"

			write(1,*), ""

		end do
		
		close(1)
	
	end program secTermThird
	
	!****************************************************************************	
	
	function compute(param) result(done)
	
		integer, intent(in) :: param

		double precision h, a, b, xArray(0:160), yArray(0:160), s(0:4*160), xArrayEx(0:4*160)
		
		common /coeff/ h, a, b, xArray, yArray, s, xArrayEx
		
		logical :: done

		integer :: i
	
		h = (b-a)/dble(4*param)

		do i = 0,4*param

			xArrayEx(i) = a + dble(i)*h

		end do

		do i = 0,param

			xArray(i) = a + dble(i*4)*h

			yArray(i) = xArray(i) * dsin(xArray(i))
			
		end do

		CALL DCSIEZ (param + 1, xArray, yArray, 4*param + 1, xArrayEx, s)

		done = .TRUE.
	
	end function compute
