	program secTermFourth
	
		implicit none
	
		integer, dimension(0:4) :: params

		double precision :: h, a, b, datan, dcos, dsin, x, dble, xArray(0:160), yArray(0:160), s(0:160), max(0:4), tmp, prev(0:4), coeff(0:3, 0:180),  sVal(0:3,0:4*160), xArrayEx(0:4*160), DCSDER, DCSITG, abs

		integer :: i, j, k

		params = (/10, 20, 40, 80, 160 /)
	
		open (unit = 1, file = "Z:\\res\\secTermFourth.txt")
		
		write(1,*), "Zhiltsov Nikita 23413/1"

		write(1, *), "f(x) = x * sin(x)"

		write(1, *), "a = 0, b = Pi"
		
		do i=0,4
		
			a = 0.0d0

			b = 4.D0*DATAN(1.D0)
	
			h = (b-a)/dble(4*params(i))

			do j = 0, 4*params(i)

				xArrayEx(j) = a + dble(j)*h

			end do

			do j = 0, params(i)

				xArray(j) = a + dble(j*4)*h

				yArray(j) = xArray(j) * dsin(xArray(j))
			
			end do

			CALL DCSINT (params(i) + 1, xArray, yArray, s, coeff)
		
			do j = 0, 4*params(i)

				do k = 0, 3

					max(k) = 0

					sVal(k, j) = DCSDER(k, xArrayEx(j), params(i) + 1, s, coeff)
			
				end do
		
			end do
			
			write(1,'(A, I3.0)'), "param (p)= ", params(i)
			
			do j = 0, 4*params(i)

				tmp = abs(sVal(0, j) - (xArrayEx(j)*dsin(xArrayEx(j))))

				if (tmp > max(0)) then

					max(0) = tmp

				end if

				tmp = sVal(1, j) - (dsin(xArrayEx(j)) + xArrayEx(j)*dcos(xArrayEx(j)))

				if (tmp > max(1)) then

					max(1) = tmp
		
				end if

				tmp = sVal(2, j) - (2.0d0 * dcos(xArrayEx(j)) - xArrayEx(j)*dsin(xArrayEx(j)))

				if (tmp > max(2)) then

					max(2) = tmp
		
				end if

				tmp = sVal(3, j) - (-3.0d0 * dsin(xArrayEx(j)) - xArrayEx(j)*dcos(xArrayEx(j)))

				if (tmp > max(3)) then

					max(3) = tmp
		
				end if
			
			end do
			
			tmp = DCSITG(a, b, params(i) + 1, s, coeff)

			write(1, '(A, F25.23)'), "Integrate S from a to b = ", tmp

			max(4) = tmp - b

			write(1, '(A, F15.12)'), "Integration error = ", max(4)
			
			write(1,'(A, F10.8 )'),"Max Error f = " , max(0)

			write(1,'(A, F8.6 )'),"Max Error firs derivative = " , max(1)

			write(1,'(A, F8.6 )'),"Max Error second derivative = " , max(2)

			write(1,'(A, F8.6 )'),"Max Error third derivative = " , max(3)

			do j = 1,3
				
				max(j) = 0

			end do

			do j = 1, 4*params(i) - 1

				tmp = sVal(1, j) - (dsin(xArrayEx(j)) + xArrayEx(j)*dcos(xArrayEx(j)))

				if (tmp > max(1)) then

					max(1) = tmp
		
				end if

				tmp = sVal(2, j) - (2.0d0 * dcos(xArrayEx(j)) - xArrayEx(j)*dsin(xArrayEx(j)))

				if (tmp > max(2)) then

					max(2) = tmp
		
				end if

				tmp = sVal(3, j) - (-3.0d0 * dsin(xArrayEx(j)) - xArrayEx(j)*dcos(xArrayEx(j)))

				if (tmp > max(3)) then

					max(3) = tmp
		
				end if
			
			end do

			write(1,'(A, F8.6 )'),"Max Error firs derivative(edge points excluded) = " , max(1)

			write(1,'(A, F6.4 )'),"Max Error second derivative(edge points excluded) = " , max(2)

			write(1,'(A, F6.4 )'),"Max Error third derivative(edge points excluded) = " , max(3)

			if (prev(0) .NE. 0.0d0) then
			
				write(1, '(A,F5.2)'),"Relation f = ", prev(0)/max(0)

				write(1, '(A,F5.2)'),"Relation firs derivative(edge points excluded) = ", prev(1)/max(1)

				write(1, '(A,F5.2)'),"Relation second derivative(edge points excluded) = ", prev(2)/max(2)

				write(1, '(A,F5.2)'),"Relation third derivative(edge points excluded) = ", prev(3)/max(3)

				write(1, '(A,F5.2)'),"Relation integral = ", prev(4)/max(4)

			end if
			
			do j = 0, 4
					
				prev(j) = max(j)

			end do
				
			write(1,*), "_______________________________________________"

			write(1,*), ""

		end do
		
		close(1)
	
	end program secTermFourth
