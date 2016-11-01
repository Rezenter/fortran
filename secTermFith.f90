	program secTermSixth
	
		implicit none
	
		integer, dimension(0:4) :: params

		double precision :: h, a, b, datan, dcos, dsin, dble, xArray(0:160), fArray(0:160, 0:160), max(0:2), tmp(0:2), prev(0:2), coeff(0:160, 0:160),  sVal(0:2, 0:8*160, 0:8*160), xArrayEx(0:8*160), knots(0:160 + 3), abs, DBS2DR

		integer :: i, j, k, m

		params = (/10, 20, 40, 80, 160 /)
	
		open (unit = 1, file = "Z:\\res\\secTermSixth.txt")
		
		write(1,*), "Zhiltsov Nikita 23413/1"

		write(1, *), "f(x, y) = x^2 * y * sin(x*y)"

		write(1, *), "a = c = 0, b = d = Pi"
		
		do i=0,4

			write(1, *), "N = ", params(i)
		
			a = 0.0d0

			b = 4.D0*DATAN(1.D0)
	
			h = (b-a)/dble(8*params(i))

			do j = 0, 8*params(i)

				xArrayEx(j) = a + dble(j)*h

			end do

			do j = 0, params(i)

				xArray(j) = a + dble(j*8)*h
			
			end do

			do j = 0, params(i)

				do k = 0, params(i)

					fArray(j, k) = (xArray(j) ** (2)) * xArray(k) * dsin(xArray(j) * xArray(k))
				
				end do
			
			end do

			CALL DBSOPK (params(i) + 1, xArray, 3, knots)
			
			do j = 1, 4

				max(j) = 0

				t(j) = 0

			end do

			CALL DBS2IN (params(i) + 1, xArray, params(i) + 1, xArray, fArray, params(i) + 1, 3, 3, knots, knots, coeff)

			do j = 0, 8*params(i)

				do k = 0, 8*params(i)

					sVal(0, j, k) = DBS2DR(0, 0, xArrayEx(j), xArrayEx(k), 3, 3, knots, knots, params(i) + 1, params(i) + 1, coeff)

					sVal(1, j, k) = DBS2DR(1, 0, xArrayEx(j), xArrayEx(k), 3, 3, knots, knots, params(i) + 1, params(i) + 1, coeff)

					sVal(2, j, k) = DBS2DR(0, 1, xArrayEx(j), xArrayEx(k), 3, 3, knots, knots, params(i) + 1, params(i) + 1, coeff)

					tmp(0) = abs(sVal(0, j, k) - (xArrayEx(j) ** (2)) * xArrayEx(k) * dsin(xArrayEx(j) * xArrayEx(k)))

					tmp(1) = abs(sVal(1, j, k) - (xArrayEx(j) * xArrayEx(k) * (2 * dsin(xArrayEx(j) * xArrayEx(k) + xArrayEx(j) * xArrayEx(k) * dcos(xArrayEx(j) * xArrayEx(k))))))

					tmp(2) = abs(sVal(2, j, k) - ((xArrayEx(j) ** (2)) * (dsin(xArrayEx(j) * xArrayEx(k)) + xArrayEx(j) * xArrayEx(k) * dcos(xArrayEx(j) * xArrayEx(k)))))

					do m = 0, 2

						if(max(m) < tmp(m)) then

							max(m) = tmp(m)
					
						end if
				
					end do

				end do
			
			end do

			write(1, '(A,  F15.12)')," Max error(f) = ", max(0)

			write(1, '(A,  F15.12)')," Max error(x derivative) = ", max(1)

			write(1, '(A,  F15.12)')," Max error(y derivative) = ", max(2)

			if(prev(0) .NE. 0) then

				write(1, '(A, F5.1)'), " relation(f) = ", prev(0)/max(0)

				write(1, '(A, F5.1)'), " relation(x derivative) = ", prev(1)/max(1)

				write(1, '(A, F5.1)'), " relation(y derivative) = ", prev(2)/max(2)

			end if

			do j = 0, 2

				prev(j) = max(j)

			end do
			
			write(1,*), "_______________________________________________"

			write(1,*), ""

		end do
		
		close(1)
	
	end program secTermSixth
