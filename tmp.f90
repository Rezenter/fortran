	program thirdTermThird
	
		implicit none

		double precision :: max, prev, abs, x(0:1024), h, r, l, dble, dexp, solv(0:1, 0:1, 1024), h2, y(0:1, 0:1024), a(0:1, 0:1024), b(0:1, 0:1024), c(0:1, 0:1024), gamma(0:1, 0:1), dsin, datan, dcos, beta, f, Pi

		integer :: i, j, params(0:6), k

		common/coeff/ x

		params = (/16, 32, 64, 128, 256, 512, 1024/)
	
		open (unit = 1, file = "Z:\\res\\thirdTermThird.txt")
		
		write(1,*), "Zhiltsov Nikita 23413/1"

		write(1, *), "equation: (d(e^(-x)*(du/dx)))dx - u*e^(-x) = 1 - x"

		write(1, *), "boundary equation: k*du/dx = 1"

		write(1, *), "boundary equation: k*du/dx = 4"

		write(1, '(A)'), "analytic solution of D5 task: u = C1* e^((1 - 5^0.5)x/2) + C2* e^((1 + 5^0.5)x/2) + x*e^x"

		write(1, *), "test equation : d2y/dx^2 + 4*Pi^2 * y = 0 with boundary equations dy/dx = 2*Pi"

		write(1, *), "x from 0 to 3"

		write(1, *), ""

		r = dble(0)

		l = dble(3)

		Pi = 4.D0*DATAN(1.D0)

		gamma(0, 0) = 1

		gamma(0, 1) = 4

		gamma(1, 0) = 2*Pi

		gamma(1, 1) = 2*Pi
		
		do i=0,6
			
			write(1, '(A, I5.0)'), "interval count = ", params(i)
			
			h = (l - r)/params(i)

			h2 = h**2

			write(1, '(A, F12.10)'), "step = ", h

			write(1, *), "numerical solution of D5 task:"

			do j = 0, params(i)

				x(j) = r + h * j

				a(0, j) = -h/2 - 1

				a(1, j) = -1

				b(0, j) = 2 + h2

				b(1, j) = 2 - h2*4*Pi**2

				c(0, j) = h/2 - 1

				c(1, j) = -1

				if (j == 1) then

					do k = 0, 1

						solv(k, 0, j) = (a(k,  j)/3 - c(k,  j))/(b(k,  j) + 4*a(k,  j)/3)

						solv(k, 1, j) = -h2*(f(k,  j) - (2*a(k,  j)*gamma(k, 0))/(3*beta(k, 0, j)*h))/(b(k,  j) + 4*a(k,  j)/3)

					end do
				
				else

					if (j > 1 .AND. j < params(i) - 1) then

						do k = 0, 1

							solv(k, 0, j) = (- c(k,  j))/(b(k,  j) + a(k,  j)*solv(k, 0, j - 1))

							solv(k, 1, j) = -(a(k,  j)*solv(k, 1, j - 1) + h2*f(k,  j))/(b(k,  j) + a(k,  j)*solv(k, 0, j - 1))

						end do
				
					end if

				end if
	
			end do

			do j = params(i) - 1, 1, -1

				if(j == params(i) - 1) then

					do k = 0, 1

						y(k,  j) = -((a(k,  j) - c(k,  j)/3)*solv(k, 1, j - 1) + h2*(f(k,  j) + (2*gamma(k, 1)*c(k,  j))/(3*beta(k, 1, j)*h)))/(b(k,  j) + 4*c(k,  j)/3 + solv(k, 0, j - 1)*(a(k,  j) - c(k,  j)/3))

					end do

				else

					do k = 0, 1
					
						y(k,  j) = solv(k, 0, j)*y(k, j+1) + solv(k, 1, j)

					end do

				end if

			end do

			do j = 0, 1

				y(j, 0) = (4*y(j, 1) - y(j, 2) - (2*h*gamma(j, 0)/ beta(j, 0, 0)))/3

				y(j, params(i)) = ((2*h*gamma(j, 1))/beta(j, 1, 0) + 4*y(j, params(i) - 1) - y(j, params(i) - 2))/3

			end do

			do j = 0, params(i)

				write(1, '(A, F10.7, A, F12.10)'), "y = ", y(0,  j), " at x = ", x(j)

				if(max < abs(y(1,  j) - sin(2*Pi*x(j)))) then

					max = abs(y(1,  j) - sin(2*Pi*x(j)))

				end if

			end do

			write(1, '(A, F7.5)'), "residual of test equations solution = ", max

			if(prev > 0) then

				write(1, '(A, F3.1)'), "residual relation = ", prev/max

			end if

			prev = max

			max = 0
		
			write(1, *), "_______________________________________________"

			write(1, *), ""

		end do

		write(1, *), "end of file."
		
		close(1)
	
	end program thirdTermThird

	

	function f(i, n) result(res)

		integer, intent(in) :: i, n

		double precision res, dexp, x(0:1024)

		common/coeff/ x

		if(i == 0) then
			
			res = (1 - x(n))*dexp(x(n))

		else

			res = 0

		end if

	end function f



	function beta(i, j, n) result(res)

		integer, intent(in) :: i, j, n

		double precision :: x(0:1024), res
			
		common/coeff/ x

		if(i == 0) then
			
			res = dexp(-x(n))

		else

			res = 1

		end if

	end function beta

