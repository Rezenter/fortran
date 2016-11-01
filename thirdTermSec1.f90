	program thirdTermSecond
	
		implicit none

		double precision :: params(0:6), max, tmp(0:1), abs, x(0:1), xGuess(0:1), scal

		external func, jac

		integer :: i, count(0:1)

		common /coeff/ count

		params = (/1.0d-5, 1.0d-6, 1.0d-7, 1.0d-8, 1.0d-9, 1.0d-10, 1.0d-11/)
	
		open (unit = 1, file = "Z:\\res\\thirdTermSec.txt")
		
		write(1,*), "Zhiltsov Nikita 23413/1"

		write(1, *), "sin(x + 0.5) - y - 1 = 0"

		write(1, *), "x + cos(y - 2) = 0"

		xGuess = (0.538 , -0.139)
		
		do i=0,6

			count(0) = 0

			count(1) = 0
			
			write(1, '(A, F15.13)'), "precision = ", params(i)
			
			CALL DNEQNF (func, params(i), 2, 200, xGuess, x, scal)

			write(1, *), "Auto Jacobian"
			
			write(1, '(A, f18.15, A, f20.17, A, I3.1)'), "root = ", x(0), " ; ", x(1), " iterations = ", count(0)	

			call func(x, tmp, 0)

			write(1, '(A, f20.18, A, f20.18)'), "residual = ", abs(tmp(0)), " ; ", abs(tmp(1))
			
			count(0) = 0

			CALL DNEQNJ (func, jac, params(i), 2, 200, xGuess, x, scal)

			write(1, *), "manual Jacobian"
			
			write(1, '(A, f18.15, A, f20.17, A, I3.1, A, I3.1)'), "root = ", x(0), " ; ", x(1), " function iterations = ", count(0), " jacobian iterations = ", count(1)

			call func(x, tmp, 0)

			write(1, '(A, f20.18, A, f20.18)'), "residual = ", abs(tmp(0)), " ; ", abs(tmp(1))

			write(1, *), "_______________________________________________"

			write(1, *), ""

		end do
		
		close(1)
	
	end program thirdTermSecond

	subroutine func(x, f, n)
	
		integer :: n, count(0:1)

		double precision :: f(0:1), x(0:1), dcos, dsin

		common /coeff/ count

		count(0) = count(0) + 1

		f(0) = dsin(x(0) + 0.5) - x(1) - 1

		f(1) = x(0) + dcos(x(1) - 2)

		return

	end

	subroutine jac(n, x, f)
	
		integer :: n, count(0:1)

		double precision :: f(0:1, 0:1), x(0:1), dcos, dsin

		common /coeff/ count

		count(1) = count(1) + 1

		f(0, 0) = dcos(x(0) + 0.5)

		f(0, 1) = -1

		f(1, 0) = 1

		f(1, 1) = -dsin(x(1) - 2)

		return

	end