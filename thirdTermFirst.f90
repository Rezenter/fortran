	program thirdTermSecond
	
		implicit none

		double precision :: params(0:8), max, tmp, abs, f, x(0:1), xGuess(0:1), scal

		external f

		integer :: i

		params = (/1.0d-5, 1.0d-6, 1.0d-7, 1.0d-8, 1.0d-9, 1.0d-10, 1.0d-11, 1.0d-12, 1.0d-13/)
	
		open (unit = 1, file = "Z:\\res\\thirdTermSec.txt")
		
		write(1,*), "Zhiltsov Nikita 23413/1"

		write(1, *), "sin(x + 0.5) - y - 1 = 0"

		write(1, *), "x + cos(y - 2) = 0"

		xGuess = (0.538 , -0.139)
		
		do i=0,8
			
			write(1, '(A, F15.13)'), "precision = ", params(i)
			
			CALL DNEQNF (f, params(i), 2, 200, xGuess, x, scal)
			
			write(1, '(A, f17.15, A, f19.17)'), "root = ", x(0), " ; ", x(1)		

			write(1, *), "_______________________________________________"

			write(1, *), ""

		end do
		
		close(1)
	
	end program thirdTermSecond

	
	subroutine f(x, f, n)
	
		integer :: n

		double precision :: f(0:1), x(0:1), dcos, dsin

		f(0) = dsin(x(0) + 0.5) - x(1) - 1

		f(1) = x(0) + dcos(x(1) - 2)

		return

	end