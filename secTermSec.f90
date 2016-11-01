	program secTermSec
	
		implicit none

		INTEGER, PARAMETER::LDA=14,N=14,LDEVEC=14

		double precision :: matrix(0:LDA-1, 0:N-1)

		double complex :: doubleEVEC(0:LDEVEC-1, 0:N-1)

		real :: realMatrix(0:LDA-1, 0:N-1)
		
		complex	:: EVEC(0:LDEVEC-1,0:N-1), complexMatrix(0:LDA-1, 0:N-1)
	
		double precision, dimension(0:2) :: params

		complex EVAL(0:N-1), tmp(0:N-1), anotherTMP(0:N-1)

		double complex :: doubleEVAL(0:N-1)
	
		common /coeff/ matrix, realMatrix
	
		logical :: readFiles
	
		integer :: i,j,k

		double precision :: EPIRG
	
1		format(14F7.1)

2		format(14F10.5)

3		format(14F15.10)

4		format(1F10.5)

5		format(1F15.10)
		
		params = (/1.0000000 , 0.0156250 , 0.0000000 /)
	
		open (unit = 1, file = "Z:\\res\\secTermSec.txt")


		!open (unit = 1, file = "X:\\Computational Mathematics\\Results\\23413.1\\Zhiltsov\\thirdTask.txt")
		
		write(1, *), "Zhiltsov Nikita 23413/1"
		
		do i=0,2

			print*,readFiles(params(i))
			
			write(1, '(A,F12.10)'), "param (p)= ", params(i)
			
			write(1, *), "current matrix (A):"
	
			do j=0,N-1

				write(1, 1),(matrix(j, k),k=0,N-1)
	
			end do

			CALL EVCRG (N, realMatrix, LDA, EVAL, EVEC, LDEVEC)

			CALL DEVCRG (N, matrix, LDA, doubleEVAL, doubleEVEC, LDEVEC)

			write(1, *), "eigenvalues and corresponding eigenVectors :"

			complexMatrix = realMatrix

			do j =0,N-1

				write(1, *), "single:"
				
				write (1,'(A,F10.5,A,F10.5,A)'),"(", real(EVAL(j)), ' + i ' , aimag(EVAL(j)), ")"

				write(1,'(A,14(F10.5,A,F10.5,A))'),"(", (real(EVEC(k, j)), ' + i ' , aimag(EVEC(k, j)), '); (', k=0,13)

				write(1, *), "double:"

				write(1, '(A, F15.10,A,F15.10, A)'), "(", real(doubleEVAL(j)), ' + i ' , aimag(doubleEVAL(j)), ")"

				write(1,'(A,14(F15.10,A,F15.10,A))'), "(", (real(doubleEVEC(k, j)), ' + i ' , aimag(doubleEVEC(k, j)), "); (", k=0,13)
				write(1, *), "residual :"

				do k = 0,N-1

					tmp(k) = EVEC(k, j)


				end do
				anotherTMP = tmp

				!CALL CGEMV (TRANS, M, N, CALPHA, CA, LDA, CX, INCX, CBETA, CY, INCY)
				!For all data types, A is an NxN matrix. These subprograms set y <- aAx + beta*y

				CALL CGEMV ("N", N, N, (-1.0,0.0), complexMatrix, LDA, tmp, 1, EVAL(j), anotherTMP, 1)

				write(1, '(A, 14(F15.10,A,F15.10, A2))'), "(", (real(anotherTMP(k)), ' + i ' , aimag(anotherTMP(k)), "); (", k=0,13)

				write(1, *), "*"
			
			end do

			write(1, *), "performance index:"

			write(1, *), EPIRG(LDA, N, matrix, N, EVAL, EVEC, LDEVEC)

			write(1, *), "_______________________________________________"

			write(1, *), ""

		end do
		
		close(1)
	
	end program secTermSec
	
	!****************************************************************************	
	
	function readFiles(param) result(done)
	
		double precision, intent(in) :: param
		
		integer, parameter :: N = 14, LDA = 14

		double precision :: matrix(0:LDA-1, 0:N-1)

		real :: realMatrix(0:LDA-1, 0:N-1)
		
	    double precision :: vec(0:N-1)

		common /coeff/ matrix, realMatrix
		
		logical :: done
	
		integer :: i,j
	
		open (unit = 2, file = "X:\\Computational Mathematics\\Tasks\\23413.1\\Zhiltsov\\RGEmtr067.dat")
	
		do i=0,13
	
			read(2,*) vec
	
			do j=0,13
	
				matrix(i, j) = vec(j)

				if (i==0 .AND. J == 0)  then
					
					matrix(i, j) = matrix(i, j) + param
	
				end if

				realMatrix(i, j) = real(matrix(i, j))

			end do
	
		end do
	
		close(2)

		done = .TRUE.
	
	end function readFiles
