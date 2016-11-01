	program first
	
		implicit none


		double precision, dimension(0:13,0:13) :: matrix, facMat
		
		real, dimension(0:13,0:13)	:: realMatrix, realFacMat
	
		double precision, dimension(0:13) :: vec, res, residual, IPVT

		real, dimension(0:13) :: realVec, realRes, realResidual, realIPVT
	
		double precision, dimension(0:2) :: params

		double precision :: rcond, r, b

		real :: realR, realB
	
		common /coeff/ matrix, realMatrix, vec, realVec, b, realB
		
		logical :: readFiles
	
		integer :: i,j,k
	
1		format(14F7.1)

2		format(14F10.3)

3		format(14F15.6)
		
		params = (/1.0000000 , 0.0156250 , 0.0000000 /)
	
		open (unit = 1, file = "Z:\\res\\firstTask.txt")

		!open (unit = 1, file = "X:\\Computational Mathematics\\Results\\23413.1\\Zhiltsov\\firstTask.txt")
		
		write(1,*), "Zhiltsov Nikita 23413/1"
		
		do i=0,2

			b = 0;
		
			print*,readFiles(params(i))
			
			write(1,*), "param (p)= "

			write(1,3), params(i)
			
			write(1,*), "current matrix (A):"
	
			do j=0,13

				write(1, 1),(matrix(j, k),k=0,13)
	
			end do
	
			write(1,*),"current vector (b): "

			write(1,1), vec

	
			!CALL LFCRG (N, A, LDA, FAC, LDFAC, IPVT, RCOND)
			!Compute the LU factorization of a real general matrix and estimate its L1=condition number.
			CALL LFCRG (14, realMatrix, 14, realFacMat, 14, realIPVT, rcond)

			!CALL LFSRG (N, FAC, LDFAC, IPVT, B, IPATH, X)
			!IPATH — Path indicator. (Input)
			!IPATH = 1 means the system AX = B is solved.
			!IPATH = 2 means the system A7X = B
			CALL LFSRG (14, realFacMat, 14, realIPVT, realVec, 1, realRes)

			write(1,*),"result vector (real) : "
			
			write(1,2), realRes
		
			!up
			CALL DLFCRG (14, matrix, 14, facMat, 14, IPVT, rcond)

			!up
			CALL DLFSRG (14, facMat, 14, IPVT, vec, 1, res)

			write(1,*),"result vector (double) : "
			
			write(1,3), res
			
			write(1,*), "condition:"

			write(1,*), 1/rcond

			write(1,*), "reverse condition:"

			write(1,*), rcond

			realResidual = realVec

			!Matrix–Vector Multiply (TRANS, M, N, SALPHA, SA, LDA, SX, INCX, SBETA, SY, INCY)
			CALL SGEMV ("N", 14, 14, -1.0, realMatrix, 14, realRes, 1, 1.0, realResidual, 1)
	
			write(1,*), "residual (r):"

			write(1,*), realResidual

			write(1,*), "error ceiling:"

			realR = 0;

			do j=0,13
				
				realR = realR + ABS(realResidual(j))

			end do

			write(1,*), realR/(realB*rcond)
			
			residual = vec

			!Matrix–Vector Multiply (TRANS, M, N, DALPHA, DA, LDA, DX, INCX, DBETA, DY, INCY)
			!DA should be given as double
			CALL DGEMV ("N", 14, 14, -1.0d0, matrix, 14, res, 1, 1.0d0, residual, 1)
	
			write(1,*), "residual (r) (double):"

			write(1,*), residual

			write(1,*), "error ceiling (double):"

			r = 0;

			do j=0,13
				
				r = r + ABS(residual(j))

			end do

			write(1,*), r/(b*rcond)
				
			write(1,*), "_______________________________________________"

			write(1,*), ""

		end do
		
		close(1)
	
	end program first
	
	!****************************************************************************	
	
	function readFiles(param) result(done)
	
		double precision, intent(in) :: param
		
		double precision, dimension(0:13,0:13) :: matrix
		
		real, dimension(0:13,0:13) :: realMatrix	
	
	    double precision, dimension(0:13) :: vec

		real, dimension(0:13) :: realVec

		double precision :: b

		real :: realR
	
		common /coeff/ matrix, realMatrix, vec, realVec, b, realB
		
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
	
		open (unit = 2, file = "X:\\Computational Mathematics\\Tasks\\23413.1\\Zhiltsov\\RGErhs067.dat")
	
		do i=0,13
	
			read(2,*) vec(i)
	
			if (i == 0) then
	
				vec(i) = vec(i) + 5*param
	
			end if

			realVec(i) = real(vec(i))
	
		end do
	
		close(2)
		
		do j=0,13
				
			b = b + ABS(vec(j))
			realB = realB + ABS(realVec(j));

		end do
	
		done = .TRUE.
	
	end function readFiles
