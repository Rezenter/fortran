	program third
	
		implicit none


		double precision, dimension(0:26,0:26) :: matrix, facMat
		
		real, dimension(0:26,0:26)	:: realMatrix, realFacMat
	
		double precision, dimension(0:26) :: vec, res, residual

		real, dimension(0:26) :: realVec, realRes, realResidual
	
		double precision, dimension(0:2) :: params

		double precision :: rcond, r, b

		real :: realR, realB
	
		common /coeff/ matrix, realMatrix, vec, realVec, b, realB
		
		logical :: readFiles
	
		integer :: i,j,k
	
1		format(27F7.1)

2		format(27F10.3)

3		format(27F15.6)
		
		params = (/1.0000000 , 0.0156250 , 0.0000000 /)
	
		open (unit = 1, file = "Z:\\res\\thirdTask.txt")

		!open (unit = 1, file = "X:\\Computational Mathematics\\Results\\23413.1\\Zhiltsov\\thirdTask.txt")
		
		write(1,*), "Zhiltsov Nikita 23413/1"
		
		do i=0,2

			b = 0;
		
			print*,readFiles(params(i))
			
			write(1,*), "param (p)= "

			write(1,3), params(i)
			
			write(1,*), "current matrix (A):"
	
			do j=0,26

				write(1, 1),(matrix(j, k),k=0,26)
	
			end do
	
			write(1,*),"current vector (b): "

			write(1,1), vec

	
			!CALL LFCRG (N, A, LDA, FAC, LDFAC, IPVT, RCOND)
			!Compute the LU factorization of a real general matrix and estimate its L1=condition number.
			!CALL LFCRG (27, realMatrix, 27, realFacMat, 27, realIPVT, rcond)
			!same for symm positive-defined
			!CALL LFCDS (N, A, LDA, FAC, LDFAC, RCOND)
			CALL LFCDS (27, realMatrix, 27, realFacMat, 27, rcond)

			!CALL LFSRG (N, FAC, LDFAC, IPVT, B, IPATH, X)
			!IPATH — Path indicator. (Input)
			!IPATH = 1 means the system AX = B is solved.
			!IPATH = 2 means the system A7X = B
			!CALL LFSRG (27, realFacMat, 27, realIPVT, realVec, 1, realRes)
			!same for symm positive-defined 
			CALL LFSDS (27, realFacMat, 27, realVec, realRes)

			write(1,*),"result vector (real) : "
			
			write(1,2), realRes
			
			!up
			!CALL DLFCRG (27, matrix, 27, facMat, 27, IPVT, rcond)
			CALL DLFCDS (27, matrix, 27, facMat, 27, rcond)

			!up
			!CALL DLFSRG (27, facMat, 27, IPVT, vec, 1, res)
			CALL DLFSDS (27, facMat, 27, vec, res)

			write(1,*),"result vector (double) : "
			
			write(1,3), res
			
			write(1,*), "condition:"

			write(1,*), 1/rcond

			write(1,*), "reverse condition:"

			write(1,*), rcond

			realResidual = realVec

			!Matrix–Vector Multiply (TRANS, M, N, SALPHA, SA, LDA, SX, INCX, SBETA, SY, INCY)
			!CALL SGEMV ("N", 27, 27, -1.0, realMatrix, 27, realRes, 1, 1.0, realResidual, 1)
			CALL SSYMV ("U", 27, -1.0, realMatrix, 27, realRes, 1, 1.0, realResidual, 1)
	
			write(1,*), "residual (r):"

			write(1,*), realResidual

			write(1,*), "error ceiling:"

			realR = 0;

			do j=0,26
				
				realR = realR + ABS(realResidual(j))

			end do

			write(1,*), realR/(realB*rcond)
			
			residual = vec

			!Matrix–Vector Multiply (TRANS, M, N, DALPHA, DA, LDA, DX, INCX, DBETA, DY, INCY)
			!DA should be given as double
			!CALL DGEMV ("N", 27, 27, -1.0d0, matrix, 27, res, 1, 1.0d0, residual, 1)
			CALL DSYMV ("U", 27, -1.0d0, matrix, 27, res, 1, 1.0d0, residual, 1)
	
			write(1,*), "residual (r) (double):"

			write(1,*), residual

			write(1,*), "error ceiling (double):"

			r = 0;

			do j=0,26
				
				r = r + ABS(residual(j))

			end do

			write(1,*), r/(b*rcond)
				
			write(1,*), "_______________________________________________"

			write(1,*), ""

		end do
		
		close(1)
	
	end program third
	
	!****************************************************************************	
	
	function readFiles(param) result(done)
	
		double precision, intent(in) :: param
		
		double precision, dimension(0:26,0:26) :: matrix
		
		real, dimension(0:26,0:26) :: realMatrix	
	
	    double precision, dimension(0:26) :: vec

		real, dimension(0:26) :: realVec

		double precision :: b

		real :: realR
	
		common /coeff/ matrix, realMatrix, vec, realVec, b, realB
		
		logical :: done
	
		integer :: i,j
	
		open (unit = 2, file = "X:\\Computational Mathematics\\Tasks\\23413.1\\Zhiltsov\\RPOmtr067.dat")
	
		do i=0,26
	
			read(2,*) vec
	
			do j=0,26
	
				matrix(i, j) = vec(j)

				if (i==0 .AND. J == 0)  then
					
					matrix(i, j) = matrix(i, j) + param
	
				end if

				realMatrix(i, j) = real(matrix(i, j))

			end do
	
		end do
	
		close(2)
	
		open (unit = 2, file = "X:\\Computational Mathematics\\Tasks\\23413.1\\Zhiltsov\\RPOrhs067.dat")
	
		do i=0,26
	
			read(2,*) vec(i)
	
			if (i == 0) then
	
				vec(i) = vec(i) + 8*param
	
			end if

			realVec(i) = real(vec(i))
	
		end do
	
		close(2)
		
		do j=0,26
				
			b = b + ABS(vec(j))
			realB = realB + ABS(realVec(j));

		end do
	
		done = .TRUE.
	
	end function readFiles
