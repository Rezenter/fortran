	program fourth
	
		implicit none
		

		double precision, dimension(0:4,0:24) :: matrix, facMat
		
		real, dimension(0:4,0:24)	:: realMatrix, realFacMat
	
		double precision, dimension(0:24) :: vec, res, residual, IPVT

		real, dimension(0:24) :: realVec, realRes, realResidual, realIPVT
	
		double precision, dimension(0:2) :: params

		double precision :: rcond, r, b

		real :: realR, realB

		integer :: i,j,k, NLCA, NUCA
	
		common /coeff/ matrix, realMatrix, vec, realVec, b, realB, NLCA, NUCA
		
		logical :: readFiles
	
1		format(25F7.1)

2		format(25F10.3)

3		format(25F15.6)
		
		params = (/1.0000000 , 0.0156250 , 0.0000000 /)

		NLCA = 3

		NUCA = 1
	
		open (unit = 1, file = "Z:\\res\\fourthTask.txt")

		!open (unit = 1, file = "X:\\Computational Mathematics\\Results\\23413.1\\Zhiltsov\\thirdTask.txt")
		
		write(1,*), "Zhiltsov Nikita 23413/1"
		
		do i=0,2

			b = 0;
		
			print*,readFiles(params(i))
			
			write(1,*), "param (p)= "

			write(1,3), params(i)
			
			write(1,*), "current matrix (A):"
	
			do j=0,24

				do k=0,24

					if (k < j + NUCA + 1 .AND. k > j - NLCA - 1) then

						if(k == j) then

							res(k) = matrix(NUCA, k)
							
						else 
							if(k > j) then

								res(k) = matrix(k - j - NUCA, 2*k - j - NUCA)

							else

								res(k) = matrix(j - k + NUCA, k)

							end if

						end if
						
					else

					res(k) = 0

					end if

				end do
			
				write(1, 3),(res(k),k=0,24)

			end do
	
			write(1,*),"current vector (b): "

			write(1,3), vec

	
			!CALL LFCRG (N, A, LDA, FAC, LDFAC, IPVT, RCOND)
			!Compute the LU factorization of a real general matrix and estimate its L1=condition number.
			CALL LFCRB (5, realMatrix, 25, NLCA, NUCA, realFacMat, 25, IPVT, rcond)


			!CALL LFSRG (N, FAC, LDFAC, IPVT, B, IPATH, X)
			!IPATH — Path indicator. (Input)
			!IPATH = 1 means the system AX = B is solved.
			!IPATH = 2 means the system A7X = B
			CALL LFSRB (5, realFacMat, 25, NLCA, NUCA, realIPVT, realVec, 1, realRes)
			
			write(1,*),"result vector (real) : "
			
			write(1,2), realRes
			
			
			CALL DLFCRB (5, matrix, 25, NLCA, NUCA, facMat, 25, IPVT, rcond)

			CALL DLFSRB (5, facMat, 25, NLCA, NUCA, IPVT, vec, 1, res)
			
			write(1,*),"result vector (double) : "
			
			write(1,3), res
			
			write(1,*), "condition:"

			write(1,*), 1/rcond

			write(1,*), "reverse condition:"

			write(1,*), rcond

			realResidual = realVec

			!Matrix–Vector Multiply (TRANS, M, N, SALPHA, SA, LDA, SX, INCX, SBETA, SY, INCY)
			!(TRANS, M, N, NLCA, NUCA SALPHA, SA, LDA, SX,            INCX, SBETA,SY, INCY)
			CALL SGBMV ("N", 25, 25, NLCA, NUCA, -1.0, realMatrix, 25, realRes, 1, 1.0, realResidual, 1)
	
			write(1,*), "residual (r):"

			write(1,*), realResidual

			write(1,*), "error ceiling:"

			realR = 0;

			do j=0,24
				
				realR = realR + ABS(realResidual(j))

			end do

			write(1,*), realR/(realB*rcond)
			
			residual = vec

			!Matrix–Vector Multiply (TRANS, M, N, DALPHA, DA, LDA, DX, INCX, DBETA, DY, INCY)
			!DA should be given as double
			CALL DGBMV ("N", 25, 25, NLCA, NUCA, -1.0d0, matrix, 25, res, 1, 1.0d0, residual, 1)
		
			write(1,*), "residual (r) (double):"

			write(1,*), residual

			write(1,*), "error ceiling (double):"

			r = 0;

			do j=0,24
				
				r = r + ABS(residual(j))

			end do

			write(1,*), r/(b*rcond)
				
			write(1,*), "_______________________________________________"

			write(1,*), ""

		end do
		
		close(1)
	
	end program fourth
	
	!****************************************************************************	
	
	function readFiles(param) result(done)
	
		double precision, intent(in) :: param
		
		double precision, dimension(0:4,0:24) :: matrix
		
		real, dimension(0:4,0:24) :: realMatrix	
	
	    double precision, dimension(0:24) :: vec

		real, dimension(0:24) :: realVec

		double precision :: b

		real :: realR

		integer :: i,j, NLCA, NUCA
	
		common /coeff/ matrix, realMatrix, vec, realVec, b, realB, NLCA, NUCA
		
		logical :: done
	
		open (unit = 2, file = "X:\\Computational Mathematics\\Tasks\\23413.1\\Zhiltsov\\RGBmtr067.dat")
	
		do i=0, (NLCA + NUCA)
	
			read(2,*) vec
	
			do j=0,24
	
				matrix(i, j) = vec(j)

				if (i == NUCA .AND. j == 0)  then
					
					matrix(i, j) = matrix(i, j) + param
	
				end if

				realMatrix(i, j) = real(matrix(i, j))

			end do
	
		end do
	
		close(2)
	
		open (unit = 2, file = "X:\\Computational Mathematics\\Tasks\\23413.1\\Zhiltsov\\RGBrhs067.dat")
	
		do i=0,24
	
			read(2,*) vec(i)
	
			if (i == 0) then
	
				vec(i) = vec(i) + 2*param
	
			end if

			realVec(i) = real(vec(i))
	
		end do
	
		close(2)
		
		do j=0,24
				
			b = b + ABS(vec(j))
			realB = realB + ABS(realVec(j));

		end do
	
		done = .TRUE.
	
	end function readFiles
