	program secTermFirst
	
		implicit none


		double precision, dimension(0:26,0:26) :: matrix, EVEC
		
		real, dimension(0:26,0:26)	:: realMatrix, realEVEC
	
		double precision, dimension(0:2) :: params

		real, dimension (0:26) :: realEVAL

		double precision, dimension (0:26) :: EVAL, tmp, anotherTMP
	
		common /coeff/ matrix, realMatrix
		
		logical :: readFiles
	
		integer :: i,j,k

		double precision :: DEPISF
	
1		format(27F7.1)

2		format(27F10.5)

3		format(27F15.10)

4		format(1F10.5)

5		format(1F15.10)
		
		params = (/1.0000000 , 0.0156250 , 0.0000000 /)
	
		open (unit = 1, file = "Z:\\res\\secTermFirst.txt")

		!open (unit = 1, file = "X:\\Computational Mathematics\\Results\\23413.1\\Zhiltsov\\thirdTask.txt")
		
		write(1, *), "Zhiltsov Nikita 23413/1"
		
		do i=0,2

			print*,readFiles(params(i))
			
			write(1, *), "param (p)= "

			write(1, 3), params(i)
			
			write(1, *), "current matrix (A):"
	
			do j=0,26

				write(1, 1),(matrix(j, k),k=0,26)
	
			end do

			CALL EVCSF (27, realMatrix, 27, realEVAL, realEVEC, 27)

			CALL DEVCSF (27, Matrix, 27, EVAL, EVEC, 27)

			write(1, *), "eigenvalues and corresponding eigenVectors :"

			do j =0,26

				write(1, *), "single:"

				write(1, 4), realEVAL(j)

				write(1, 2), (realEVEC(k, j), k=0,26)

				write(1, *), "double"

				write(1, 5), EVAL(j)

				write(1, 3), (EVEC(k, j), k=0,26)

				write(1, *), "residual :"

				do k = 0,26

					tmp(k) = realEVEC(k, j)

				end do

				anotherTMP = tmp

				CALL DSYMV ("U"  , 27, -1.0d0, matrix, 27 , tmp, 1   , EVAL(j), anotherTMP, 1)

				write(1, 3), anotherTMP

				write(1, *), "*"
			
			end do

			write(1, *), "performance index:"

			write(1, *), DEPISF(27, 27, matrix, 27, EVAL, EVEC, 27)

			CALL MXTXF (27, 27, realEVEC, 27, 27, realMatrix, 27)

			write(1, *), "orthogonal check:"
	
			do j=0,26

				write(1, 2),(realMatrix(j, k),k=0,26)
	
			end do

			write(1, *), "_______________________________________________"

			write(1, *), ""

		end do
		
		close(1)
	
	end program secTermFirst
	
	!****************************************************************************	
	
	function readFiles(param) result(done)
	
		double precision, intent(in) :: param
		
		double precision, dimension(0:26,0:26) :: matrix
		
		real, dimension(0:26,0:26) :: realMatrix	
	
	    double precision, dimension(0:26) :: vec

		common /coeff/ matrix, realMatrix
		
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

		done = .TRUE.
	
	end function readFiles
