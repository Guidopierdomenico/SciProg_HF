module LinearSystem

   implicit none
   private
   public :: linearsolver

   interface linearsolver
      module procedure linearsolver
   end interface

contains

   subroutine linearsolver (C, a, epsilon)
      !inout and output variables of the subroutine
      real(8), intent(in)     :: C(:,:), a(:)
      real(8), intent(out)    :: epsilon(:)
      !local variables needed for DGESV
      real(8)                 :: matrix(size(C,1),size(C,2)), b_vector(size(a))
      integer                 :: n, nrhs, lda, ldb, info
      integer                 :: ipiv(size(C,1))


      ! Error checking
      if (size(C,1) /= size(C,2)) then
         print*," linearsolver assumes square matrices"
         stop "error in linearsolver routine"
      endif

      if (size(C,1) /= size(a)) then
         print*," number of rows of the matrix does not match with dimension of the vector "
         stop "error in linearsolver routine"
      endif

      if (size(C,1) > size(epsilon)) then
         print*," dimension of result array too small "
         stop "error in linearsolver routine"
      endif
      
      !arguments that DGESV takes are n(num equations), nrhs (num rhs vectors), a(matrix), lda(leading dimension of matrix), ipiv(integer array of size n, used to track row swaps) , b(b vector), ldb(leading dimension b vector), info(scalar integer that outputs 0 if successful)
      n = size(C,1)
      lda = n
      ldb = size(a)
      nrhs = 1
      !creating copies of the input matrix and vector, since the DGESV call needs to overwrite them during the calculation
      matrix = C
      b_vector = a ! vector to copy the input vector, which will be overwritten with the result by DGESV

      !calling the DGESV routine
      call DGESV(n, nrhs, matrix, lda, ipiv, b_vector, ldb, info)

      ! if the DGESV runs without error the info variable will be zero and in that case we can assign the result to the epsilon otherwise we print an error message
      if (info == 0) then
         epsilon = b_vector
      else
         print *, "error in linearsolver routine"
         stop "error in linearsolver routine"
      endif

   end subroutine

end module
