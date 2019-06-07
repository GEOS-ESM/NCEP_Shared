subroutine gfsio_solve_axb ( A, lda, m, b, ipvt, j0 )
!$$$ document block
!
! module:   gfsio_solve_axb     API for general linear system solve 
!
! Abstract: Uses either ESSL or LAPACK to solve Ax=b
!
! Program history log
!    2009-01-07    Todling
!
!$$$ end document block

implicit none
integer(4) lda,m,j0
integer(4) ipvt(lda)
real(8) a(lda,m),b(lda)
integer(4) info
#ifdef ibm_sp
        CALL DGEF(A,lda,m,IPVT)
        CALL DGES(A,lda,m,IPVT,b,J0)
#else /* ibm_sp */
        CALL DGETRF(m,m,A,lda,IPVT,INFO)
          if(INFO/=0)then
            print*,'trouble factorizing matrix, info= ', INFO
            call exit(99)
          endif
        CALL DGETRS('N',m,1,A,lda,IPVT,b,m,INFO)
          if(INFO/=0)then
            print*,'trouble solving system of eqs, info= ', INFO
            call exit(99)
          endif
#endif /* ibm_sp */
return
end subroutine gfsio_solve_axb
