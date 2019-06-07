	subroutine imax_verify_(imax)
	  implicit none
	  integer,intent(in) :: imax
	  logical :: invalid
	  integer :: m,n
	  m=imax
	  call factor_get_(m,2,n)
	    if(n<1 .or. n>25) call exitn_(2,n,imax)
	  call factor_get_(m,3,n)
	    if(n>2) call exitn_(3,n,imax)
	  call factor_get_(m,5,n)
	    if(n>1) call exitn_(5,n,imax)
	  call factor_get_(m,7,n)
	    if(n>1) call exitn_(7,n,imax)
	  call factor_get_(m,11,n)
	    if(n>1) call exitn_(11,n,imax)
	    if(m>1) call exitm_(m,imax)
	end subroutine imax_verify_
	subroutine exitn_(mp,n,im)
	  implicit none
	  integer,intent(in) :: mp,n,im
	  integer*4 :: ierr=2

	  write(*,'(a,3(a,i5))') 'spffte(): -- ERROR -- ',
     &      'invalid imax =',im,', power of factor ',mp, ' is ',n
          call exit(ierr)
	end subroutine exitn_
	subroutine exitm_(m,im)
	  implicit none
	  integer,intent(in) :: m,im
	  integer*4 :: ierr=2
	  write(*,'(a,3(a,i5))') 'spffte(): -- ERROR -- ',
     &      'invalid imax =',im,', non-factorized number =',m
          call exit(ierr)
	end subroutine exitm_
	subroutine factor_get_(m,mp,n)
	  implicit none
	  integer,intent(inout) :: m
	  integer,intent(in) :: mp
	  integer,intent(out) :: n
	  integer :: k
	  n=0
	  do while(mod(m,mp)==0)
	    m=m/mp
	    n=n+1
	  enddo
	end subroutine factor_get_
