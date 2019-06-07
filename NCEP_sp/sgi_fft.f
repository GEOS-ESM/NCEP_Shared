!--------------------------------------------------------------
      SUBROUTINE SGI_FFT (isign,nlon,lot,x,ldx,y,ldy,afft)

      implicit none

! !INPUT PARAMETERS:

      integer, intent(in) :: isign
      integer, intent(in) :: nlon
      integer, intent(in) :: lot
      integer, intent(in) :: ldx
      integer, intent(in) :: ldy

      integer  :: mype

! !OUTPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:

      real(8)  :: x(ldx,lot)
      complex(8) :: y(ldy,lot)
      real(8) :: afft (*)

! !DESCRIPTION:
!  This routine calls an FFT from a librrary (-lscs at Goddard).

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm
!  04Oct2004  R. Errico  Module development
!  12July2005 R. Errico  For the QG model
!  29Mar2007  R. Yang  modified
! REMARK: man page about zdfftm/dzfftm is very confusing
!
!-------------------------------------------------------------------------
!
!    local variables

      external dzfftm, zdfftm
      integer                    :: isys(0:1)
      real(8)                    :: work(nlon+2)
      real(8)                   :: scale
!  just follow Ron's dimension value--why multiplying 2 ??
      integer                    :: i

      isys(0)=1

      if (isign == 0) then      ! set a table for later use
         scale=0.d0
        call dzfftm (isign, nlon, lot,scale,x,ldx,y,ldy,
     &               afft, work, isys )
!       write (6,*) 'INITIAL: tablefft '
!...........................
      elseif (isign == 1) then  ! transform zonal coefs to fields
        scale=1.d0
! fill tablefft
        call zdfftm (isign, nlon, lot, scale, y, ldy, x, ldx,
     &               afft, work, isys )
      elseif (isign == -1) then ! transform fields to zonal coefs
        scale=1.d0/real(nlon )
        call dzfftm (isign, nlon, lot, scale, x, ldx, y, ldy,
     &               afft, work, isys )

      endif
      return

      END SUBROUTINE SGI_FFT
