       use sfcio_module
       implicit none
       type(sfcio_head):: head
       type(sfcio_data):: data
       integer ii,jj,iret
       real(4),allocatable,dimension(:,:)::fld
       call sfcio_srohdc(11,'sfcf06',head,data,iret)
       print *, 'idate ', head%fhour,head%idate
       print *, 'lat, lon, ivs = ', head%latb, head%lonb, head%ivs
       call sfcio_swohdc(12,'sfcf06.new',head,data,iret)
       call baopenwt(22,'ncepsfc.grd',iret)
       allocate(fld(size(data%vtype,1),size(data%vtype,2)))
       call sp2np_(data%vfrac)
       call wryte(22,4*head%latb*head%lonb,fld)
       call sp2np_(data%vtype)
       print *, 'min,max vtype: ',minval(data%vtype), maxval(data%vtype)
       call wryte(22,4*head%latb*head%lonb,fld)
       call sp2np_(data%stype)
       print *, 'min,max vtype: ',minval(data%stype), maxval(data%stype)
       call wryte(22,4*head%latb*head%lonb,fld)
       call sp2np_(data%tsea)
       call wryte(22,4*head%latb*head%lonb,fld)
       call baclose(22,iret)
       deallocate(fld)
       contains
       subroutine sp2np_ (fldin)
       real(4) fldin(:,:)
       integer i,j
       integer im,jm
       im=size(fld,1)
       jm=size(fld,2)
       fld(:, 1)=fldin(:,jm)               ! add South Pole points
       do j=2,jm-1
                                           ! for j :=    2,   3, ..., jm-2,jm-1
                                           !  jm-j == jm-2,jm-3, ...,    2,   1
          fld(:,j)=fldin(:,jm-j)
       end do
       fld(:,jm)=fldin(:,   1)             ! add North Pole points
       end subroutine sp2np_
       end

