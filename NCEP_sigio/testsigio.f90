       use sigio_module
       type(sigio_head):: head
       type(sigio_data):: data
       call sigio_srohdc(11,'sigf06',head,data,iret)
       call sigio_swohdc(12,'sigf06.new',head,data,iret)
       print *, 'lat, lon, ivs = ', head%latb, head%lonb, head%ivs
       end

