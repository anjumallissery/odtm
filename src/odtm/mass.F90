subroutine mass_balance
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!   A subroutine to ensure mass balance of total thickness
!
!
!
!
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        use size_mod, only : gdx, gdy, h, i, j, k
        use size_mod, only : taup, rkmh, dz
        use size_mod, only :  isc, iec, jsc, jec, km
        use param_mod, only : deg2rad, Re

        use mpp_domains_mod, only : mpp_update_domains, domain2d
        use mpp_mod, only : mpp_sum, mpp_error, NOTE, WARNING, FATAL


        implicit none
        
        real :: rmass_g, rmass_b, rgrid_g

        integer :: im, ip, jm, jp

        character(len=1024) :: msg

        rmass_g = 0; rmass_b = 0; rgrid_g = 0

        do i=isc,iec
           do j=jsc,jec
               if (rkmh(i,j) == 0.) cycle
               do k=1,km-1
                ip = i + 1
                im = i - 1
                jp = j + 1
                jm = j - 1
     rmass_g = rmass_g + h(i,j,k,taup)* deg2rad* &
     Re*Re*(gdx(i) - gdx(im))*          &
     (sin(gdy(j)*deg2rad) - sin(gdy(jm)*deg2rad))

     rmass_b = rmass_b + dz(k)* deg2rad* &
     Re*Re*(gdx(i) - gdx(im))* &
     (sin(gdy(j)*deg2rad) - sin(gdy(jm)*deg2rad))

     rgrid_g = rgrid_g + rkmh(i,j)*deg2rad* &
     Re*Re*(gdx(i) - gdx(im))* &
     (sin(gdy(j)*deg2rad) - sin(gdy(jm)*deg2rad))

               enddo
           enddo
        enddo 
        call mpp_sum(rmass_g)
        call mpp_sum(rmass_b)
        call mpp_sum(rgrid_g)
!        write(msg,*) "Global Sum of mass difference", &
!        rmass_g-rmass_b
!        call mpp_error(NOTE,msg)

        do i=isc,iec
          do j=jsc,jec
            if (rkmh(i,j) == 0.) cycle
            do k=1,km-1
                h(i,j,k,taup) = h(i,j,k,taup) - &
        ( (rmass_g-rmass_b)/ (rgrid_g*(km-1)) ) *rkmh(i,j)
            enddo
          enddo
        enddo

        return
end subroutine mass_balance
