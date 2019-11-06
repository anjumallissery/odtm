        function h_mix (i1,j1,k1)
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

        use size_mod, only : h, i, j, k, rkmh
        use size_mod, only : taum

        use mpp_domains_mod, only : mpp_update_domains, domain2d
        use mpp_mod, only : mpp_sum, mpp_error, NOTE, WARNING, FATAL


        implicit none
        
        real :: h_mix, rKH

        integer :: im, ip, jm, jp, i1, j1, k1

        character(len=1024) :: msg

        h_mix = 0.0

                ip = i + 1
                im = i - 1
                jp = j + 1
                jm = j - 1

        rKH = 10.0e-4

        if (k .gt. 1 .and. rkmh(i1,j1) .eq. 1) then
         h_mix =                               &
        ( ((h(i,j,k-1,taum) - h(i,j,k,taum))/  &
        (h(i,j,k-1,taum)*0.5 + h(i,j,k,taum)*0.5)) - &
          ((h(i,j,k,taum) - h(i,j,k+1,taum))/  &
        (h(i,j,k,taum)*0.5 + h(i,j,k+1,taum)*0.5)) )*rKH / &
         max(1.0,(h(i,j,k,taum)))
        endif

        return
end function h_mix
