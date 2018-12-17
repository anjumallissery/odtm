
subroutine advec_mld(fld,u,v,rkmh,rkmu,rkmv,rdxh,rdyh,advec)
    implicit none
    real, intent(in), dimension(:,:,:) :: fld, u, v
    real, intent(in), dimension(:,:  ) :: rkmh, rkmu, rkmv, rdx, rdy
    real, intent(out), dimension(:,:,:) :: advec

    integer :: i, j, k, im, ip, jm, jp
    real :: rtemp1, rtemp2, rtemp

    do k = 1, size(fld,3) !--> Prajeesh
        do j = 2, size(fld,2)-1
            do i = 2, size(fld,2)-1
                im = i-1; ip = i+1
                jm = j-1; jp = j+1
                rtemp1 = (fld(i,j,k) + fld(ip,j,k))/
                       & max(1.0,rkmh(i,j)+rkmh(ip,j))
                rtemp2 = (fld(i,j,k) + fld(im,j,k))/
                       & max(1.0,rkmh(i,j)+rkmh(im,j))
                rtemp = (u(i,j,k)+u(ip,j,k))/
                      & max(1.0,rkmu(i,j)+rkmu(ip,j))
                advec = (rtemp1 - rtemp2) * rtemp * rdx(i,j)
                
                rtemp1 = (fld(i,j,k) + fld(i,jp,k))/
                       & max(1.0,rkmh(i,j)+rkmh(i,jp))
                rtemp2 = (fld(i,j,k) + fld(i,jm,k))/
                       & max(1.0,rkmh(i,j)+rkmh(i,jm))
                rtemp = (v(i,j,k)+v(i,jp,k))/
                      & max(1.0,rkmv(i,j)+rkmv(i,jp))
                advec = advec + (rtemp1 - rtemp2) * rtemp * rdy(i,j)
            end do
        end do
    end do 
    return
end subroutine advec_mld 


subroutine diffuse_mld (fld,rdxu,rkmh,rdyv,rdxh,rdyh,diffuse)
    implicit none
    real, dimension(:,:,:), intent(in) :: fld
    real, dimension(:,:), intent(in) :: rdxu, rkmh, rdyv, rdxh, rdyh
    real, dimension(:,:,:), intent(out) :: diffuse

    integer :: i, j, k
    real :: diffuse
    integer :: ip, im, jp, jm
    real :: rtemp1, rtemp2, rtemp3, rtemp4, rtemp5, rtemp6

    do k = 1, size(fld,3) !--> Prajeesh
        do j = 2, size(fld,2)-1
            do i = 2, size(fld,2)-1
                ip = i+1
                im = i-1
                jp = j+1
                jm = j-1

                rtemp1 = rkmh(ip,j)* (fld(ip,j,k,1) - fld(i,j,k,1) )*rdxu(ip,j)
                rtemp2 = rkmh(im,j)* (fld(i,j,k,1) - fld(im,j,k,1) )*rdxu(i,j)
                rtemp3 = (rtemp1 - rtemp2)*rdxh(i,j) !--->Prajeesh

                rtemp4 = rkmh(i,jp)* (fld(i,jp,k,1) - fld(i,j,k,1))*rdyv(i,jp)
                rtemp5 = rkmh(i,jm)* (fld(i,j,k,1) - fld(i,jm,k,1))*rdyv(i,j)
                rtemp6 = (rtemp4 - rtemp5)*rdyh(i,j)

                diffuse(i,j,k) = rtemp3 + rtemp6
                if (rtemp3 .eq. 0.0 .or. rtemp6 .eq. 0.0) diffuse(i,j,k) = 0.0 !Prajeesh 
            end do
        end do
    end do

    return

end subroutine diffuse_mld

