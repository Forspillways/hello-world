# hello-world
!Just a handsome boy
!Hi world
!My name is max I will be you dady of the Hydraulics!
program D2Q9LBM
implicit none
include "head.inc"
integer t_max,time,k,x,y,i,j,h,flag
integer obst(lx,ly)
real*8 u_x(lx,ly),u_y(lx,ly),Pixel(lx,ly)
real*8 rho(lx,ly),a1,tw,con(lx,ly)
real visc,s3,s5,Ds,Dsp
real*8 upx(lx,ly),u(lx,ly),upy(lx,ly),u0(lx,ly),v0(lx,ly),stM(0:8,0:8),tM(0:8,0:8),S(0:8),M(0:8,0:8)
real*8 ff(0:8,lx,ly),geq(0:8,lx,ly),g(0:8,lx,ly),sum0,sum1 !u0,u1,v0,v1分别为上一时步和当前时步的x,y速度
data xc/0.d0, 1.d0, 0.d0, -1.d0, 0.d0, 1.d0, -1.d0, -1.d0, 1.d0/,&
&yc/0.d0, 0.d0, 1.d0, 0.d0, -1.0d0, 1.d0, 1.d0, -1.d0, -1.d0/
data ex/0, 1, 0, -1, 0, 1, -1, -1, 1 /,&
&ey/0, 0, 1, 0, -1, 1, 1, -1, -1/
data opp/3,4,1,2,7,8,5,6/

cc = 1.d0
c_squ = cc *cc / 3.d0
!D2Q9 weight
t_k(0) = 4.d0 / 9.d0
do k =1,4
   t_k(k) = 1.d0 / 9.d0
end do
do k =5,8
   t_k(k) = 1.d0 / 36.d0
end do

!D2Q5 weight
gt_k(0)=4.d0 / 9.d0!1.d0/3.d0
do k =1,4
   gt_k(k)=1.d0 / 9.d0!1.d0/6.d0
end do
do k =5,8
   gt_k(k) =1.d0 / 36.d0
end do

!MRT变量
M(0,:)=(/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
M(1,:)=(/-4.0,-1.0,-1.0,-1.0,-1.0,2.0,2.0,2.0,2.0/)
M(2,:)=(/4.0,-2.0,-2.0,-2.0,-2.0,1.0,1.0,1.0,1.0/)
M(3,:)=(/0.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,1.0/)
M(4,:)=(/0.0,-2.0,0.0,2.0,0.0,1.0,-1.0,-1.0,1.0/)
M(5,:)=(/0.0,0.0,1.0,0.0,-1.0,1.0,1.0,-1.0,-1.0/)
M(6,:)=(/0.0,0.0,-2.0,0.0,2.0,1.0,1.0,-1.0,-1.0/)
M(7,:)=(/0.0,1.0,-1.0,1.0,-1.0,0.0,0.0,0.0,0.0/)
M(8,:)=(/0.0,0.0,0.0,0.0,0.0,1.0,-1.0,1.0,-1.0/)
a1=1.0/36.0
tM(0,:)=(/4.0*a1,-4.0*a1,4.0*a1,0,0,0,0,0,0/)
tM(1,:)=(/4.0*a1,-a1,-2.0*a1,6.0*a1,-6.0*a1,0,0,9.0*a1,0/)
tM(2,:)=(/4.0*a1,-a1,-2.0*a1,0,0,6.0*a1,-6.0*a1,-9.0*a1,0/)
tM(3,:)=(/4.0*a1,-a1,-2.0*a1,-6.0*a1,6.0*a1,0,0,9.0*a1,0/)
tM(4,:)=(/4.0*a1,-a1,-2.0*a1,0,0,-6.0*a1,6.0*a1,-9.0*a1,0/)
tM(5,:)=(/4.0*a1,2.0*a1,a1,6.0*a1,3.0*a1,6.0*a1,3.0*a1,0,9.0*a1/)
tM(6,:)=(/4.0*a1,2.0*a1,a1,-6.0*a1,-3.0*a1,6.0*a1,3.0*a1,0,-9.0*a1/)
tM(7,:)=(/4.0*a1,2.0*a1,a1,-6.0*a1,-3.0*a1,-6.0*a1,-3.0*a1,0,9.0*a1/)
tM(8,:)=(/4.0*a1,2.0*a1,a1,6.0*a1,3.0*a1,-6.0*a1,-3.0*a1,0,-9.0*a1/)

write (6,*) 'lattice size lx = ',lx
write (6,*) 'ly = ',ly

call read_parameters(t_max,visc,Ds)
s3=1.0
s5=1.0
!S(:)=(/1.0,1.4,1.4,s3,1.2,s5,1.2,2.0/(1+6.0*visc),2.0/(1+6.0*visc)/)
S(:)=(/1.0,1.1,1.1,1.0,8.0*(2.0-2.0/(1.0+6.0*visc))/(8.0-2.0/(1.0+6.0*visc)),1.0,8.0*(2.0-2.0/(1.0+6.0*visc))/(8.0-2.0/(1.0+6.0*visc)),2.0/(1.0+6.0*visc),2.0/(1.0+6.0*visc)/)
do i=0,8
    do j=0,8
        stm(i,j)=tM(i,j)*S(j)
    end do
end do

call read_obstacles(obst)
call init_density(obst,u_x,u_y,rho,ff,g,con,Pixel)

error=1.0  !给误差赋初值
flag=0
do time = 1, t_max
  ! if(error .ge. 0.0001 .or. flag==1) then 
    !   u0=u_x  !计算误差上一步的速度
    !   v0=u_y  !计算误差上一步的速度     
    !   call inletoutlet(obst,ff,u_x,u_y) 
   !    call getuv(u_x,u_y,u,rho,ff,obst)   
    !   call collision(obst,u_x,u_y,rho,ff,tM,M,S,stM)
    !   call boundaries(u_x,u_y,obst,rho,ff)
    !   call stream(obst,ff)
      !计算误差
    !   sum0=0.0
    !   sum1=0.0
    !   do x = 1, lx
    !      do y = 1, ly
    !         sum0=sum0+(u_x(x,y)-u0(x,y))*(u_x(x,y)-u0(x,y))+(u_y(x,y)-v0(x,y))*(u_y(x,y)-v0(x,y))
    !         sum1=sum1+u_x(x,y)*u_x(x,y)+u_y(x,y)*u_y(x,y)
    !      end do
    !   end do
    !   error=sqrt(sum0)/sqrt((sum1+1e-30)) !分母加上数字，防止分母为0
   ! end if
    !扩散
   ! if(error .lt. 0.0001) then
    !   write(*,*) error
       call g_inletoutlet(obst,con,g) 
       call getcon(con,obst,Ds,g,Pixel,flag)
       call g_collision(u_x,u_y,con,g,geq)
       call g_boundaries(obst,g,geq)
       call g_stream(g)
  !  end if
    !输出结果
    if(mod(time, Nwri) .eq. 0 .or. time .eq. 1) then
       write(*,*) time
       call write_results2(obst,rho,u_x,u_y,u,time,con,Pixel)
    end if

end do

write (6,*) '======= end ========'
16 format(i8,f15.8)
end
