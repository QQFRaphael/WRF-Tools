program read_surf
implicit none
character*6 :: year,mon,day,hr
integer :: i,j,k,ii,jj ,iseq_num
integer,parameter :: total = 254160
character*10 mdate
character*16 infile_sat(2) 


logical bogus
real,dimension(720,353) ::speed,dir1
real  lat(total),lon(total),z(total),ter(total)    
real  p(total),slp(total),spd(total),dir(total)  
real  t(total),td(total)    

data infile_sat /'speed_200481100','dir_200481100'/ 
k=1
open(10,file='F:\wl\quikscat\readwind\'//infile_sat(1),status='old' )    
open(11,file='F:\wl\quikscat\readwind\'//infile_sat(2),status='old')     

read(10,'(720e16.7)') ((speed(i,j),i=1,720),j=1,353)

read(11,'(720e16.7)') ((dir1(i,j),i=1,720),j=1,353)
spd =reshape( speed ,(/254160/))
dir =reshape( dir1 ,(/254160/))


do i=1,total
ter(i)=10   !高度为10米
ii=int(i/720.)
jj=i-ii*720
lat(i)= ii*0.5-88.
lon(i)= jj*0.5
  
call write_obs(p(i),z(i),t(i),td(i),spd(i),dir(i),slp(i),&
&              ter(i),lat(i),lon(i),mdate,k,&
&        '99001  Maybe more site info             ',&
&        'SURFACE DATA FROM ????????? SOURCE       ',&
&        'FM-12 SYNOP                             ',&
&        '                                        ',&
&              bogus , i, 20 ) 

enddo

write(*,*)">> END READ SURFACE DATA "

close(10)
close(11)
end
 
subroutine write_obs(p,z,t,td,spd,dir,slp,ter,&
& lat,lon,mdate,kx,string1,string2,string3,string4,&
& bogus,id,iunit)
implicit none
integer k,iunit,kx
real p(kx),z(kx),t(kx),td(kx),spd(kx),dir(kx)
real lat,lon,ter,slp
integer id
character *20 date_char
character *40 string1,string2,string3,string4
character *84 rpt_format
character *22 meas_format
character *14 end_format
character *10 mdate
character *20 outfile
logical bogus
rpt_format ='(2f20.5 , 2a40 , 2a40 , 1f20.5 ,' &
& // '5i10 , 3L10 ,2i10 , a20 ,13( f13.5 , i7 ))'
meas_format ='( 10( f13.5 , i7 ) )'
end_format = ' ( 3( i7 ) ) '
mdate='2004081100'
write(date_char(7:16),fmt='(a10)')mdate
date_char(17:20)='0000'
date_char(1:6)='      '
outfile='surface_obs'//mdate(7:10)
open(iunit,file='F:\wl\little-r\'//outfile,form='formatted')
write(iunit,err=100,fmt=rpt_format) &
&  lat,lon,string1,string2,string3,string4, &
&  ter,kx*6,0,0,id,0,.false.,bogus,.false., &
&  -888888,-888888,date_char, &
& -888888.,0,-888888.,0, -888888.,0, -888888.,0, -888888.,0, &
&  -888888.,0, &
&  -888888.,0, -888888.,0

!do k=1,kx
write(iunit,err=100,fmt=meas_format) &
&  -888888.,0,-888888.,0,-888888.,0,-888888.,0, &
&  spd,0,dir,0, &
&  -888888.,0, -888888.,0, -888888.,0, -888888.,0

!enddo 
write(iunit,err=100,fmt=meas_format) &
&  -777777.,0, -777777.,0, float(kx),0, &
&  -888888.,0, -888888.,0, -888888.,0, &
&  -888888.,0, -888888.,0, -888888.,0, &
&  -888888.,0
write(iunit,err=100,fmt=end_format)kx*6,0,0
100continue 
!print*,'troubles writing a sounding'
!stop 100
end
