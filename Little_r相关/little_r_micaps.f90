PROGRAM MAIN
IMPLICIT NONE
CHARACTER*40:: fdate,UTC_date
CHARACTER*40:: fdate_char,UTC_date_char
CHARACTER*40:: cTemp
CHARACTER*80:: plot_path
CHARACTER*40:: id
CHARACTER*80:: st_name='AGADEZ / NIGER'
CHARACTER*80:: plat_form='FM-35 TEMP'
CHARACTER*80:: source='GTS (ROHK) ULNR01 DRRN 051100'
INTEGER i,j,k,qc
REAL g,dew
INTEGER  yr,mon,date,hr,level,stat_num  
REAL lev(11)
INTEGER,ALLOCATABLE:: stat_id(:,:),stat_class(:,:)											  
REAL,ALLOCATABLE:: stat_elev(:,:),height(:,:),temp(:,:),longi(:,:),lati(:,:),temp_dew(:,:),win_dir(:,:),V(:,:)
REAL,ALLOCATABLE:: ht(:,:) 
REAL a(13)
INTEGER b(13)
LOGICAL is_sound,bogus,discard
CHARACTER*40:: stat_name='aaa'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

fdate='${BJ_time}'   !beijing time
UTC_date='${UTC}' !UTC
plot_path='${MPS_PATH}' !PATH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


qc=0
is_sound=.true.
bogus=.false.
discard=.false.
fdate_char='20'//trim(fdate)//'0000'
UTC_date_char='20'//trim(UTC_date)//'0000'
a=-888888.0
b=0
DATA lev/1000.0,925.0,850.0,700.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0/

OPEN (101,FILE=trim(plot_path)//'plot/850/'//trim(adjustl(fdate))//'.000')
READ(101,*)
READ(101,*) yr,mon,date,hr,level,stat_num
CLOSE(101)

ALLOCATE (stat_id(stat_num,size(lev)))
ALLOCATE (longi(stat_num,size(lev)))
ALLOCATE (lati(stat_num,size(lev)))
ALLOCATE (stat_elev(stat_num,size(lev)))
ALLOCATE (stat_class(stat_num,size(lev)))
ALLOCATE (height(stat_num,size(lev)))
ALLOCATE (temp(stat_num,size(lev)))
ALLOCATE (temp_dew(stat_num,size(lev)))
ALLOCATE (win_dir(stat_num,size(lev)))
ALLOCATE (V(stat_num,size(lev)))
ALLOCATE (ht(stat_num,size(lev)))

DO j=1,size(lev)
	WRITE(ctemp,'(i4)') int(lev(j))
	OPEN(101,FILE=trim(plot_path)//'plot/'//trim(adjustl(ctemp))//'/'//trim(adjustl(fdate))//'.000')
		READ(101,*)
		READ(101,*)
		DO i=1,stat_num
			READ(101,*) stat_id(i,j),longi(i,j),lati(i,j),stat_elev(i,j),stat_class(i,j),&
									height(i,j),temp(i,j),temp_dew(i,j),win_dir(i,j),V(i,j)
		ENDDO
	CLOSE(101)
ENDDO

OPEN(101,FILE='obs.20'//trim(adjustl(UTC_date))//'',ACCESS='append')
	DO i=1,stat_num
		g=9.7803185*(1+0.005278895*(sin(lati(i,1)))**2+0.000023462*(sin(lati(i,1)))**4)
		ht(i,:)=height(i,:)*10*9.80665/g			
		WRITE(ctemp,'(i5)') stat_id(i,1)	
		WRITE(id,'(i5)') stat_id(i,1)
		WRITE(101,'(2f20.5,4A40,f20.5,5I10,3L10,2I10,A20,13(f13.5,I7))') longi(i,1),lati(i,1),adjustl(id),&
									st_name,plat_form,source,stat_elev(i,1),-888888,-888888,-888888,&
									i,-888888,is_sound,bogus,discard,-888888,-888888,trim(UTC_date_char),&  
									-888888.,0,((a(k),b(k)),k=1,12)
		DO j=1,size(lev)
				dew=(temp(i,j)+273.16-temp_dew(i,j))
				IF(temp(i,j)==9999.) temp(i,j)=-888888.-273.16
				IF((temp_dew(i,j)==9999.).or.(temp(i,j)==9999.)) dew=-888888.0
				IF(win_dir(i,j)==9999.) win_dir(i,j)=-888888.
				IF(V(i,j)==9999.) V(i,j)=-888888.
				WRITE(101,'(10(F13.5,I7))') lev(j)*100,qc,ht(i,j),qc,(temp(i,j)+273.16),qc,dew,qc,&
																		V(i,j),qc,win_dir(i,j),qc,((a(k),b(k)),k=1,4)
		ENDDO
		WRITE(101,'(10(F13.5,I7))') -777777.,0,-777777.,0,((a(k),b(k)),k=1,8)
		WRITE(101,*) 50,0,0
	ENDDO
CLOSE(101)

DEALLOCATE(stat_id)
DEALLOCATE(longi)
DEALLOCATE(lati)
DEALLOCATE(stat_elev)
DEALLOCATE(stat_class)
DEALLOCATE(height)
DEALLOCATE(temp)
DEALLOCATE(temp_dew)
DEALLOCATE(win_dir)
DEALLOCATE(V)

END 