
!***** module containing functions used in main program

module strings
implicit none
contains
  
function findxiplus(r,m)
double precision,intent(in)::r,m
double precision::findxiplus
findxiplus=m*sqrt((1.d0-r)/(r*(m**2)-1.d0))
end function findxiplus


function ustring(str,u,N)
integer, intent(in)::N,str(1:N)
double precision, intent(in)::u(1:N)
integer::i
double precision::ux=0.d0
double precision::ustring


ux=0.d0
do i=1,N
ux=ux+1.d0*str(i)*u(i)
enddo

ustring=ux

end function ustring

 
function vstring(str,v,N)
integer, intent(in)::N,str(1:N)
double precision, intent(in)::v(1:N)
integer::i
double precision::vx=0.d0
double precision::vstring

vx=0.d0
do i=1,N
vx=vx+1.d0*str(i)*v(i)
enddo

vstring=vx
end function vstring

 
function lnfstring(str,u,v,y,N)
integer, intent(in)::N,str(1:N)
double precision, intent(in)::u(1:N),v(1:N),y
integer::i
double precision::ux=0.d0,vx=0.d0,z=0.d0,z1=0.d0
double precision::lnfstring


ux=0.d0
vx=0.d0
do i=1,N
ux=ux+1.d0*str(i)*u(i)
vx=vx+1.d0*str(i)*v(i)
enddo

z=y-vx
z1=exp(z)
lnfstring=ux-log(1.d0+(z1)**2)

end function lnfstring


function intxy(r1,m1,r2,m2)
double precision,intent(in)::r1,m1,r2,m2
double precision::intxy
intxy=sqrt((r1-r2)/(r2/m1**2-r1/m2**2))
end function intxy


function bitstr(i,N)
  integer,intent(in)::N,i
  integer::k,j,bitstr(1:N)
    
  do j=1,N
     k=j
     bitstr(j)=(mod(i,2**(k))-mod(i,2**(k-1)))
     bitstr(j)=bitstr(j)/(2**(k-1))
  enddo


end function bitstr



function insertbit(str,i,ii,N)
  integer,intent(in)::N,i,ii
  integer::str(1:N),insertbit(1:N),k,j

  insertbit=str
  
  do j=2,i
     k=str(j)
     insertbit(j-1)=k
  enddo

  insertbit(i)=ii

end function insertbit


function dec(str,N)
  integer,intent(in)::N,str(1:N)
  integer::i,dec
!common N

  dec=0
  do i=1,n
     dec=dec+str(i)*(2**(i-1))
  enddo
  
end function dec

  function dec2(str,N)
  integer,intent(in)::N,str(1:N)
  integer::i,dec2
!common N

  dec2=0
  do i=1,n
     dec2=dec2+str(n-i+1)*(2**(i-1))
  enddo


end function dec2



function neighbor(str,i,N)
    integer,intent(in)::N,i,str(1:N)
    integer::neighbor(1:N)

    neighbor=str
    neighbor(i)=abs(1-str(i))


  end function neighbor

function bitflip(str,i,N)
    integer,intent(in)::N,i,str(1:N)
    integer::bitflip(1:N)

    bitflip=str
    bitflip(i)=abs(1-str(i))

  end function bitflip


    function strsize(str,N)
    integer,intent(in)::N
    integer::str(1:N),i
    double precision::strsize
    strsize=0.d0
    do i=1,N
       strsize=strsize+1.d0*str(i)
    enddo

  end function strsize


	function hamdist(g1,g2,N)
	integer,intent(in)::N,g1,g2
    	integer::str1(1:N),str2(1:N),j
    	double precision::hamdist
   	
	hamdist=0.d0
	str1=bitstr(g1,N)
	str2=bitstr(g2,N)	
	do j=1,n
	hamdist=hamdist+1.d0*abs(1.d0*str1(j)-1.d0*str2(j))
	enddo
  
	end function hamdist
    
       
FUNCTION ran2(idum)
IMPLICIT NONE
INTEGER, PARAMETER :: im1=2147483563,im2=2147483399,imm1=im1-1,&
     ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,&
     ir2=3791,ntab=32,ndiv=1+imm1/ntab
DOUBLE PRECISION , PARAMETER ::   am=1.d0/im1,eps=1.d-14,rnmx=1.d0-eps
DOUBLE PRECISION :: ran2
INTEGER, DIMENSION(ntab) :: iv
INTEGER :: idum,idum2,j,k,iy

save iv,iy,idum2
data idum2/123456789/,iv /ntab*0/,iy /0/
      
if(idum.le.0) then
  idum=max(-idum,1)
  idum2=idum
  do j=ntab+8,1,-1
    k=idum/iq1
    idum=ia1*(idum-k*iq1)-ir1*k
    if(idum.lt.0) idum=idum+im1
    if(j.le.ntab) iv(j)=idum
   end do
   iy=iv(1)
endif

k=idum/iq1
idum=ia1*(idum-k*iq1)-ir1*k
if(idum.lt.0) idum=idum+im1
k=idum2/iq2
idum2=ia2*(idum2-k*iq2)-ir2*k

if(idum2.lt.0) idum2=idum2+im2

j=1+iy/ndiv
iy=iv(j)-idum2
iv(j)=idum
if (iy.lt.1)iy=iy+imm1
ran2=min(am*iy,rnmx)

END FUNCTION ran2


end module


!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM
!**********************************************MAIN PROGRAM


program main
  
  use strings
  implicit none
  
  !***** model parameters
  integer(kind=4),parameter::N=100,rea=10000
  double precision, parameter::mu=0.d0
  !***** simultaion paremeters
  integer(kind=4),parameter::tup=10000
  !*****
  
  integer(kind=4),parameter::ns=200,npop=1000,xnmup=50
  double precision,parameter::vavg=.797,uavg=-0.518
  integer::seed=-8591077
  
  double precision::reacount(1:xnmup)=0.d0,xpts(1:xnmup)=0.d0,foldinc=0.d0,ubv=0.d0
  integer(kind=4)::i,j,jk,ii,m,ik,rcntr=0,ini=0,flg=0,xnm=1
  double precision::mnum(1:xnmup,0:tup)=0.d0,fparms(1:xnmup,0:tup,1:3)=0.d0
  double precision::tradeoff(1:xnmup,0:tup)=0.d0,slcoeffx(1:xnmup,0:tup)=0.d0
  double precision::xiplus(1:N)=0.d0,ximinus(1:N)=0.d0,x0plus=0.d0,x0minus=0.d0
  integer(kind=4)::str(1:N)=0,npk=0,nstr(1:N)=0,str1(1:N)=0,str2(1:N)=0,fnbrstr(1:N)=0,strx(1:N)=0
  integer(kind=4)::flownbrstr(1:N)=0,nbr,sflip=1,tlast=0,oldsize=0,newsize=0,flg2=0,timestep=0
  double precision::fcurrent=0.d0,fcompare=0.d0
  double precision::selc=0.d0,xx=0.d0,xdmax=0.d0
 
  double precision::rsingle(1:N)=1.d0,micsingle(1:N)=1.d0
  double precision::usingle(1:N)=0.d0,vsingle(1:N)=0.d0,lnr=0.d0
  double precision::x1,x2,x
  double precision::sltemp=0.d0,tdtemp=0.d0,tdreacounts(1:xnmup,0:tup)=0.d0
  double precision::xin=0.d0,xinc=2.d0*vavg,fcurrentlow=0.d0,fori=0.d0,fittercounts=0.d0
 
  double precision::ftns(1:xnmup)=0.d0
  double precision::vsig(1:xnmup)=0.d0,vsig2(1:xnmup)=0.d0,usig(1:xnmup)=0.d0
  double precision::ftns2(1:xnmup)=0.d0,usig2(1:xnmup)=0.d0,wlength(1:xnmup)=0.d0,wdiff(1:xnmup)=0.d0
  double precision::uvsorted(1:N)=0.d0,uvsavg(1:N)=0.d0,uvtemp=0.d0,uvscum(1:N)=0.d0
  double precision::pi=acos(-1.d0)

  integer(kind=4)::conccounter=0,wl=0,fnbrsaved(1:N)=0
  double precision::irrevup(1:xnmup)=0.d0,irrevdown(1:xnmup)=0.d0,mindist=0.d0,lnm=0.d0,dirpaths(1:xnmup)=0.d0



30 format(f48.9,3x,f16.7)
character(len=12)::flnm

!***** setting concentration values
xin=0.25*N*vavg
xdmax=1.d0*N*vavg
xinc=0.25d0*vavg*N
xnm=int((xdmax-xin)/xinc)+1

do ik=1,xnm
xpts(ik)=xin+xinc*(1.d0*ik-1.d0)
enddo

!***** realization loop begins
do rcntr=1,rea,1
reacount=0.d0
	!*** setting u_i, v_i
   do i=1,N
     x1=ran2(seed)
     x2=ran2(seed)
     vsingle(i)=(1.d0*abs(1.d0*sqrt(-2.d0*log(x1))*sin(2.d0*pi*x2)))
     flg=0
      do while(flg==0)
      	x1=ran2(seed)
      	x2=ran2(seed)
      	usingle(i)=-mu-1.d0*abs(1.d0*sqrt(-2.d0*log(x1))*cos(2.d0*pi*x2))
	if (-1.d0*usingle(i)<2.d0*vsingle(i)) flg=1
     enddo
     lnm=lnm+vsingle(i)/(1.d0*rea*N)
	lnr=lnr+usingle(i)/(1.d0*rea*N)
	ubv=ubv+usingle(i)/vsingle(i)*(1.d0/(1.d0*rea*N))
  enddo

!*** sorting u/v
do i=1,N
uvsorted(i)=-usingle(i)/vsingle(i)
enddo

do i=1,N
do j=2,N
if (uvsorted(j)<uvsorted(j-1)) then 
uvtemp=uvsorted(j)
uvsorted(j)=uvsorted(j-1)
uvsorted(j-1)=uvtemp
endif
enddo
enddo

do i=1,N
uvsavg(i)=uvsavg(i)+uvsorted(i)/(1.d0*rea)
enddo

!***** concentration loop begins
do conccounter=1,xnm
x=xpts(conccounter)


ini=0

str=0

wl=0
oldsize=0
flg2=0

timestep=0
tlast=-1


mnum(conccounter,timestep)=mnum(conccounter,timestep)+1.d0*strsize(str,N)
fparms(conccounter,timestep,1)=fparms(conccounter,timestep,1)+lnfstring(str,usingle,vsingle,x,N)
fparms(conccounter,timestep,2)=fparms(conccounter,timestep,2)+ustring(str,usingle,N)
fparms(conccounter,timestep,3)=fparms(conccounter,timestep,3)+vstring(str,vsingle,N)
tdtemp=0.d0

do jk=1,N
tdtemp=tdtemp+1.d0*str(jk)*(abs(usingle(jk))/vsingle(jk))
enddo
if (strsize(str,N)>0) then
tradeoff(conccounter,timestep)=tradeoff(conccounter,timestep)+tdtemp/(1.d0*strsize(str,N))
tdreacounts(conccounter,timestep)=tdreacounts(conccounter,timestep)+1.d0
endif

!*********************** adaptive walk begins
do
fcurrent=lnfstring(str,usingle,vsingle,x,N)
fori=fcurrent
fcompare=0.d0
fcurrentlow=1.1d0
flownbrstr=str
flg=0

fittercounts=0.d0
sltemp=0.d0
do ik=1,N
nstr=neighbor(str,ik,N)
fcompare=lnfstring(nstr,usingle,vsingle,x,N)
if (fcompare.ge.fcurrent) then
selc=exp(fcompare-fcurrent)-1.d0
sltemp=sltemp+1.d0-exp(-2.d0*selc)
fittercounts=fittercounts+1.d0
endif
enddo
if (fittercounts>0.d0) then 
slcoeffx(conccounter,timestep)=slcoeffx(conccounter,timestep)+sltemp/fittercounts
endif

do ik=1,N
nstr=neighbor(str,ik,N)
fcompare=lnfstring(nstr,usingle,vsingle,x,N)
if (fcompare.ge.fcurrent) then
flg=1
fnbrstr=nstr
fcurrent=fcompare
endif

if ((fcompare.ge.fori).and.(fcompare.le.fcurrentlow)) then
flownbrstr=nstr
fcurrentlow=fcompare
endif
enddo
 
   if (flg==0) then 
	selc=0.d0
	 exit
	endif

if (tlast.ge.tup) exit

fnbrsaved=fnbrstr

fcurrent=lnfstring(str,usingle,vsingle,x,N)

!***** adaptive step begins
do 
	   timestep=timestep+1
       sflip=ran2(seed)*1.d0*N
	   sflip=sflip+1
       strx=str
	   nstr=neighbor(str,sflip,N)
         
           fcompare=lnfstring(nstr,usingle,vsingle,x,N)
           selc=exp(fcompare-fcurrent)-1.d0
  
   	   xx=ran2(seed)
	     if (xx<1.d0-exp(-2.d0*selc)) then 
	     !if (selc>0.d0) then
		 str=nstr
		tlast=timestep
                if (timestep.le.tup) then
		mnum(conccounter,timestep)=mnum(conccounter,timestep)+1.d0*strsize(str,N)
  		fparms(conccounter,timestep,1)=fparms(conccounter,timestep,1)+lnfstring(str,usingle,vsingle,x,N)
  		fparms(conccounter,timestep,2)=fparms(conccounter,timestep,2)+ustring(str,usingle,N)
  		fparms(conccounter,timestep,3)=fparms(conccounter,timestep,3)+vstring(str,vsingle,N)
	    tdtemp=0.d0
	    do jk=1,N
	    tdtemp=tdtemp+1.d0*str(jk)*(abs(usingle(jk))/vsingle(jk))
	    enddo
	    if (strsize(str,N)>0) then
	    tradeoff(conccounter,timestep)=tradeoff(conccounter,timestep)+tdtemp/(1.d0*strsize(str,N))
	    tdreacounts(conccounter,timestep)=tdreacounts(conccounter,timestep)+1.d0
	    endif
	        
            endif
		exit       
	   endif
 
if (timestep.le.tup) then
 mnum(conccounter,timestep)=mnum(conccounter,timestep)+1.d0*strsize(str,N)
fittercounts=0.d0
sltemp=0.d0
do ik=1,N
nstr=neighbor(str,ik,N)
fcompare=lnfstring(nstr,usingle,vsingle,x,N)
if (fcompare.ge.fcurrent) then
selc=exp(fcompare-fcurrent)-1.d0
sltemp=sltemp+1.d0-exp(-2.d0*selc)
fittercounts=fittercounts+1.d0
endif
enddo
if (fittercounts>0.d0) then 
slcoeffx(conccounter,timestep)=slcoeffx(conccounter,timestep)+sltemp/fittercounts
endif
fparms(conccounter,timestep,1)=fparms(conccounter,timestep,1)+lnfstring(str,usingle,vsingle,x,N)
fparms(conccounter,timestep,2)=fparms(conccounter,timestep,2)+ustring(str,usingle,N)
fparms(conccounter,timestep,3)=fparms(conccounter,timestep,3)+vstring(str,vsingle,N)
tdtemp=0.d0
do jk=1,N
tdtemp=tdtemp+1.d0*str(jk)*(abs(usingle(jk))/vsingle(jk))
enddo
if (strsize(str,N)>0.d0) then
tradeoff(conccounter,timestep)=tradeoff(conccounter,timestep)+tdtemp/(1.d0*strsize(str,N))
tdreacounts(conccounter,timestep)=tdreacounts(conccounter,timestep)+1.d0
endif

endif
!print*,rcntr,timestep!,strsize(str,N)
if (timestep>tup) exit
       
enddo
!***** adaptive step ends


if (timestep>tup) exit
enddo
!***** adaptive walk ends

if (timestep<tup) then 

do i=timestep+1,tup
  mnum(conccounter,i)=mnum(conccounter,i)+1.d0*strsize(str,N)
  fparms(conccounter,i,1)=fparms(conccounter,i,1)+lnfstring(str,usingle,vsingle,x,N)
  fparms(conccounter,i,2)=fparms(conccounter,i,2)+ustring(str,usingle,N)
  fparms(conccounter,i,3)=fparms(conccounter,i,3)+vstring(str,vsingle,N)
  tdtemp=0.d0
    do jk=1,N
	tdtemp=tdtemp+1.d0*str(jk)*(abs(usingle(jk))/vsingle(jk))
	enddo
	   if (strsize(str,N)>0.d0) then
	tradeoff(conccounter,i)=tradeoff(conccounter,i)+tdtemp/(1.d0*strsize(str,N))
	tdreacounts(conccounter,i)=tdreacounts(conccounter,i)+1.d0
	endif
enddo

endif

enddo
!***** concentration loop ends

enddo
!***** realization loop ends


do i=1,N
uvscum(i)=0.d0
do j=1,i
uvscum(i)=uvscum(i)+uvsavg(j)
enddo
uvscum(i)=uvscum(i)/(1.d0*i)
enddo

open (unit = 10, file = "output.d")

do ik=1,xnm
do wl=0,tup
write(10,*) wl,mnum(ik,wl)/(1.d0*rea),fparms(ik,wl,1)/(1.d0*rea),fparms(ik,wl,2)/(1.d0*rea),&
fparms(ik,wl,3)/(1.d0*rea),slcoeffx(ik,wl)/(1.d0*rea),tradeoff(ik,wl)/(1.d0*tdreacounts(ik,wl))
enddo

write(10,*) " "


enddo

!write(10,*) lnr, lnm

close(10)


end program main




  
  

  
