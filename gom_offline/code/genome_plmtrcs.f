subroutine plmtrcs(rl,rc,rr,a,s,n)
implicit none
c
integer,intent(in)  :: n
real,   intent(in)  :: rl(n),rc(n),rr(n),a(n)
real,   intent(out) :: s(n)
c
c**********
c*
c  1) generate slopes for monotonic piecewise linear distribution
c
c  2) input arguments:
c       rl   - left grid spacing ratio
c       rc   - center grid spacing ratio
c       rr   - right grid spacing ratio
c       a    - scalar field zone averages
c       n    - number of zones
c
c  3) output arguments:
c       s    - zone slopes
c
c  4) Tim Campbell, Mississippi State University, September 2002.
c*
c**********
c
integer,parameter :: ic=2, im=1, imax=100
real,parameter :: fracmin=1e-6, dfac=0.5
c
integer i,j
real    sl,sc,sr
real    dnp,dnn,dl,dr,ds,frac
c
c Compute zone slopes
c Campbell Eq(15) -- nonuniform grid
c
s(1)=0.0
do j=2,n-1
  sl=rl(j)*(a(j)-a(j-1))
  sr=rr(j)*(a(j+1)-a(j))
  if (sl*sr.gt.0.) then
    s(j)=sign(min(abs(sl),abs(sr)),sl)
  else
    s(j)=0.0
  endif
enddo
s(n)=0.0
c
c Minimize discontinuities between zones
c Apply single pass discontinuity minimization: Campbell Eq(19)
c
do j=2,n-1
  if(s(j).ne.0.0) then
    dl=-0.5*(s(j)+s(j-1))+a(j)-a(j-1)
    dr=-0.5*(s(j+1)+s(j))+a(j+1)-a(j)
    ds=sign(min(abs(dl),abs(dr)),dl)
    s(j)=s(j)+2.0*ds
  endif
enddo
return
end subroutine plmtrcs
c
c
c> Revision history:
c>
c> Aug  2002 - new routine to put all tracer interactions in one place
c> Dec. 2003 - inforce non-negative bio-tracers
c> Aug. 2005 - interpolate trwall to actual layer structure
