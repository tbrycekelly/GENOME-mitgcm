subroutine plmtrcx(pt, s,ss,ki,ks)
implicit none
c
integer ki,ks
real    pt(ki+1),s(ki,ks),ss(ki,ks)
c
c**********
c*
c  1) generate a monotonic PLM interpolation of a layered field
c
c  2) input arguments:
c       pt    - layer interface thicknesses (non-zero)
c       s     - scalar fields in layer space
c       ki    - 1st dimension of s (number of layers)
c       ks    - 2nd dimension of s (number of fields)
c
c  3) output arguments:
c       ss    - scalar field slopes for PLM interpolation
c
c  4) except at data voids, must have:
c           pi(   1) == zero (surface)
c           pi( l+1) >= pi(:,:,l)
c           pi(ki+1) == bathymetry
c
c  5) Tim Campbell, Mississippi State University, September 2002.
c*
c**********
c
integer l
real    ql(ki),qc(ki),qr(ki)
c
!compute grid spacing ratios for slope computations
ql(1)=0.0
qc(1)=0.0
qr(1)=0.0
do l=2,ki-1
  ql(l)=2.0*pt(l)/(pt(l-1)+pt(l))
  qc(l)=2.0*pt(l)/(pt(l-1)+2.0*pt(l)+pt(l+1))
  qr(l)=2.0*pt(l)/(pt(l)+pt(l+1))
enddo
ql(ki)=0.0
qc(ki)=0.0
qr(ki)=0.0
!compute normalized layer slopes
do l=1,ks
  call plmtrcs(ql,qc,qr,s(1,l),ss(1,l),ki)
enddo
return
end subroutine plmtrcx
