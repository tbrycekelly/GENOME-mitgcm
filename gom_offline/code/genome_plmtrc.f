subroutine plmtrc(si,pi,ki,ks, so,po,ko)
implicit none
c
integer ki,ks,ko
real    si(ki,ks),pi(ki+1),
&        so(ko,ks),po(ko+1),flag
c
c**********
c*
c  1) remap from one set of vertical cells to another.
c     method: piecewise linear across each input cell
c             the output is the average of the interpolation
c             profile across each output cell.
c
c  2) input arguments:
c       si    - scalar fields in pi-layer space
c       pi    - layer interface depths (non-negative m)
c                 pi(   1) is the surface
c                 pi(ki+1) is the bathymetry
c       ki    - 1st dimension of si     (number of  input layers)
c       ks    - 2nd dimension of si,so  (number of fields)
c       po    - target interface depths (non-negative m)
c                 po(k+1) >= po(k)
c       ko    - 1st dimension of so     (number of output layers)
c       flag  - data void (land) marker
c
c  3) output arguments:
c       so    - scalar fields in po-layer space
c
c  4) except at data voids, must have:
c           pi(   1) == zero (surface)
c           pi( l+1) >= pi(l)
c           pi(ki+1) == bathymetry
c           0 <= po(k) <= po(k+1)
c      output layers completely below the bathymetry inherit values
c      from the layer above.
c
c  5) Tim Campbell, Mississippi State University, October 2002.
C     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2005.
c*
c**********
c
real,parameter :: thin=1.e-6  !minimum layer thickness
c
integer i,k,l,lf
real    q,qc,zb,zc,zt,sok(ks)
real    sis(ki,ks),pit(ki+1)
c
c ---   compute PLM slopes for input layers
  do k=1,ki
    pit(k)=max(pi(k+1)-pi(k),thin)
  enddo
  call plmtrcx(pit,si,sis,ki,ks)
c ---   compute output layer averages
  lf=1
  zb=po(1)
  do k= 1,ko
    zt = zb
    zb = po(k+1)
*         WRITE(6,*) 'k,zt,zb = ',k,zt,zb
    if     (zb-zt.lt.thin .or. zt.ge.pi(ki+1)) then
c
c ---       thin or bottomed layer, values taken from layer above
c
      do i= 1,ks
        so(k,i) = so(k-1,i)
      enddo !i
    else
c
c           form layer averages.
c
      if     (pi(lf).gt.zt) then
        WRITE(6,*) 'bad lf = ',lf
        stop
      endif
      do i= 1,ks
        sok(i) = 0.0
      enddo !i
      do l= lf,ki
        if     (pi(l).gt.zb) then
*               WRITE(6,*) 'l,lf= ',l,lf,l-1
          lf = l-1
          exit
        elseif (pi(l).ge.zt .and. pi(l+1).le.zb) then
c
c               the input layer is completely inside the output layer
c
          q   = max(pi(l+1)-pi(l),thin)/(zb-zt)
          do i= 1,ks
            sok(i) = sok(i) + q*si(l,i)
          enddo !i
*               WRITE(6,*) 'L,q = ',l,q
        else
c
c               the input layer is partially inside the output layer
c               average of linear profile is its center value
c
          q   = max( min(pi(l+1),zb)-max(pi(l),zt), thin )/(zb-zt)
          zc  = 0.5*(min(pi(l+1),zb)+max(pi(l),zt))
          qc  = (zc-pi(l))/pit(l) - 0.5
          do i= 1,ks
            sok(i) = sok(i) + q*(si(l,i) + qc*sis(l,i))
          enddo !i
*               WRITE(6,*) 'l,q,qc = ',l,q,qc
        endif
      enddo !l
      do i= 1,ks
        so(k,i) = sok(i)
      enddo !i
    endif
  enddo !k
return
end subroutine plmtrc
