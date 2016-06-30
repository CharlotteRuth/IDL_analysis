
dirs=['10','15','20','25','30']
ned=n_elements(dirs)
m=2.325e5
t=1e9
bin=1e8
per=fltarr(ned)
maxsfh=fltarr(ned)
avesf=fltarr(ned)
variation=fltarr(ned)

xmargin=[7.2,0.8]
ymargin=[3.0,0.2]
;!p.multi=[0,2,2]
paperps,file='wheresf.eps',charsize=1.5
;for j=1,ned-1 do begin
j=1
    ;files=file_search(dirs[j]+'kms/'+dirs[j]+'kms.01000')
    ;file=files[n_elements(files)-1]
    file=dirs[j]+'kms/'+dirs[j]+'kms.01000'
print,file
    rtipsy,file,h,g,d,s
    ntot=h.ngas+h.ndark+h.nstar
    rf=rbvector(file+'.rform')
    rf = {x:rf.x[h.ngas+h.ndark:ntot-1],y:rf.y[h.ngas+h.ndark:ntot-1],z:rf.z[h.ngas+h.ndark:ntot-1]}
    rsqr = s.x^2. + s.y^2. +s.z^2.
    rs = sqrt(rsqr)
    rfs = sqrt(rf.x^2. + rf.y^2. +rf.z^2.)
    ;plot,rfs,rs,/xlog,/ylog,ytit='Current Radius[kpc]',xtit='Formation Radius[kpc]',psym=1
    contour_plus,rfs,rs,ytit='Current Radius[kpc]', $
      xtit='Formation Radius[kpc]',threshold=20,xrange=[0.1,2.2], $
      xstyle=1,/ylog,yrange=[0.1,25],ystyle=1
    q=findgen(100)*3./100.
    oplot,q,q
ppsclose

; Information about the velocities of these particles when they formed
paperps,file='velsf.eps',charsize=2.0
    vf=rbvector(file+'.vform')
    ngd = h.ngas+h.ndark
    vf = {x:vf.x[ngd:ntot-1],y:vf.y[ngd:ntot-1],z:vf.z[ngd:ntot-1]}
    irs = where(rs GT 2)
    vfs = (vf.x*rf.x +vf.y*rf.y +vf.z*rf.z)/rfs

apo = apocenters(g,d,s)
; Find apocenters of stars  (B+T p. 107-8)
;   phis = [g.phi,d.phi,s.phi]
;   allrs= [sqrt(g.x^2. + g.y^2. +g.z^2.), $
;            sqrt(d.x^2. + d.y^2. +d.z^2.), $
;            sqrt(s.x^2. + s.y^2. +s.z^2.)]
;   nea = n_elements(allrs)-1
;   isr = sort(allrs)
;   apo = fltarr(n_elements(s))
;   for i =0,n_elements(s)-1 do begin
;    Larr =crossp([s[i].x,s[i].y,s[i].z],[s[i].vx,s[i].vy,s[i].vz])
;    Lsqr=transpose(Larr)#Larr
;    phi_effi = s[i].phi + Lsqr[0]/(2.0*rsqr[i])
;    goal = 0.5*vfs[i]^2. + phi_effi
;    temp = where(allrs[isr] EQ rs[i])
;    lowi = temp[0]
;    highi = nea
;    if(2*lowi LT highi) then highi = 2*lowi
;    loval = phis[isr[lowi]] + Lsqr[0]/(2.*allrs[isr[lowi]]^2.)
;    hival = phis[isr[highi]] + Lsqr[0]/(2.*allrs[isr[highi]]^2.)
;   addstep = highi - lowi
;    while(hival LT goal) do begin
;      if ( highi + addstep GE nea ) then begin
;        highi = nea - 1
;        break
;      endif else highi = highi + addstep ;$
;      hival = phis[isr[highi]] + Lsqr[0]/(2.*allrs[isr[highi]]^2.)
;    endwhile
;      ;if(2*highi LT nea) then highi = 2*highi $
;      ;else highi =nea
;    halfi = floor((highi-lowi)/2.)
;    newval = phis[isr[halfi]] + Lsqr[0]/(2.*allrs[isr[halfi]]^2.)
;    frdiff =abs((newval - goal)/goal)
;    its = 0
;    while (frdiff GT 1e-3) do begin
;     if ( hival LT loval) then begin ; should never do this part
;      if ( newval GT goal ) then lowi = halfi $
;      else highi = halfi
;     endif else begin
;      if ( newval GT goal ) then highi = halfi $
;      else lowi = halfi
;     endelse
;      halfi = floor((highi-lowi)/2)
;      halfi = lowi + halfi
;      newval = phis[isr[halfi]] + Lsqr[0]/(2.*allrs[isr[halfi]]^2.)
;      frdiff =abs((allrs[isr[halfi]] - allrs[isr[lowi]])/allrs[isr[halfi]])
;     if(its GT 17) then begin
;      if (highi - halfi LT 3) then continue
;      print,'Not converging'
;      STOP
;     endif
;     its = its+1
;    endwhile
;    apo[i] = allrs[isr[halfi]]
; Keplerian system only
    ;GM = -s[i].phi * rs[i]
    ;e = Lsqr[0]/GM
    ;a = Lsqr[0]/GM/(1.-e^2.)
    ;apo[i] = a*(1+e)
   ;endfor
    ;vfdotrfs = (vf.x[irs]*rf.x[irs] +vf.y[irs]*rf.y[irs] +vf.z[irs]*rf.z[irs])/rfs[irs]
    contour_plus,vfs,apo,ytit='Apocenter [kpc]', $
      xtit='Formation v!lr!n [km/s]',threshold=20,nlevels=6, $
      /ylog,yrange=[0.3,25], ystyle=1, xbinsize=0.5,ybinsize=0.5
    ;plot,rs[irs],vfdotrfs,psym=1,xtit='Current R [kpc]',ytit='r!lform!n * v!lform!n'
ppsclose

paperps,file='apowheresf.eps',charsize=1.5
    contour_plus,rfs,apo,ytit='Apocenter [kpc]', $
      xtit='Formation Radius [kpc]',threshold=20,xrange=[0.1,2.2], $
      xstyle=1,/ylog,yrange=[0.1,25],ystyle=1
    oplot,q,q
ppsclose
end
