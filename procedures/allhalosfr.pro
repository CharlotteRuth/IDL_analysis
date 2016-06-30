pro allhalosfr,tipsbase, MASSUNIT=massunit,KPCUNIT=kpcunit,MAXROWS=maxrows,ALLHALOS=allhalos,EARLY=early

!P.FONT=1
!P.CHARSIZE=1.9
!P.THICK=4
!p.CHARTHICK=4
!X.THICK=4
!Y.THICK=4
openw,lun,tipsbase+'.sfrtimes',/get_lun
if (keyword_set(maxrows) EQ 0) then maxrows=4
if (keyword_set(massunit) EQ 0) then massunit=2.73e12
if (keyword_set(kpcunit) EQ 0) then kpcunit=2.739e3
sysvel=73*kpcunit/2.894405e3
set_plot,'ps'
device,file=tipsbase+'.halosfr.ps',/landscape,SET_FONT='Times',/tt
  rtipsy,tipsbase,h,g,d,s
  rtipsy,tipsbase+'.gtp',gh,gg,gd,gs
  ;rtipsy,tipsbase+'.sogtp',soh,sog,sod,sos
  rdfloat,tipsbase+'.grp',gr,skip=1
;  mass = read_ascii_array(tipsbase+'.massform')
;  s.mass = mass[h.ngas+h.ndark:h.n-1]
  restore,'/home/stinson/bd/trace/amigastattemp.sav'
  b=read_ascii(tipsbase+'.stat',template=temp)
  ;restore,'/net/grads-1/stinson/trace/sovcirctemp.sav'
  ;so=read_ascii(tipsbase+'.sovcirc',template=temp)
  ;b.mtot=b.mtot*massunit
 ; b.mstar=b.mstar*massunit
 ; b.mgas=b.mgas*massunit
  b.mtot=b.mtot
  b.mstar=b.mstar
  b.mgas=b.mgas
 ; istar = where(b.mstar GT 7e7)
 ; lmtot=min(alog10(b.mtot[istar]))
 ; hmtot=MAX(alog10(b.mtot[istar]))
 ; dmtot=hmtot-lmtot
 ; if (keyword_set(early)) then starmasslim=1e-6 else starmasslim = 1e-5
  ;bg = where(so.rvir EQ max(so.rvir[where(so.mstar GT starmasslim)]))
  ;bg = 105
  ;r = so.rvir[bg]
  ;if (keyword_set(allhalos)) then r=1e20
r=1e20
  ;bgx=sos[bg].x
  ;bgy=sos[bg].y
  ;bgz=sos[bg].z
bgx=0
bgy=0
bgz=0
  chalos = where(((gs.x - bgx)^(2.0)+(gs.y - bgy)^(2.0)+(gs.z - bgz)^(2.0))^(0.5)  LE r,nhalos) +1
if (nhalos LT 1) then begin
  print,'No halos are within '+strtrim(string(r),2)+' (the virial radius of SO halo ID '+strtrim(string(bg),2)+')'
endif
;where(b.mtot EQ MAX(b.mtot))+1.))
  todo = chalos[where(b.mstar[chalos] GT 5e5,nh)]
if (nh GT 256) then nh = 256
  halos = todo[reverse(sort(b.mtot[todo]))]
  rows = floor(nh^(0.5))
  if (rows GT maxrows) then rows = maxrows
  !p.multi=[0,rows,rows]
  ;for i=0,n_elements(halos)-1 do begin
  for i=0,nh-1 do begin
    ih=where(gr EQ halos[i]+1)
    ;ih=where(gr EQ i+1)
    neg=ih-h.ndark-h.ngas
    negind=where(neg GT 0,num)
    if (num GT 10) then begin
      ind=neg[negind]
      hi = halos[i]+1
      ;hi = i+1
      ;himo = i
      himo = halos[i]
      sfrcosmo,s[ind],bin=5e7,xrange=[0,3],MASS=massunit;,MAXTIME=maxtime,HALFTIME=halftime,AVETIME=avetime,SANDTAU=sandtau;TOFFSET=toffset,color=floor((alog10(b.mtot[i-1])-lmtot)*255/dmtot)
      ;legend,[strtrim(hi,2)],charsize=0.5,box=0
      ;legend,[string(sysvel*b.maxvc[himo],format='(I)')+' km/s'],charsize=0.5*!P.CHARSIZE,box=0,/right,/bottom
      legend,[string(b.mstar[himo],format='(e11.2)')+' M'+sunsymbol(),string(halos[i]+1)],charsize=0.5*!P.CHARSIZE,box=0,/right
;      printf,lun,string(hi)+string(maxtime)+string(halftime)+string(avetime)+string(sandtau)+string(b.mtot[himo])+string(sysvel*b.maxvc[himo])+string(sysvel*b.outvc[himo])+string(b.mstar[himo])+string(b.mgas[himo])+string(((gs[himo].x-bgx)^(2.0)+(gs[himo].y-bgy)^(2.0)+(gs[himo].z-bgz)^(2.0))^(0.5))
    endif
  endfor
  ;colorbar,range=[lmtot,hmtot],title="log [halo DM mass (M!lsol!n)]",/vertical
device,/close
close,lun
set_plot,'x'
!p.multi=[0,1,1]

end
