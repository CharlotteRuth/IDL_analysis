pro metsfrallhalo,tipsbase, MASSUNIT=massunit

if (keyword_set(maxrows) EQ 0) then maxrows=7
if (keyword_set(massunit) EQ 0) then massunit=3.171e15
set_plot,'ps'
device,file=tipsbase+'.halometsfr.ps',/landscape
  rtipsy,tipsbase,h,g,d,s
  rtipsy,tipsbase+'.gtp',gh,gg,gd,gs
  rtipsy,tipsbase+'.sogtp',soh,sog,sod,sos
  rdfloat,tipsbase+'.grp',gr,skip=1
  restore,'skidstattemp.sav'
  b=read_ascii(tipsbase+'.stat',template=temp)
  restore,'sovcirctemp.sav'
  so=read_ascii(tipsbase+'.sovcirc',template=temp)
  b.mtot=b.mtot*massunit
  b.mstar=b.mstar*massunit
  b.mgas=b.mgas*massunit
  istar = where(b.mstar GT 7e7)
  lmtot=min(alog10(b.mtot[istar]))
  hmtot=MAX(alog10(b.mtot[istar]))
  dmtot=hmtot-lmtot
  bg = where(so.mstar EQ max(so.mstar))
  ;r = so.rvir[bg]
r=1
  halos = where(((gs.x - sos[bg].x)^(2.0)+(gs.y - sos[bg].y)^(2.0)+(gs.z - sos[bg].z)^(2.0))^(0.5)  LE r) +1
;where(b.mtot EQ MAX(b.mtot))+1.))
  nh = n_elements(where(b.mstar[halos] GT 0.0))
  rows = ceil(nh^(0.5))
  if (rows GT maxrows) then rows = maxrows
  !p.multi=[0,rows,rows]

  minm = MIN(s.metals)
  maxm = MAX(s.metals)
  for i=0,n_elements(halos)-1 do begin
    ih=where(gr EQ halos[i])
    neg=ih-h.ndark-h.ngas
    negind=where(neg GT 0,num)
    if ((num GT 0) AND (n_elements(negind) GT 5)) then begin
      ind=neg[negind]
      sfrmet,s[ind],bint=1e8,MASS=massunit,minmet=minm,maxmet=maxm;,color=floor((alog10(b.mtot[i-1])-lmtot)*255/dmtot)
      legend,[strtrim(halos[i],2)],charsize=0.5,box=0,/right
    endif
  endfor
  ;colorbar,range=[lmtot,hmtot],title="log [halo DM mass (M!lsol!n)]",/vertical
device,/close
set_plot,'x'

end
