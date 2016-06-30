function sigsfrgas,file,tclip=tCLIP,MASSUNIT=massunit,TIMEUNIT=timeunit,_EXTRA=_extra
if(keyword_set(tclip) EQ 0) then tclip=0.1 ;max to clip at

IF (keyword_set (timeunit) EQ 0) THEN timeunit=1e9
IF (keyword_set(massunit) EQ 0) THEN massunit=2.362e5

addpath, '/astro/net/grads-1/stinson/trace'
rtipsy,file,h,g,d,s
rform=rbvector(file+'.rform')

xbin=max(rform.x)-min(rform.x)
ybin=max(rform.y)-min(rform.y)
mgas=better_hist2d(g.x,g.y,g.mass, binsize1=xbin, binsize2=ybin, $
                 obin1=xobin, obin2=yobin)

siggas = max(mgas)/(xbin*ybin)*massunit/1e6
print,'Size (kpc^2):'+string(xbin*ybin)

mform=rbarray(file+'.massform')
ind=where(s.tform GT h.time-tclip)
ind=ind+h.ngas+h.ndark
mstar=better_hist2d(rform.x[ind],rform.y[ind],mform[ind], binsize1=xbin, binsize2=ybin, $
                 obin1=xobin, obin2=yobin)
sigsfr = max(mstar)/(xbin*ybin)/tclip*massunit/timeunit
return,{sfr:sigsfr,gas:siggas}
END
