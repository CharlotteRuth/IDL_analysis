
dirs=['10','15','20','25','30']
ned=n_elements(dirs)
m=2.325e5
t=1e9
bin=1e8
per=fltarr(ned)
maxsfh=fltarr(ned)
avesf=fltarr(ned)
variation=fltarr(ned)
medb=fltarr(ned)

xmargin=[7.2,0.8]
ymargin=[3.0,0.2]
!p.multi=[0,2,2]
paperps,file='sfhs.eps',charsize=1.5
for j=1,ned-1 do begin
    ;files=file_search(dirs[j]+'kms/'+dirs[j]+'kms.01000')
    ;file=files[n_elements(files)-1]
    file=dirs[j]+'kms/'+dirs[j]+'kms.00500'
print,file
    rtipsy,file,h,g,d,s
    ;sfr,s,sfh=sfh,bin=5e7,maxt=maxt,massu=m,time=t,xmargin=[7.5,0.3],ymargin=[3.0,1.8],tit=dirs[j]+' km s!u-1!n'
    sfr,s,sfh=sfh,maxt=maxt,massu=m,time=t,tarray=tarray,xmargin=xmargin, $
       ymargin=ymargin
    legend,[dirs[j]+' km s!u-1!n'],/right
    nbins=maxt/1e9
    spb = 1e9/bin
    thing = fltarr(nbins)
    weights=fltarr(n_elements(tarray))+1.0
    a=[1e4,2e9]
    sandfit=curvefit(tarray,sfh,weights,a,sigma,function_name='sandagetau')
    for k=1,nbins-1 do begin
      ; 20 sfr bins per Gyr
      startbin=k*spb
      stopbin=(k+1)*spb
      temarray = sfh[startbin:stopbin] - sandfit[startbin:stopbin]
      thing[k] = max(temarray)-min(temarray)
    endfor
    maxsfh[j]=max(sfh)
    avesf[j]=mean(sfh)
    variation[j]=mean(thing)/avesf[j]
print,'Mean SFH',avesf[j]
print,'SFH Variation',variation[j]
    f=dofft(sfh)
    i=where(f eq max(f))
    per[j]=4.*maxt/i
print,'Period',per[j]
endfor
ppsclose

paperps,file='varvmass.eps',charsize=2
plot,dirs,variation,xtit='Virial Velocity [km/s]',ytit='Amplitude Variation',xmargin=[6.5,0.2],ymargin=[3.0,0.2],xrange=[5,35],psym=6
ppsclose

paperps,file='pervmass.eps',charsize=2
plot,dirs,per/1e6,xtit='Virial Velocity [km/s]',ytit='Period [Myr]',xmargin=[6.5,0.2],ymargin=[3.0,0.2],xrange=[5,35],psym=6
ppsclose

!p.multi=[0,2,2]
paperps,file='stuffvmass.eps',charsize=2
plot,dirs,per/1e6,xtit='Virial Velocity [km/s]',ytit='Period [Myr]',xmargin=[6.5,0.2],ymargin=[3.0,0.2],xrange=[5,35],psym=6

plot,dirs,maxsfh,xtit='Virial Velocity [km/s]',ytit='Max SFH',/ylog,xmargin=[7.0,0.4],ymargin=[3.0,0.2],xrange=[5,35],psym=6
plot,dirs,avesf,xtit='Virial Velocity [km/s]',ytit='<SFH>',/ylog,xmargin=[6.0,0.2],ymargin=[3.0,0.2],xrange=[5,35],psym=6
ppsclose


!p.multi=[0,2,2]
paperps,file='bs.eps',charsize=1.5
for j=0,ned-1 do begin
    ;files=file_search(dirs[j]+'kms/'+dirs[j]+'kms.[0-9]????')
    file=dirs[j]+'kms/'+dirs[j]+'kms.01000'
print,file
    rtipsy,file,h,g,d,s
    sfr,s,/dob,sfh=sfh,maxt=maxt,massu=m,time=t,xmargin=xmargin, $
       ymargin=ymargin,xtit=xtit,ytit=ytit,xtickname=xtickname, $
       ytickname=ytickname
     medb[j]=median(sfh)
print,"Median b:",medb[j]
legend,[dirs[j]+' km s!u-1!n'],/right
endfor
ppsclose

!p.multi=[0,2,2]
paperps,file='ages.eps',charsize=1.5
for j=1,ned-1 do begin
    file=dirs[j]+'kms/'+dirs[j]+'kms.01000'
print,file
    rtipsy,file,h,g,d,s
    sfr,s,sfh=sfh,maxt=maxt,massu=m,time=t,tarray=tarray,/xlog,/age, $
          xrange=[-0.5,9.5],xstyle=1
          ;xmargin=xmargin, ymargin=ymargin,
print,"Total Star Mass:",total(s.mass)*m
vxmom=moment(s.vx)
print,"X Velocity dispersion:",sqrt(vxmom[1])
endfor
ppsclose

END
;    if (j LT 3) then begin
;     ymargin=[0.5,0.2]
;     xtit=''
;     xtickname=replicate(' ',4)
;    endif else begin
;     ymargin=[3.0,0.2]
;     xtickname=replicate('',4)
;     xtit='Time [Gyr]'
;    endelse
;    ; right column
;    if (j EQ 2 OR j EQ 4) then begin
;     ytit = ''
;     xmargin=[4.4,0.4]
;    endif else begin
;     xmargin=[7.5,0.4]
;     ytit = 'SFR [M'+sunsymbol()+' yr!u-1!n]'
;     ytickname=replicate('',6)
;    endelse
;    ; top row
;    if (j LT 3) then begin
;     ymargin=[0.5,0.2]
;     xtit=''
;     xtickname=replicate(' ',4)
;    endif else begin
;     ymargin=[3.0,0.2]
;     xtickname=replicate('',4)
;     xtit='Time [Gyr]'
;    endelse
;    ; right column
;    if (j EQ 2 OR j EQ 4) then begin
;     ytit = ''
;     ytickname=replicate(' ',5)
;     xmargin=[0.4,0.4]
;    endif else begin
;     xmargin=[5.5,0.4]
;     ytit = 'b'
;     ytickname=replicate('',5)
;    endelse
