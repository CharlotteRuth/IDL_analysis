
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
paperps,file='rsfhs.eps',charsize=1.5
;for j=1,ned-1 do begin
j=1
    ;files=file_search(dirs[j]+'kms/'+dirs[j]+'kms.01000')
    ;file=files[n_elements(files)-1]
    file=dirs[j]+'kms/'+dirs[j]+'kms.01000'
print,file
    rtipsy,file,h,g,d,s
    rs = sqrt(s.x^2. + s.y^2. +s.z^2.)
    ;sfr,s,sfh=sfh,bin=5e7,maxt=maxt,massu=m,time=t,xmargin=[7.5,0.3],ymargin=[3.0,1.8],tit=dirs[j]+' km s!u-1!n'
    sr = s[where(rs GE 0. AND rs LT 1.)]
    cumsf,sr,/norm
    for k=1,4 do begin
      sr = s[where(rs GE k AND rs LT k+1.)]
      cumsf,sr,/norm,/over,line=k
    endfor
    legend,['0 < r < 1 kpc','1 < r < 2 kpc','2 < r < 3 kpc','3 < r < 4 kpc','r > 4 kpc'],line=indgen(5),/bottom,/right

ppsclose
end
