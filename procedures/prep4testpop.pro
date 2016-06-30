pro prep4testpop,file,outfile,timeunit=timeunit,massunit=massunit

if (keyword_set(timeunit) EQ 0) then timeunit = 1e9
if (keyword_set(massunit) EQ 0) then massunit = 2.325e5
rtipsy,file,h,g,d,s
sfr,s,massu=massunit,timeunit=timeunit,sfh=sfh
openw,lun,outfile,/get_lun
for i=0,h.nstar-1 do printf,lun,1,s[i].metals,alog10((h.time-s[i].tform)*timeunit),1

END
