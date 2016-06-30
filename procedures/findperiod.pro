function findperiod,file,MASSUNIT=m,TIMEUNIT=t

IF (keyword_set(massunit) EQ 0) THEN m=2.325e5
IF (keyword_set (timeunit) EQ 0) THEN t=1e9
rtipsy,file,h,g,d,s

sfr,s,sfh=sfh,maxt=maxt,massu=m,timeu=t

f=dofft(sfh)
i=where(f eq max(f))

RETURN,4.*maxt/i
END
