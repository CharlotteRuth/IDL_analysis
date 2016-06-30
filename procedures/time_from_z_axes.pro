function time_from_z_axes,axis,index,value
ageUniverse = 13.7346*1e9 ;wmap3_lookback(10000)
time = (ageUniverse - wmap3_lookback(value))/1e9
;time = (wmap3_lookback(value))/1e9
return,strtrim(ROUND(time),2)
end
