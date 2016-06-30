function z_from_time_axes,axis,index,value
ageUniverse = 13.7346*1e9 ;wmap3_lookback(10000)
z = REVERSE(findgen(100)*10.0/100.0)
t = ageUniverse - wmap3_lookback(z)
redshift = spline(t,z,value*1e9)

;if (where(value*1e9 gt ageUniverse))[0] eq -1  THEN redshift[where(value*1e9 gt ageUniverse)] = 0
return,strtrim(ROUND(redshift),2)
end
