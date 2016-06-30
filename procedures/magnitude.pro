function mag,array
mag = 0
FOR i=0, N_ELEMENTS(array) -1 DO mag = mag + array[i]*array[i]
return,sqrt(mag)
end
