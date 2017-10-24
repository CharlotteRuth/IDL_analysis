function get_timestep,mass_unit,length_unit
;This function will take the mass_unit ( in solar masses)  and
;length_unit ( in kiloparsecs) of a simulation
;and return the time unit in years


return,3.872d3*SQRT((length_unit*3.0856805d21)^3/(mass_unit*1.98892d33))/31556926.0
;4.7158332e11*SQRT(length_unit^3/mass_unit)
end
