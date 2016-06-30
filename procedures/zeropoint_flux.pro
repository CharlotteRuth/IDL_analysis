
FUNCTION zeropoint_flux,filter,frequency = freqency
;Returns the zero_point flux for a filter (Vega magnitude system)
ergs_watt = 1e7 ;Watts/erg
m_A = 1e-10 ;m/Angstrom
m_mic = 1e-6 ;m/micrometer
m_cm = 1e-2 ;m/cm

IF keyword_set(frequency) THEN BEGIN
ENDIF ELSE BEGIN
;Bessell et al 1998
    IF strcmp(filter,'U_Johnson.res',6) THEN zero_point = 417.5*1e-11 ;ergs/s/Angstrom/cm^2
    IF strcmp(filter,'B_Johnson.res',6) THEN zero_point = 632*1e-11 ;ergs/s/Angstrom/cm^2
    IF strcmp(filter,'V_Johnson.res',6) THEN zero_point = 363.1*1e-11 ;ergs/s/Angstrom/cm^2
    IF strcmp(filter,'R_Cousins.res',6) THEN zero_point = 217.7*1e-11 ;ergs/s/Angstrom/cm^2
    IF strcmp(filter,'I_Cousins.res',6) THEN zero_point = 112.6*1e-11 ;ergs/s/Angstrom/cm^2
    IF strcmp(filter,'J_Johnson.res',6) THEN zero_point = 31.47*1e-11 ;ergs/s/Angstrom/cm^2
    IF strcmp(filter,'H_Johnson.res',6) THEN zero_point = 11.38*1e-11 ;ergs/s/Angstrom/cm^2
    IF strcmp(filter,'K_Johnson.res',6) THEN zero_point = 3.961*1e-11 ;ergs/s/Angstrom/cm^2
    IF strcmp(filter,'L_Johnson.res',6) THEN zero_point = 0.708*1e-11 ;ergs/s/Angstrom/cm^2
;Cohen+, 2008
    IF strcmp(filter,'J_2MASS.res',6)   THEN zero_point = 3.129e-13*ergs_watt/m_mic*m_A ;ergs/s/Angstrom/cm^2
    IF strcmp(filter,'H_2MASS.res',6)   THEN zero_point = 1.133e-13*ergs_watt/m_mic*m_A ;ergs/s/Angstrom/cm^2
    IF strcmp(filter,'K_2MASS.res',6)   THEN zero_point = 4.283e-14*ergs_watt/m_mic*m_A ;ergs/s/Angstrom/cm^2
ENDELSE

IF NOT keyword_set(zero_point) THEN  zero_point = 0
RETURN,zero_point
END
