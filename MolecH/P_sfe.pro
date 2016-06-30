FUNCTION P_sfe,P
;Relates Pressure to SFE, equation from Leroy '08
k = 5.25e-10
a = 1.2
c = 5.7902e-6

sfe = k*c*P^a/(c*P^a + 1)
RETURN,sfe
END
