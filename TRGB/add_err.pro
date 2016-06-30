FUNCTION fit_obs_err,mag,err,fake_mag
;This function takes observed errors in magnitude and fits a 4th
;degree polynomial function to
;them so that the errors for the fake stars can be guessed at
sigmafit = poly_fit(mag,err,4)
fake_err = sigmafit[0] + sigmafit[1]*fake_mag + sigmafit[2]*fake_mag^2.0+sigmafit[3]*fake_mag^3 + sigmafit[4]*fake_mag^4 
return,fake_err
END

;add_err,'/astro/net/scratch2/christensen/TRGBDust/NGC0253-WIDE1/NGC0253-WIDE1.vi2.fits','/astro/net/scratch2/christensen/TRGBDust/Fake1/fake.vi','/astro/net/scratch2/christensen/TRGBDust/Fake1/fake.vi2'
PRO add_err,obsfile,fakefile,outfile
;This reads in an observational file and a file of fake magnitudes and
;outputs a file with the errors added to the fake magnitudes
infile=obsfile
columns = ['mag1_acs', 'mag1_err','mag2_acs','mag2_err']
data = mrdfits(infile,1,columns=columns)
readcol,fakefile,v,i
v_err = fit_obs_err(data.mag1_acs,data.mag1_err,v)
i_err = fit_obs_err(data.mag2_acs,data.mag2_err,i)
results = TRANSPOSE([[v],[v_err],[i],[i_err]])
OPENW,1,outfile
printf,1,results,FORMAT='(F,F,F,F)'
close,1
END
