
FUNCTION ch_dem_gauss_compute_dem, ltemp, p, scl=scl, temp_log=temp_log

IF n_elements(scl) EQ 0 THEN scl=fltarr(3)+1.

t=10.^ltemp

;
; Compute DEM (phi). This expression comes from Eq. 3 of Warren &
; Brooks (2009, ApJ, 700, 762).
;
IF keyword_set(temp_log) THEN BEGIN
  phi=p[0] * scl[0] * exp( -(ltemp-p[1]*scl[1])^2/2./(p[2]*scl[2])^2 ) / (p[2]*scl[2]) / sqrt(2.*!pi) $
      / (t* alog(10.) )
ENDIF ELSE BEGIN
  phi=p[0] * scl[0] * exp( -(t-p[1]*scl[1])^2/2./(p[2]*scl[2])^2 ) / (p[2]*scl[2]) / sqrt(2.*!pi)
ENDELSE

return,phi

END
