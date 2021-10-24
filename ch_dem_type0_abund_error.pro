
FUNCTION ch_dem_type0_abund_error, dem_str, int, err, contrib, dlogt, dem_err=dem_err


;+
; NAME:
;     CH_DEM_TYPE0_ABUND_ERROR
;
; PURPOSE:
;     Computes the error bar for "type-0" abundance measurements by
;     propagating the errors from the DEM parameters.
;
; CATEGORY:
;     CHIANTI; differential emission measure (DEM); abundance.
;
; CALLING SEQUENCE:
;     Result = CH_DEM_TYPE0_ABUND_ERROR( Dem_Str )
;
; INPUTS:
;     Dem_Str:  A structure in the format returned by
;              ch_dem_linear_fit and ch_dem_gauss_fit.
;
; OPTIONAL INPUTS:
;     Int:    The intensity of the line for which the abundance error
;             is being derived.
;     Err:    The intensity error of the line for which the abundance error
;             is being derived.
;     Contrib: The contribution function of the line for which the
;              abundance error is being derived.
;     Dlogt:  The step size in log-T for the DEM calculation.
;
; OUTPUTS:
;     Returns the error on the abundance measurement for the line.
;
; OPTIONAL OUTPUTS:
;     Dem_Err: The 1-sigma errors on output.dem function. Note this
;              can be returned without doing the abundance error
;              calculation. 
;
; EXAMPLE:
;     Please check the routine ch_dem_linear_fit.pro to see how this
;     routine is called.
;
; MODIFICATION HISTORY:
;     Ver.1, 23-Jul-2019, Peter Young.
;-


;
; For the DEM function (in dem_str.dem) I need to derive a
; corresponding error array that is determined from the errors of the
; node points. This is an exercise in error propagation!
;

n_nodes=dem_str.n_nodes
ltemp_nodes=dem_str.ltemp_nodes
ltemp=dem_str.ltemp
dem=dem_str.dem
sigmaa=dem_str.sigmaa

n=n_elements(dem_str.dem)
dem_err=make_array(n,/double)

;
; In the code below y is the log of the DEM (phi), and we have:
;     alog10(phi) = y = y2*k + y1*(1-k)
; where k=(x-x1)/(x2-x1).
;
; Thus,  phi = (phi1)^(1-k) * (phi2)^k
;
; Propagating errors for these gives:
;
;      err_a = error( (phi1)^(1-k) ) = (phi1)^(1-k) * (1-k) * sig1 / phi1
;      err_b = error( (phi2)^k ) = (phi2)^k * k * sig2 / phi2
;
FOR i=0,n_nodes-2 DO BEGIN
  getmin=min(abs(ltemp_nodes[i]-ltemp),i0)
  getmin=min(abs(ltemp_nodes[i+1]-ltemp),i1)
  x=ltemp[i0:i1]
  x1=ltemp[i0]
  x2=ltemp[i1]
  y1=alog10(dem[i0])
  y2=alog10(dem[i1])
  k=(x-x1)/(x2-x1)
 ;
  a=dem[i0]^(1-k)
  b=dem[i1]^k
 ;
  err_a=(dem[i0])^(1-k) * (1-k) * sigmaa[i] / dem[i0]
  err_b=(dem[i1])^k * k * sigmaa[i+1] / dem[i1]
 ;
  i_temp=findgen(i1-i0+1) + i0
  dem_err[i_temp]= dem[i_temp] * $
                   sqrt( ( err_a/a )^2 + $
                         ( err_b/b )^2 )
ENDFOR 


;
; If INT wasn't input, then return. This is intended for when
; the user wants to access dem_err independently of the abundance
; calculation. 
;
IF n_elements(int) EQ 0 THEN return,-1

;
; Now compute errors for abundance.
;
; The abundance is given by int / den, where 'den' is the denominator
; given by sum( contrib*temp*dlogt*alog(10.)*dem ). The error on 'den'
; is derived from dem_err. Also need to account for the error on int.
;
a=contrib*10.^ltemp*dlogt*alog(10.)
den=total( a*dem )
err_den=sqrt( total( (a*dem_err)^2 ) )

abund_err=int/den * sqrt( (err/int)^2 + (err_den/den)^2 )

return,abund_err

END
