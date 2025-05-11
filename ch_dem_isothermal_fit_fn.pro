
FUNCTION ch_dem_isothermal_fit_fn, x, p, _extra=e


;+
; NAME:
;      CH_DEM_ISOTHERMAL_FIT_FN
;
; PURPOSE:
;      A fit function in the format used by the MPFIT procedures. To
;      be used with ch_dem_isothermal_fit.pro.
;
; CATEGORY:
;      CHIANTI; differential emission measure; DEM; fitting.
;
; CALLING SEQUENCE:
;      Result = CH_DEM_ISOTHERMAL_FIT_FN( X, P )
;
; INPUTS:
;      X:    Should be a simple index array that identifies the lines
;            to be computed. 
;      P:    A 1D array containing the fit parameters. Element 0 is
;            log T, and element 1 is log EM.
;      _EXTRA: This is expected to be a structure containing auxiliary
;              information needed to compute the intensities. It is
;              defined in ch_dem_isothermal_fit. 
;
; OUTPUTS:
;      Computes line intensities (units: erg/cm2/s/sr) for a set of
;      N(X) emission lines. 
;
; MODIFICATION HISTORY:
;      Ver.1, 11-May-2025, Peter Young
;-



line_data=e.line_data
abstr=e.abstr
ltemp=e.ltemp
ab_ref=e.ab_ref

n=n_elements(line_data)
int=dblarr(n)

ab_type=abstr[line_data.ab_ind].type
init_ab_ind=abstr[line_data.ab_ind].init_ab_ind
abund=abstr[line_data.ab_ind].abund
FOR i=0,n-1 DO BEGIN
  CASE ab_type[i] OF
    1: abfac=1.0 
    2: abfac=p[init_ab_ind[i]]
    3: abfac=abund[i]/ab_ref
  ENDCASE
  contrib=line_data[i].contrib
  k=where(contrib NE 0.)
  
  y2=spl_init(ltemp[k],alog10(contrib[k]))
  yi=spl_interp(ltemp[k],alog10(contrib[k]),y2,p[0])
  int[i]=abfac*ab_ref*10.^yi*10.^p[1]

ENDFOR 


return,int

END

