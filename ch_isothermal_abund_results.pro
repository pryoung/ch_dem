

PRO ch_isothermal_abund_results, abstr, line_data, aa, sigmaa, ltemp, interr_scale=interr_scale

;+
; NAME:
;     CH_ISOTHERMAL_ABUND_RESULTS
;
; PURPOSE:
;     Updates the abundance data within the CH_DEM abundance structure
;     with the results from the DEM fitting. This routine is not
;     intended to be called directly by the user.
;
; CATEGORY:
;     CHIANTI; DEM; results.
;
; CALLING SEQUENCE:
;     CH_ISOTHERMAL_ABUND_RESULTS, Abstr, Line_Data, Aa, Sigmaa, Ltemp
;
; INPUTS:
;     Abstr:  An IDL structure in the format returned by
;             CH_DEM_PROCESS_ABUND.
;     Line_Data: An IDL structure in the format returned by
;                CH_DEM_READ_LINE_IDS.
;     Aa:     The fit parameters returned by the CH_DEM fitting
;             routine. 
;     Sigmaa: The 1-sigma errors on AA.
;     Ltemp:  The log temperature values at which the EM is calculated.
;             calculated. 
;
; OPTIONAL INPUTS:
;     Interr_Scale:  An additional fractional error to be added in
;                    quadrature to the fittng errors. The fraction is
;                    multiplied by the intensity.
;	
; OUTPUTS:
;     The "abund", "error" and "ratio" tags of ABSTR are updated based
;     on the fit results.
;
; EXAMPLES:
;     See the routines CH_DEM_ISOTHERMAL_FIT.
;
; MODIFICATION HISTORY:
;     Ver.1, 11-May-2025, Peter Young
;-



i=where(abstr.type EQ 1)
ab_ref=abstr[i[0]].abund

abstr_out=abstr

;
; 'aa' contains updated abundances for 'type 2' elements, so load
; these and their errors into abstr.
;
k=where(abstr.init_ab_ind NE -1,nk)
IF nk NE 0 THEN BEGIN
  abstr_out[k].abund=aa[abstr[k].init_ab_ind]*ab_ref
  abstr_out[k].error=sigmaa[abstr[k].init_ab_ind]*ab_ref
  abstr_out[k].ratio=abstr_out[k].abund/ab_ref
ENDIF

;
; 'type 0' elements are not part of the minimization procedure but
; their updated abundances are obtained algebraically from the DEM.
;
; Note that I have to use ld here (not ld_all or ld_fit)
;
k=where(abstr.type EQ 0,nk)
IF nk NE 0 THEN BEGIN
  FOR i=0,nk-1 DO BEGIN
    j=k[i]
    ii=where(line_data.element_num EQ abstr[j].elt_num)
    contrib=line_data[ii].contrib
    line_int=line_data[ii].int
    line_err=line_data[ii].err
    IF n_elements(interr_scale) NE 0 THEN line_err=sqrt(line_err^2 + (line_int*interr_scale)^2)

    ic=where(contrib NE 0.)
    y2=spl_init(ltemp[ic],alog10(contrib[ic]))
    yi=spl_interp(ltemp[ic],alog10(contrib[ic]),y2,aa[0])
    
    abstr_out[j].abund=line_int / 10.^yi / 10.^aa[1]
    abstr_out[j].error=abstr_out[j].abund * line_err / line_int
    abstr_out[j].ratio=abstr_out[j].abund/ab_ref
  ENDFOR 
ENDIF

abstr=temporary(abstr_out)

END
