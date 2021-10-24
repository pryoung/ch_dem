
PRO ch_dem_write_results, ld_fit, abstr


;+
; NAME:
;     CH_DEM_WRITE_RESULTS
;
; PURPOSE:
;     Prints the results of the ch_dem optimization procedures to the
;     IDL input window. Usually called from within the CH_DEM fitting
;     routine (e.g., ch_dem_linear_fit). 
;
; CATEGORY:
;     CHIANTI; DEM; output.
;
; CALLING SEQUENCE:
;     CH_DEM_WRITE_RESULTS, Ld_Fit, Abstr
;
; INPUTS:
;     Ld_Fit:  A structure in the format returned by
;              ch_dem_read_line_ids but containing only those spectral
;              features that were used for the fit.
;     Abstr:   A structure in the format returned by
;              ch_dem_process_abund containing the abundance results
;              from the DEM method.
;
; OUTPUTS:
;     Prints text to the IDL input window giving the results from the
;     CH_DEM procedure. 
;
; MODIFICATION HISTORY:
;     Ver.1, 16-Jun-2021, Peter Young
;-


;
; Print out the comparison between observed and fitted
; intensities.
;
nfit=n_elements(ld_fit)
print,'    Ion                 Wvl   Obs_Int   Obs_Err   Model_Int'
FOR i=0,nfit-1 DO BEGIN
  print,format='(a7,a20,2f10.2,f12.2)', ld_fit[i].ion,ld_fit[i].label, $
        ld_fit[i].int,ld_fit[i].err,ld_fit[i].model_int
ENDFOR

;
; Now print abundance results
;
nab=n_elements(abstr)
print,''
print,' Element   Abundance (x10^6)  Log Ab     Type'
FOR i=0,nab-1 DO BEGIN
  z2element,abstr[i].elt_num,name,/symbol
  CASE abstr[i].type OF
    0: type='algebraic'
    1: type='reference'
    2: type='optimized'
    3: type='fixed'
  ENDCASE
  IF abstr[i].type EQ 2 OR abstr[i].type EQ 0 THEN BEGIN
    errstr=' +/-'+string(format='(f6.2)',abstr[i].error*1e6)
  ENDIF ELSE BEGIN
    errstr='          '
  ENDELSE 
  print,format='(6x,a2,f10.2,a10,f8.2,a12)',name,abstr[i].abund*1e6,errstr, $
        alog10(abstr[i].abund)+12.,type
ENDFOR


END
