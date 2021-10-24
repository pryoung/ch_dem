
FUNCTION ch_dem_init_abund, init, abstr, line_data=line_data

;+
; NAME:
;     CH_DEM_INIT_ABUND
;
; PURPOSE:
;     Takes the array of initial parameters used for the DEM
;     minimization process, and adds the initial element abundance
;     parameters. 
;
; CATEGORY:
;     CHIANTI; differential emission measure (DEM).
;
; CALLING SEQUENCE:
;     Result = CH_DEM_INIT_ABUND( Init, Abstr )
;
; INPUTS:
;     Init:   The array of initial parameters.
;     Abstr:  A structure in the format returned by
;             ch_dem_process_abund, containing the abundance
;             information for the elements used in the fit process.
;
; OPTIONAL INPUTS:
;     Line_Data: A structure in the format returned by
;                ch_dem_read_line_ids.pro. 
;	
; OUTPUTS:
;     The INIT array is returned, but with the additional abundance
;     parameters added.
;
; MODIFICATION HISTORY:
;     Ver.1, 18-Jul-2019, Peter Young
;-

IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> output=ch_dem_init_abund(init,abstr [ ,line_data= ])'
  return,init
ENDIF 

ninit=n_elements(init)

;
; Only "type 2" elements are added to the initial parameter array. If
; none, then just return the original array.
;
k=where(abstr.type EQ 2,nk)
IF nk EQ 0 THEN return,init

init_ab=[init,abstr[k].ratio]

;
; Add the initial abundance values to the LINE_DATA structure.
;
IF n_elements(line_data) NE 0 THEN BEGIN
  FOR i=0,nk-1 DO BEGIN
    j=where(line_data.element_num EQ abstr[k[i]].elt_num)
    line_data[j].init_ab_ind=ninit+i
  ENDFOR 
ENDIF 

return,init_ab

END
