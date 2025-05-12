

FUNCTION ch_dem_compute_ints, data

;+
; NAME:
;     CH_DEM_COMPUTE_INTS
;
; PURPOSE:
;     Computes a set of line intensities from a DEM function. See the
;     routines ch_dem_gauss_fit_fn and ch_dem_linear_fit_fn for examples
;     of how to use it.
;
; CATEGORY:
;     CHIANTI; differential emission measure; DEM; fitting.
;
; CALLING SEQUENCE:
;     Result = CH_DEM_COMPUTE_INTS( Data )
;
; INPUTS:
;     Data:  A structure containing various parameters related to the DEM.
;            See ch_dem_linear_fit_fn for an example of how to define it.
;            DATA is set to zero after the routine is called.
;
; OUTPUTS:
;     An array of intensities for all emission lines contained in
;     data.line_data.
;
; MODIFICATION HISTORY:
;     Ver.1, 12-May-2025, Peter Young
;-

;
; Unpack the data structure
;
phi=data.phi
line_data=data.line_data
ab_ref=data.ab_ref
abstr=data.abstr
ltemp=data.ltemp
p=data.p

dlogt=ltemp[1]-ltemp[0]
t=10.^ltemp

nl=n_elements(line_data)
int=dblarr(nl)

ab_type=abstr[line_data.ab_ind].type
init_ab_ind=abstr[line_data.ab_ind].init_ab_ind
abund=abstr[line_data.ab_ind].abund
FOR i=0,nl-1 DO BEGIN
  CASE ab_type[i] OF
    1: abfac=1.0 
    2: abfac=p[init_ab_ind[i]]
    3: abfac=abund[i]/ab_ref
  ENDCASE
  contrib=line_data[i].contrib
  int[i]=total(phi * abfac * ab_ref * contrib * dlogt * t* alog(10.) )
ENDFOR 

;
; Delete data since it's not needed any more
;
data=0.

return,int

END
