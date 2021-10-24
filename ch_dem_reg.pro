
FUNCTION ch_dem_reg, line_data, lpress=lpress, ldens=ldens, interr_scale=interr_scale


;+
; This calls the regularization method of Hannah & Kontar to obtain
; the DEM.
;
; PRY, 6-Aug-2019
;   I tried this for one data-set but it crashed when giving the /pos
;   keyword to data2dem_reg (it's meant to force a positive
;   solution). Note that it doesn't look like abundances can be
;   varied with this method.
;-

nl=n_elements(line_data)

;
; Set up temperatures for DEM
;
ltemp=findgen(21)/20+5.5
nt=n_elements(ltemp)


;
; Make a copy of LINE_DATA. From this point on, line_data is not
; used. 
;
ld_all=line_data

;
; Make sure ld_all has observed intensities.
;
nl=n_elements(ld_all)
k=where(ld_all.int NE -1.,nk)
IF nk EQ 0 THEN BEGIN
  print,'% CH_DEM_GAUSS_FIT: The input LINE_DATA structure does not contain observed intensities. Please run '
  print,'                     ch_dem_read_line_ids to load intensities from an input file. Returning...'
  return,-1
ENDIF



;
; Note on the LINE_DATA structures.
;   LINE_DATA  Input to this routine; not modified.
;   LD_ALL     Identical to LINE_DATA, but is modified by this routine
;              and then placed in the output structure (tag:
;              line_data).
;   LD         This is obtained from LD_ALL by using
;              ch_dem_process_blends, which sums any blends (if
;              necessary).
;   LD_FIT     The sub-set of LD that will actually be fitted. For
;              example, if there's only one line from an
;              element then this line is not used to compute
;              the DEM (it's only used to derive the abundance) 
;
; The following call modifies LD_ALL (by adding the "contrib" tag
; containing the contribution functions), and then sums any blends to
; create the LD structure.
;
ld=ch_dem_add_contrib(ld_all, ltemp, avalfile=avalfile, $
                      log_press=lpress, log_dens=ldens, $
                      ioneq_file=ioneq_file, dir_lookup=dir_lookup)
IF n_tags(ld) EQ 0 THEN BEGIN
  print,'% CH_DEM_LINEAR_FIT: error found. Returning...'
  return,-1
ENDIF


;
; Restrict to iron initially.
;
k=where(ld.element_num EQ 26,k)
ld=ld[k]


cfmatrix=ld.contrib
line_in=ld.int
eline_in=ld.err
;
; Add an additional component to the error if interr_scale is given.
;
IF n_elements(interr_scale) NE 0 THEN BEGIN
  IF interr_scale LT 0 OR interr_scale GT 1 THEN BEGIN
    print,'% CH_DEM_GAUSS_FIT:  the input interr_scale should take a value between 0 and 1. Returning...'
    return,-1
  ENDIF ELSE BEGIN
    eline_in=sqrt(eline_in^2 + (line_in*interr_scale)^2 )
  ENDELSE
ENDIF



; ;order of regularization, default is 0th
order=0
; ;control the regularization parameter/chisq of result in DEM space: reg_tweak=1
reg_tweak=1
; ;Use guess solution in final regularization? default is no, guess=0.
guess=1
;; Use the min of the EM loci curves as the initial guess solution
;; used to weight/create the constraint matrix and possibly in the regularization itself (if guess=1)
gloci=0


; run the regularization
reg=data2dem_reg(ltemp, CFmatrix, line_in, eline_in,$
 	mint=5.7, maxt=6.4, nt=50, $
	order=order,reg_tweak=reg_tweak, guess=guess, $
	channels=ld.label,gloci=gloci,/pos)
	
return,reg

END
