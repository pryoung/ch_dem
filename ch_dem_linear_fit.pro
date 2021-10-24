
FUNCTION ch_dem_linear_fit, line_data, ldens=ldens, lpress=lpress, $
                            abund_file=abund_file, $
                            interr_scale=interr_scale, dir_lookup=dir_lookup, $
                            log=log, initial_params=initial_params,quiet=quiet, $
                            ab_elt_fix=ab_elt_fix, fixed_abund=fixed_abund, $
                            dlogt=dlogt, ltemp_nodes=ltemp_nodes, $
                            em=em, ioneq_file=ioneq_file, swtch_ab=swtch_ab

;+
; NAME:
;      CH_DEM_LINEAR_FIT
;
; PURPOSE:
;      Fits a multi-linear DEM to a set of observed line
;      intensities. The lines are defined in logT-logDEM space,
;      must be contiguous, and have fixed temperature nodes.
;
; CATEGORY:
;      CHIANTI; differential emission measure; DEM.
;
; CALLING SEQUENCE:
;      Result = CH_DEM_LINEAR_FIT(Line_Data)
;
; INPUTS:
;      Line_Data:   A structure containing details of the lines used
;                   to derive the DEM. Note that the model_int and
;                   logt_eff tags will be modified by
;                   ch_dem_gauss_int. 
;
; OPTIONAL INPUTS:
;      Ldens:  Specifies the logarithm of the electron number density
;              (units: cm^-3) to be used for the calculation. Default
;              is 9.0.
;      Abund_File:  The CHIANTI abundance file to be used for the
;                   calculation. The default is
;                   sun_coronal_1992_feldman_ext.abund.
;      Interr_Scale: A number between 0 and 1 that adds an additional
;                    amount to the lines' intensity errors that
;                    is interr_scale*intensity. For example,
;                    interr_scale=0.1 will add 10% of the intensity to
;                    the error (added in quadrature). This can be useful if
;                    the intensity errors are very small. 
;      Dir_Lookup:  Directory containing the pre-calculated lookup tables
;      Initial_Params:  A 3-element array specifying the initial parameters
;               for the fit. Only use this if the function fails to
;               return a solution. 
;      Quiet:   Suppresses some text produced by program.
;      Dlogt:   Specifies the logT intervals for the temperature
;               range. Default is 0.025
;      Logt_Min:  Specifies the minimum logT value for the temperature
;               range. Default is 5.0
;      Logt_Max:  Specifies the maximum logT value for the temperature
;               range. Default is 7.0
;      Ab_Elt_Fix:  Either an integer or a string specifying the
;                   reference element. For example, ab_elt_fix=26 or
;                   'fe' specifies iron. If not specified, then the
;                   routine uses the element with the most number of
;                   lines in LINE_DATA.
;      Ioneq_File:  The name of an ionization equilibrium file. The
;                   default is to use the CHIANTI file (!ioneq_file).
;      Swtch_Ab:  An array of same size as the number of unique
;                 elements in LINE_DATA. A zero indicates the
;                 element's abundance should be fixed, and a
;                 one indicates the abundance should be a
;                 variable. This is expected to be used only in
;                 special cases - see CHIANTI Technical Report
;                 No. 22. 
;
; KEYWORD PARAMETERS:
;      LOG:  If set, then the model Gaussian is defined in logT space
;            rather than T space.
;      QUIET: If set, text information to screen is not printed.
;      EM:   If set, then the linear segments are defined in
;            logT-logEM space rather than logT-logDEM.
;      FIXED_ABUND:  If set, then the abundances are not varied.
; 
; OUTPUTS:
;      A structure containing the results. The tags are:
;        .aa   3-element array containing fit params (EM0,T0,sigmaT).
;        .sigmaa Errors on fit params.
;        .chi2 Reduced chi^2 for fit.
;        .yfit The line intensities from the model.
;        .ltemp The log temperatures at which the DEM is defined.
;        .dem   The DEM function.
;        .lteff The effective temperatures for the lines.
;        .int_ratio  Ratios of observed to model intensities.
;        .int_ratio_err Errors on the ratios.
;        .log    This is keyword_set(log).
;        .coldepth  The column depth (cm) corresponding to the DEM. 
;        .abund  Structure containing information about the
;                abundances. 
;
; CALLS:
;      READ_ABUND, READ_IONEQ, 
;      CH_DEM_GOFNT, CH_DEM_PROCESS_INDEX_STRING, Z2ELEMENT,
;      CH_DEM_PROCESS_ABUND, CH_DEM_INIT_ABUND
;
; MODIFICATION HISTORY:
;      Ver.0.1, 17-Jul-2019, Peter Young.
;         Using ch_dem_gauss_fit (v.9) as a starting point.
;      Ver.0.2, 19-Jul-2019, Peter Young.
;         Fixed bug for ldens and lpress inputs; fixed bug for lpress
;         and /em implementations.
;      Ver.0.3, 23-Jul-2019, Peter Young
;         Now checks for $CH_LOOKUP_DEM; changed default location of
;         the aval file (pop_lookup_line_list.txt).
;      Ver.0.4, 23-Sep-2019, Peter Young
;         Now calls ch_dem_write_results.
;      Ver.0.5, 16-Jun-2021, Peter Young
;         /FIXED_ABUND is now implemented; the fit values are now
;         loaded into ld_fit; removed stop statement if fit
;         doesn't work; now checks to see if there are enough
;         free parameters. 
;      Ver.0.6, 17-Jun-2021, Peter Young
;         Added SWTCH_AB optional input; now exits gracefully if
;         ch_dem_process_abund fails. 
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> output=ch_dem_linear_fit(line_data [, ldens=, /log, initial_params=, abund_file='
  print,'                                   interr_scale=, dir_lookup=, dlogt=, logt_min= '
  print,'                                   logt_max=, /fixed_abund, initial_params=, /quiet,'
  print,'                                   swtch_ab= ])'
  return,-1
ENDIF 



; Check the ldens and lpress inputs. We need one or the other, but not
; both (or neither).
;
IF n_elements(ldens) NE 0 AND n_elements(lpress) NE 0 THEN BEGIN
  print,'% CH_DEM_LINEAR_FIT: Please specify either LDENS or LPRESS, but not both. Returning...'
  return,-1
ENDIF 
;
IF n_elements(ldens) EQ 0 AND n_elements(lpress) EQ 0 THEN BEGIN
  print,'% CH_DEM_LINEAR_FIT: Please specify either LDENS or LPRESS. Returning...'
  return,-1
ENDIF

nl=n_elements(line_data)
k=where(line_data.int NE -1.,nk)
IF nk EQ 0 THEN BEGIN
  print,'% CH_DEM_LINEAR_FIT: The input LINE_DATA structure does not contain observed intensities. Please run '
  print,'                     ch_dem_read_line_ids to load intensities from an input file. Returning...'
  return,-1
ENDIF 



IF n_elements(dlogt) EQ 0 THEN dlogt=0.05

;
; Handle the temperature nodes.
;
n_nodes=n_elements(ltemp_nodes) 
IF n_nodes EQ 0 THEN BEGIN
  ltemp_nodes=[5.0,6.0,7.0]
ENDIF

;
; Re-define the nodes if they don't line up with the temperature grid.
;
ltmp_nodes=round(ltemp_nodes/dlogt)*dlogt

;
; Now create the temperature array at which the contribution functions
; are defined.
;
ltemp=findgen(round((max(ltmp_nodes)-min(ltmp_nodes))/dlogt)+1)*dlogt + min(ltmp_nodes)
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
; INITIAL PARAMETER SETUP
; -----------------------
; The initial DEM values at the nodes are set to 10^22. However I used
; a scaling parameter (SCL) such that the actual initial values pass
; to the minimization routine are 1. If /EM has been set, then I need
; to adjust the scaling parameter due to the different definitions for
; DEM and EM.
;
scl=1d22
IF keyword_set(em) THEN scl=scl*10.^(median(ltmp_nodes))*alog(10.)*dlogt
IF n_elements(initial_params) EQ 0 THEN BEGIN
  init=make_array(n_elements(ltmp_nodes),/double,value=1.)
ENDIF ELSE BEGIN
  init=initial_params/scl
ENDELSE 



;
; ELEMENT ABUNDANCE SETUP
; -----------------------
; - Element abundance information is stored in ABSTR and AB_REF.
; - A subset of LD is extracted to LD_FIT -> these are the lines that
;   will be used in the DEM fitting procedure
;
IF n_elements(abund_file) EQ 0 THEN abund_file=!abund_file
chck=file_info(abund_file)
IF NOT chck.exists THEN BEGIN
  print,'% CH_DEM_LINEAR_FIT: could not find the element abundance file. Returning...'
  return,-1
ENDIF 
abstr=ch_dem_process_abund(ld,abund_file=abund_file,ab_elt_fix=ab_elt_fix, $
                           init=init,out_init=out_init,fixed_abund=fixed_abund, $
                           swtch_ab=swtch_ab)
IF n_tags(abstr) EQ 0 THEN return,-1
init=out_init
k=where(abstr.type EQ 1)
ab_ref=abstr[k[0]].abund
;
ab_type=abstr[ld.ab_ind].type
k=where(ab_type GE 1)
ld_fit=ld[k]
IF n_tags(ld_fit) EQ 0 THEN return,-1
nfit=n_elements(ld_fit)



;
; Modify the line intensity errors with INTERR_SCALE (if ncessary). 
;
int=ld_fit.int
err=ld_fit.err
;
; Add an additional component to the error if interr_scale is given.
;
IF n_elements(interr_scale) NE 0 THEN BEGIN
  IF interr_scale LT 0 OR interr_scale GT 1 THEN BEGIN
    print,'% CH_DEM_GAUSS_FIT:  the input interr_scale should take a value between 0 and 1. Returning...'
    return,-1
  ENDIF ELSE BEGIN
    err=sqrt(err^2 + (int*interr_scale)^2 )
  ENDELSE
ENDIF


;
; SET PARAMETER LIMITS (PARINFO structure)
; --------------------
;    - the DEM/EM values can not be zero (since I take the log later)
;      so I make the lower limit 10.^-10
;
np=n_elements(init)
parinfo=replicate({fixed: 0, limited: [0,0], limits:[0.d,0.d], step: 0.},np)
parinfo[*].limited=[1,0]
parinfo[0:n_nodes-1].limits[0]=1e-20



;
; I pass additional data to the fit function by using the structure
; 'other'. Note that 'line_data', 'abstr', 'ab_ref', 'ltemp' and 'dem'
; are essential. 
;
; The other tags are specific to the type of fit function.  
;
junk=temporary(other)
other={ line_data: ld_fit, $
        abstr: abstr, $
        ab_ref: ab_ref, $
        ltemp: ltemp, $
        scl: scl, $
        em: keyword_set(em), $
        n_nodes: n_nodes, $
        ltmp_nodes: ltmp_nodes}


;
; Check the number of free parameters to see if the fit is possible.
;
IF np GE nfit THEN BEGIN
   print,'% CH_DEM_LINEAR_FIT: there are too many free parameters! Try using /fixed_abund to'
   print,'                     fix the element abundances.'
   print,'  No. free parameters: '+trim(np)
   print,'  No. lines:           '+trim(nfit)
   print,''
   return,-1
ENDIF 

;
; This is the minimization procedure, producing fit parameters AA.
;
aa=mpfitfun('ch_dem_linear_fit_fn',findgen(nfit), int , err, init, $
            perror=sigmaa, /quiet, bestnorm=bestnorm,yfit=yfit, $
            parinfo=parinfo,status=status, functargs= other)

k=where(sigmaa EQ 0.,nk)
IF nk GT 0 THEN BEGIN
  print,'% CH_DEM_LINEAR_FIT: at least one of the parameter errors (sigmaa) is zero, therefore fit has failed.'
  print,'                    A possible solution is to adjust the initial parameters using the INIT optional input.'
  print,'                    The current values are:'
  print,format='(27x,6e12.3)',init[0:n_nodes-1]*scl
  return,-1
ENDIF 

;
; Need to re-scale the DEM fit parameters
;
aa[0:n_nodes-1]=aa[0:n_nodes-1]*scl
sigmaa[0:n_nodes-1]=sigmaa[0:n_nodes-1]*scl

;
; Compute the DEM from the fit parameters.
;
dem=ch_dem_linear_compute_dem(ltemp,ltmp_nodes,aa[0:n_nodes-1])

;
; This updates abstr with the optimized abundances.
;
ch_dem_abund_results,abstr,ld,aa,sigmaa,ltemp,dlogt,dem,interr_scale=interr_scale



chi2=bestnorm/float(nfit-np)

;
; Write out the following to the IDL window:
;  - the comparison between the observed intensity and model intensity
;    for each line
;  - the abundance results
;
ld_fit.model_int=yfit
IF NOT keyword_set(quiet) THEN ch_dem_write_results,ld_fit, abstr


IF NOT keyword_set(quiet) THEN BEGIN
  print,''
  print,'NODE VALUES: '
  print,'    Log Temp           Node            Err'
  FOR i=0,n_nodes-1 DO BEGIN
    print,format='(5x,f7.2,5x,e10.2," +/- ",e10.2)',ltmp_nodes[i],aa[i],sigmaa[i]
  ENDFOR 
  print,''
  print,format='("Reduced chi^2: ",f8.2)',chi2
ENDIF



;
; Compute the T_eff of each line in LD_ALL and populate the ab_ind
; and model_int tags of LD_ALL
;
lteff=fltarr(nl)
nab=n_elements(abstr)
FOR i=0,nl-1 DO BEGIN
  contrib_fn=ld_all[i].contrib
  func=contrib_fn*dem*10.^ltemp
  getmax=max(func,imax)
  lteff[i]=ltemp[imax]
  ld_all[i].logt_eff=lteff[i]
 ;
 ; Populate ld_all.ab_ind
 ;
  FOR j=0,nab-1 DO BEGIN
    k=where(ld_all.element_num EQ abstr[j].elt_num)
    ld_all[k].ab_ind=j
  ENDFOR
 ;
 ; Populate ld_all.model_int
 ;
  ab_i=abstr[ld_all[i].ab_ind].abund
  ld_all[i].model_int=ab_i*total(dem*contrib_fn*10.^ltemp*dlogt*alog(10.))
ENDFOR 


;
; Compute column depth distribution. 
;
IF n_elements(lpress) NE 0 THEN ldens=lpress - ltemp
h_to_e_ratio=proton_dens(ltemp,abund_file=abund_file,ioneq_file=ioneq_file,/hydrogen)
coldepth=10.^ltemp*alog(10.)*dem*dlogt/h_to_e_ratio/(10.^ldens)^2

;
; Create the output structure.
;
IF n_elements(ldens) EQ 0 THEN ldens=-1.
IF n_elements(lpress) EQ 0 THEN lpress=-1.
IF n_elements(interr_scale) EQ 0 THEN interr_scale=-1.
output={method: 'linear', $
        aa: aa, $
        sigmaa: sigmaa, $
        chi2: chi2, $
        yfit: yfit, $
        ltemp: ltemp, $
        dem: dem, $
        coldepth: coldepth, $
        line_data: ld_all, $
        abstr: abstr, $
        interr_scale: interr_scale, $
        log_dens: ldens, $
        log_press: lpress, $
        em: keyword_set(em) }

return,output

END
