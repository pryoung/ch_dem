
FUNCTION ch_dem_plot, input, no_line_ints=no_line_ints, _extra=_extra, $
                      xpos_legend=xpos_legend

;+
; NAME:
;      CH_DEM_PLOT
;
; PURPOSE:
;      Takes output from ch_dem_gauss_fit and plots the final DEM
;
; CATEGORY:
;      CHIANTI; differential emission measure; DEM; plot.
;
; CALLING SEQUENCE:
;      Result = CH_DEM_PLOT( Input )
;
; INPUTS:
;      Input:  The structure returned by CH_DEM_GAUSS_FIT.
;
; OPTIONAL INPUTS:
;      Xpos_Legend: Allows the x-position of the legend to be
;                   adjusted.
;
; KEYWORD PARAMETERS:
;      NO_LINE_INTS: If set, then the points and errors for each line
;                    intensity will not be plotted.
;
; OUTPUTS:
;       Creates an IDL plot object showing the DEM function obtained
;       from the ch_dem routines. For each
;       emission line a point with error bars is shown. If the model
;       perfectly reproduces the line intensity, then the point will
;       lie on the DEM curve. If the model intensity is too high, then
;       the point will appear below the curve, indicating that the DEM
;       would need to be lower to reproduce the line's intensity.
;
; MODIFICATION HISTORY:
;       Ver.1, 20-Jul-2018, Peter Young.
;       Ver.2, 23-Jul-2018, Peter Young.
;         Now handles the logT case correctly.
;       Ver.3, 17-Jul-2019, Peter Young.
;         DEM now interpolated on linear scale rather than log; added
;         xpos_legend; increased symbol size.
;       Ver.4, 07-Nov-2024, Peter Young.
;         Modified how MCMC DEM is plotted; fixed bug when
;         /no_line_ints was set.
;-


;
; Set some default parameters for the plot.
;
IF n_elements(dimension) EQ 0 THEN dimension=[600,450]
IF n_elements(position) EQ 0 THEN position=[0.14,0.12,0.97,0.97]
IF n_elements(thick) EQ 0 THEN th=2 ELSE th=thick
IF n_elements(font_size) EQ 0 THEN font_size=12 ELSE fs=font_size
sym_siz=1.5

scl=1d20

ltemp=input.ltemp
dem=input.dem/scl

int_ratio=input.line_data.int/input.line_data.model_int

int=input.line_data.int
model_int=input.line_data.model_int
err=input.line_data.err
IF input.interr_scale NE -1. THEN err=sqrt(err^2 + (int*input.interr_scale)^2)
int_ratio=int/model_int
int_ratio_err=err/model_int

;
; The DEM for the MCMC method is typically plotted in a "stairstep" style.
;
IF input.method EQ 'mcmc' THEN stairstep=1


w=window(dimension=dimension,background_color=bgcolor)


p=plot(ltemp,dem,thick=th,xth=th,yth=th,/current, $
       ytitle='DEM / 10!u20!n K cm!u-5!n', $
       xtitle='Log!d10!n (Temperature / K)', $
       xticklen=0.02,yticklen=0.02, $
       _extra=_extra,/xsty, $
       stairstep=stairstep)

;
; Plots the DEM uncertainties for the MCMC method.
;
IF input.method EQ 'mcmc' AND keyword_set(no_line_ints) THEN BEGIN
  nt=n_elements(ltemp)
  demerr=input.demerr/scl
  FOR i=0,nt-1 DO BEGIN
    px=plot(/overplot,ltemp[i]*[1,1],demerr[i,*],th=th)
  ENDFOR 
ENDIF 


;
; Overplot the points and error bars for each emission line. 
;
IF NOT keyword_set(no_line_ints) THEN BEGIN
  symbol_all=['X','o','s','D','tu','tl','tr','p']
  nab=n_elements(input.abstr)
  symbol=symbol_all[0:nab-1]
 ;
  lteff=input.line_data.logt_eff
  nl=n_elements(lteff)
 ;
  yy=fltarr(nl)
  ee=fltarr(nl)
  FOR i=0,nl-1 DO BEGIN
    getmin=min(abs(ltemp-lteff[i]),imin)
    yy[i]=dem[imin]*int_ratio[i]
    ee[i]=dem[imin]*int_ratio_err[i]
  ENDFOR 
 ;
  q=objarr(nab)
  FOR i=0,nab-1 DO BEGIN
    k=where(input.line_data.ab_ind EQ i,nk)
    IF nk NE 0 THEN BEGIN
      q[i]=errorplot(lteff[k],yy[k],ee[k],linestyle='none',thick=th,/overplot, $
                  symbol=symbol[i],errorbar_thick=th,sym_thick=th, $
                  sym_size=sym_siz)
    ENDIF 
  ENDFOR 

  a=objarr(nab)
  
  FOR i=0,nab-1 DO BEGIN
    z2element,input.abstr[i].elt_num,name
    a[i]=plot(p.xrange,p.yrange,symbol=symbol[i],sym_thick=th, $
              sym_size=sym_siz,hide=1,name=name,/overplot)
  ENDFOR 

  l=legend(target=a,sample_width=0,pos=[0.92,0.92],font_size=font_size,thick=th, $
           horizontal_spacing=0.08)

ENDIF

;
; The following adds text giving the Gaussian fit parameters for the 'gauss'
; method.
;
IF input.method EQ 'gauss' THEN BEGIN 
  xr=p.xrange
  yr=p.yrange
;
  IF n_elements(xpos_legend) EQ 0 THEN xpos_legend=0.95*xr[0]+0.05*xr[1]

  emstr=trim(string(format='(f8.2)',alog10(input.aa[0])))
  t1=text(xpos_legend,0.9*yr[1]+0.1*yr[0],/data,'Log EM!d0!n = '+emstr,font_size=font_size)
  IF input.temp_log EQ 0 THEN BEGIN 
    tstr='Log T!d0!n = '+trim(string(format='(f8.2)',alog10(input.aa[1])))
    sigstr='Log $\sigma_T$ = '+trim(string(format='(f8.2)',alog10(input.aa[2])))
  ENDIF ELSE BEGIN
    tstr='Log T!d0!n = '+trim(string(format='(f8.2)',input.aa[1]))
    sigstr='$\sigma_T$ = '+trim(string(format='(f8.3)',input.aa[2]))
  ENDELSE 
  t2=text(xpos_legend,0.82*yr[1]+0.18*yr[0],/data,tstr,font_size=font_size)
  t3=text(xpos_legend,0.74*yr[1]+0.26*yr[0],/data,sigstr,font_size=font_size)
ENDIF


return,w

END
