
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
; KEYWORD PARAMETERS:
;      NO_LINE_INTS: If set, then the points and errors for each line
;                    intensity will not be plotted.
;
; OUTPUTS:
;       Creates an IDL plot object showing the Gaussian DEM. For each
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
;-


w=window(dim=[600,450],background_color=bgcolor)

th=2
fs=14
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

p=plot(ltemp,dem,thick=th,xth=th,yth=th,/current, $
       ytitle='DEM / 10!u20!n K cm!u-5!n', $
       xtitle='Log!d10!n (Temperature / K)', $
       font_size=fs,xticklen=0.02,yticklen=0.02, $
       pos=[0.14,0.12,0.97,0.97],_extra=_extra,/xsty)


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
ENDIF

a=objarr(nab)

FOR i=0,nab-1 DO BEGIN
  z2element,input.abstr[i].elt_num,name
  a[i]=plot(p.xrange,p.yrange,symbol=symbol[i],sym_thick=th, $
            sym_size=sym_siz,hide=1,name=name,/overplot)
ENDFOR 

l=legend(target=a,sample_width=0,pos=[0.92,0.92],font_size=fs,thick=th, $
        horizontal_spacing=0.08)



IF input.method EQ 'gauss' THEN BEGIN 
xr=p.xrange
yr=p.yrange
;
IF n_elements(xpos_legend) EQ 0 THEN xpos_legend=0.95*xr[0]+0.05*xr[1]

emstr=trim(string(format='(f8.2)',alog10(input.aa[0])))
t1=text(xpos_legend,0.9*yr[1]+0.1*yr[0],/data,'Log EM!d0!n = '+emstr,font_size=fs)
IF input.temp_log EQ 0 THEN BEGIN 
  tstr='Log T!d0!n = '+trim(string(format='(f8.2)',alog10(input.aa[1])))
  sigstr='Log $\sigma_T$ = '+trim(string(format='(f8.2)',alog10(input.aa[2])))
ENDIF ELSE BEGIN
  tstr='Log T!d0!n = '+trim(string(format='(f8.2)',input.aa[1]))
  sigstr='$\sigma_T$ = '+trim(string(format='(f8.3)',input.aa[2]))
ENDELSE 
t2=text(xpos_legend,0.82*yr[1]+0.18*yr[0],/data,tstr,font_size=fs)
t3=text(xpos_legend,0.74*yr[1]+0.26*yr[0],/data,sigstr,font_size=fs)
ENDIF


return,w

END
