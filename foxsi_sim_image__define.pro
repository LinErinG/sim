;;;; Name the object foxsi_sim_image

function foxsi_sim_image::init

	;-- allocate memory to pointers when initializing object
	self.source = ptr_new(/allocate)  
	self.energy = ptr_new(/allocate)  
	self.maps = ptr_new(/allocate)
	self.cubes = ptr_new(/allocate)
	self.total_map = ptr_new(/allocate)
		
	return,1

end

function foxsi_sim_image::get, energy=energy, ellipse=ellipse, maps=maps, $
	total_map=total_map, cubes=cubes, total_cube=total_cube
;;
;; This routine is the master GET function.
;; It will return all the data you might desire!
;;
;; ENERGY returns the energy array (which is either the default or user-defined).
;; ELLIPSE returns a source structure with all the ellipse parameters.
;; MAPS returns a map structure containing all the maps.
;; TOTAL_MAP returns a map structure with all the maps added.

	if keyword_set( energy ) then return, *(self.energy)
	if keyword_set( ellipse ) then return, *(self.source)
	if keyword_set( maps ) then return, *(self.maps)
	if keyword_set( cubes ) then return, *(self.cubes)
;;;	if keyword_set( total_map ) then return, *(self.total_map)		; old way
	if keyword_set( total_map ) then begin
		total = (*(self.maps))[0]
		total.data = total((*(self.maps)).data, 3)
		total.id = 'X-ray sources'
		return, total
	endif
	if keyword_set( total_cube ) then begin
		total = (*(self.cubes))[0].cube
		total.data = total((*(self.cubes)).cube.data, 4)
		for i=0, n_elements(total)-1 do total[i].id = strmid(total[i].id,strpos(total[i].id,'energy'))
		return, total
	endif
	return, -1		; if nothing has been returned yet, toss back an error

end

pro foxsi_sim_image::define_energy, bin=bin, emin=emin, emax=emax, array=array
;;
;; Define the energy array.  User set or defaults.
;; Meant to be parameterized, but can set specific array using ARRAY keyword.
;; Array can be 1D or 2D, but either way should represent bin edges.
;;

	; Default energy array parameters in case none is defined
	default, bin, 0.2
	default, emin, 3.
	default, emax, 100.

	nbin = fix((alog(emax)-alog(emin))/bin)
	energy1 = findgen(nbin)*bin + alog(emin)
	energy1 = exp(energy1)
	energy2 = get_edges( energy1, /edges_2)

	; If ARRAY is set, then specify the array itself.
	; Energy array can be 1D or 2D, but either way should represent bin edges.
	if keyword_set( ARRAY ) then begin
		case (size(array))[0] of
			1: begin
					energy1 = array
					energy2 = get_edges( energy1, /edges_2)
				 end
			2: begin
					energy2 = array
					energy1 = get_edges( energy2, /edges_1)
				 end
			else: begin
					print, 'Not an acceptable energy array.'
					return
				 end
		endcase
		nbin = n_elements( energy1 ) - 1
		emin = min(energy1)
		emax = max(energy1)
		bin = -1.
	endif

	energy_mid = get_edges( energy1, /mean )
	energy_wid = get_edges( energy1, /width )

	energy = {logbin:bin, min:emin, max:emax, n:nbin, energy1:energy1, $
				energy2:energy2, energy_mid:energy_mid, energy_wid:energy_wid}

	*(self.energy) = energy

	return

end

pro foxsi_sim_image::add_spectrum, spec=spec, name=name, vth=vth, bpow=bpow

	; SPEC should be an array of X-ray fluxes corresponding to the energy array
	; (energy array is either user-specified or default)
	; Alternatively, an array of spectral parameters can be defined via VTH or BPOW,
	; using the usual array of parameters for each.  (See f_vth.pro and f_bpow.pro.)
	; Both VTH and BPOW can be used at the same time if desired; the fluxes are added.

	if not keyword_set( SPEC ) and not keyword_set( VTH ) and not keyword_set( BPOW) then begin
		print, 'Need to input either a spectrum (keyword SPEC) or an '
		print, 'array of model parameters (keyword VTH or BPOW).'
		return
	endif
	
	if (keyword_set( VTH ) or keyword_set( BPOW )) and keyword_set( SPEC ) then begin
		print, 'Input a model parameter array OR a spectrum, but not both.'
		return
	endif

	if not keyword_set( NAME ) then begin
		print, 'No source name specified; returning.'
		return
	endif
	
	; If no energy array has been defined yet, then load the default.
	if not exist(*(self.energy)) then self.define_energy
	
	if not keyword_set( SPEC ) then spec = fltarr( n_elements((*(self.energy)).energy_mid) )

	if keyword_set( VTH )  then spec += f_vth(  (*(self.energy)).energy2, VTH )
	if keyword_set( BPOW ) then spec += f_bpow( (*(self.energy)).energy2, BPOW )

	; The spectrum needs to correspond to the energy array.
	; Check that spec and energy array are the same size.
	if n_elements( (*(self.energy)).energy_mid ) ne n_elements( SPEC ) then begin
		print, 'Spectrum and energy array are not the same size. Returning.'
		return
	endif
	
	; Create map cube for this source.
	; TO DO: HANDLE CASE WHERE >1 SOURCE FITS THE NAME
	index = where( (*(self.maps)).id eq name )
	if index[0] eq -1 then begin
		print, 'Source name not found.'
		return
	endif
	temp_cube = replicate( (*(self.maps))[index], n_elements((*(self.energy)).energy_mid) )
	for i=0, n_elements(temp_cube)-1 do begin
		temp_cube[i].data /= total(temp_cube[i].data)
		temp_cube[i].data *= spec[i]
		temp_cube[i].id += ' flux, energy '+string((*(self.energy)).energy1[i],format='(f5.1)') $
							+' -'+string((*(self.energy)).energy1[i+1],format='(f5.1)')+' keV
	endfor

	; Save the cube.  If cube structure doesn't exist yet, create it.
	if not isa( *(self.cubes), 'STRUCT') then begin
		*(self.cubes) = {name:name, cube:temp_cube}
	endif else *(self.cubes) = [(*(self.cubes)), {name:name, cube:temp_cube}]
	;stop

end

pro foxsi_sim_image::add_ellipse, intensity=intensity, xy_center=xy_center, $
								  fwhm=fwhm, rotation=rotation, name=name

	default, rotation, 0.
	default, intensity, 1.
	default, name, 'source'
	
	print, 'Adding elliptical Gaussian source of intensity ', intensity

	if n_elements(xy_center) lt 2 then begin
		print, 'ERROR: No center position defined.'
		return
	endif
	
	if fwhm[0] eq -1 then begin
		print, 'ERROR: No source size defined.'
		return
	endif

	if n_elements(fwhm) eq 1 then fwhm=[fwhm,fwhm]
	
	source = {SOURCE, name:name, intensity:intensity, xy_center:xy_center, fwhm:fwhm, $
			  rotation:rotation}

	if not isa( *(self.source), 'STRUCT') then begin
		*(self.source) = source
	endif else *(self.source) = [*(self.source), source]

	; Now we've got the structure.  Create a map and add it to the total map.
	; If overall map doesn't exist yet, also generate the default map.
	if exist( *(self.total_map) ) eq 0  then self->define_map
	ref_map = *(self.total_map)
	dim = (size(ref_map.data))[1:2]
	dx = ref_map.dx
	dy = ref_map.dy
	fov_cen = [ref_map.xc, ref_map.yc]
	xy = (XY_CENTER-fov_cen)/[dx,dy] + [dim[0],dim[1]]/2.
	source  = intensity*PSF_GAUSSIAN(npix = [dim[0], dim[1]], $ 
             /double, st_dev = fwhm/2.35482,        $
             centroid = xy, /norm )
	source_map = make_map(source, dx=dx, dy=dy, xc = fov_cen[0], yc = fov_cen[1], $
						  id=name )
	; Important!
	; In order to avoid sources trailing into infinity, cut this sources off at 
	; a "full-width-tenth-max."
	factor = 0.1
	source_map.data[ where( source_map.data lt factor*max(source_map.data) ) ] = 0.
	source_map = rot_map( source_map, rotation )
	source_map.roll_angle = 0.

	coreg = coreg_map( source_map, *(self.total_map), drot=0., /resc, /no_proj )
	;;;(*(self.total_map)).data += coreg.data
	
	; Save the individual map too.  If map structure doesn't exist yet, create it.
	if not isa( *(self.maps), 'STRUCT') then begin
		*(self.maps) = coreg
	endif else *(self.maps) = [*(self.maps), coreg]

	;;stop
	return

end

pro foxsi_sim_image::remove_source, name

	; Remove from ellipse list if it's there.
	index = where( (*(self.source)).name eq name )
	if index[0] ne -1 then remove, index, *(self.source)
	
	; Remove from maps if it's there.
	index = where( (*(self.maps)).id eq name )
	if index[0] ne -1 then remove, index, *(self.maps)
	
	; Regenerate total map without this source.
	(*(self.total_map)).data = 0.
	for i=0, n_elements( *(self.maps) ) - 1 do (*(self.total_map)).data += (*(self.maps))[i].data
	
	return
end


pro foxsi_sim_image::define_map, dim=dim, _extra=_extra
;;
;; Define the overall map for the simulation.
;;

	default, dim, 256

	if size(dim, /n_dim) eq 2 then data = fltarr(dim[0],dim[1]) $
		else data = fltarr(dim,dim)

	*(self.total_map) = make_map( data, _extra=_extra )

	return

end

pro foxsi_sim_image::add_map, map, intensity=intensity, name=name
;;
;; If an overall map has not yet been defined, this will initialize it with the input map.
;; If an overall map is already defined, this will coregister it to that map and 
;; add its data to the map.
;; INTENSITY is a normalization factor.
;;

	default, intensity, 1.
	default, name, 'source'
	
	temp = map
	temp.id = name
	temp.data *= intensity

	; Check to see if the map is already defined.
	if exist( *(self.total_map) ) then begin
		temp = coreg_map( temp, *(self.total_map), drot=0., /resc, /no_proj )
		(*(self.total_map)).data += temp.data
	endif else *(self.total_map) = temp
	
	; Save the individual map too.  If map structure doesn't exist yet, create it.
	if not isa( *(self.maps), 'STRUCT') then begin
		*(self.maps) = temp
	endif else begin
		copy = (*(self.maps))[0]
		copy.data = temp.data
		copy.id = name
		*(self.maps) = [*(self.maps), copy]
	endelse

	return
	
end

;function foxsi_sim_image::foxsi_sim_cube, dir=dir

	;;;; NOT YET FUNCTIONAL.  need to find a way to run the setup script.
	;;;; Or maybe just have the user run it first.

	; This runs the total cube through the FOXSI response.
	; FOXSI simulation software must be on the user's computer.
	; DIR is required and specifies the path to the FOXSI-SMEX IDL directory.
	; e.g. ~/foxsi/smex/code/foxsi-smex/idl
	
;	if not keyword_set( DIR ) then begin
;		print, 'FOXSI-SMEX IDL directory must be specified in keyword DIR.'
;		return, -1
;	endif
	
;	curdir=curdir()
;	cd, dir
;	.run foxsi-smex-setup-script
;	cd, curdir
;
;end

function foxsi_sim_image::cleanup

	;-- free memory allocated to pointer when destroying object
	ptr_free,self.ptr
	return, 1

end

pro foxsi_sim_image__define

	void={foxsi_sim_image, 	$
		source:ptr_new(), 	$
		energy:ptr_new(), 	$
		maps:ptr_new(),		$
		cubes:ptr_new(),	$
		total_map:ptr_new()	$
		}

 return 

end