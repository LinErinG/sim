;;;; Name the object foxsi_sim_image

function foxsi_sim_image::init

	;-- allocate memory to pointers when initializing object
	self.source = ptr_new(/allocate)  
	self.energy = ptr_new(/allocate)  
	self.maps = ptr_new(/allocate)
	self.total_map = ptr_new(/allocate)
		
	return,1

end

function foxsi_sim_image::get, energy=energy, ellipse=ellipse, maps=maps, total_map=total_map
;;
;; This routine is the master GET function.
;; It will eventually return all the data you might desire!
;;
;; ENERGY returns the energy array (which is either the default or user-defined).
;; ELLIPSE returns a source structure with all the ellipse parameters.
;; MAPS returns a map structure containing all the maps.
;; TOTAL_MAP returns a map structure with all the maps added.

	if keyword_set( energy ) then return, *(self.energy)
	if keyword_set( ellipse ) then return, *(self.source)
	if keyword_set( maps ) then return, *(self.maps)
	if keyword_set( total_map ) then return, *(self.total_map)
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

	energy = {ENERGY, logbin:bin, min:emin, max:emax, n:nbin, energy1:energy1, $
				energy2:energy2, energy_mid:energy_mid}

	*(self.energy) = energy

	return

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

	coreg = coreg_map( source_map, *(self.total_map), drot=0., /resc, /no_proj )
	(*(self.total_map)).data += coreg.data
	
	; Save the individual map too.  If map structure doesn't exist yet, create it.
	if not isa( *(self.maps), 'STRUCT') then begin
		*(self.maps) = coreg
	endif else *(self.maps) = [*(self.maps), coreg]

	;;stop
	return

end

;pro foxsi_sim_image::remove_source, name
;
;	return
;end


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

pro foxsi_sim_image::add_map, map, intensity=intensity
;;
;; If an overall map has not yet been defined, this will initialize it with the input map.
;; If an overall map is already defined, this will coregister it to that map and 
;; add its data to the map.
;; INTENSITY is a normalization factor.
;;

	default, intensity, 1.

	; Check to see if the map is already defined.
	if exist( *(self.total_map) ) then begin
		coreg = coreg_map( map, *(self.total_map), drot=0., /resc, /no_proj )
		(*(self.total_map)).data += coreg.data*intensity
	endif else *(self.total_map) = map

	return
	
end

;;;;;;; NOT USED YET
pro foxsi_sim_image::plot, map

	; plot the input map
	if isa( map,'STRUCT' ) then plot_map, map
	return
	
end

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
		total_map:ptr_new()	$
		}

 return 

end