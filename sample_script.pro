; This script shows a sample of how to use the FOXSI_SIM_IMAGE to simulate FOXSI 
; images for the SMEX and MIDEX proposals.

; Start the simulation object
sim = obj_new('foxsi_sim_image')

; Define a general map of all zeroes - this is the backdrop on which we'll draw.
; Any make_map keywords will work here.  Make sure to get the FOV you want.
; This can also be used to clear the map!

sim->define_map, dim=256, xc=950, yc=-220

; Define some sources.  Right now the "intensity" keyword is just a normalization factor.

; Coronal source:
sim->add_ellipse, xy_center=[995,-235], fwhm=[15.,15.], intensity=10, name='coronal'

; Footpoints:
sim->add_ellipse, xy_center=[923.,-199.], fwhm=[5.,5.], intensity=3, name='fp1'
sim->add_ellipse, xy_center=[907.,-265.], fwhm=[5.,5.], intensity=3, name='fp2'

; Add the thermal loop, this time using data from a map saved in a file.
restore,'fh1_figure_image_aia_rhsi.dat'	,/v		; thermal loop data from Säm
; Truncate the data to eliminate wings.
h04.data[ where( h04.data lt 0.18*max(h04.data) ) ] = 0.
sim->add_map, h04, intensity=0.001, name='loop'

; See what the map looks like so far...
test = sim->get( /total_map )
window, 0, xsi=500,ysi=500
!p.multi=1
plot_map, test, grid=5

; Check out the ellipse parameters, for reference
test = sim->get( /ellipse )
help, test, /str

; Look at the individual maps, if desired.
test = sim->get( /maps )
window, 1, xsize=600, ysize=600
!p.multi=[0,2,2]
for i=0, n_elements(test)-1 do plot_map, test[i], grid=5

; Try removing a source.
sim->remove_source, 'loop'
test = sim->get( /total_map )
window, 2, xsi=500,ysi=500
!p.multi=1
plot_map, test, grid=5

; Put it back.
sim->add_map, h04, intensity=0.001, name='loop'
test = sim->get( /total_map )
window, 2, xsi=500,ysi=500
!p.multi=1
plot_map, test, grid=5

; Add spectra for some sources.  (Example does identical spectra for all sources, but 
; the user can put in whatever they like.)
; General practice: MAPS must be added first, then spectral information is added
; to each by referring to name.  Name must be exact.

; First, an energy range should be defined.
; If this is not done, then the default energy range will be used.

; Add a thermal spectrum to the loop.
;sim->define_energy		; default energy bins
sim->define_energy, bin=0.05, emin=2., emax=20.	; special binning
energy = sim->get( /energy )
flux = f_vth( energy.energy2, [2.,1.0,1.] )

sim->add_spectrum, spec=flux, name='loop'
sim->add_spectrum, spec=flux, name='fp1'

; Get the individual source cubes to see what they look like.
; Each element of the CUBES array has a name and a cube for one source.
sources = sim->get( /cubes )
help, sources, /str
movie_map, sources[0].cube, /log

; Get the total cube that has all of the sources together.
; This is just a map array.
total = sim->get( /total_cube )
movie_map, total, /log

; Add the same spectrum to the other two sources and retrieve the total cube again.
sim->add_spectrum, spec=flux, name='fp2'
sim->add_spectrum, spec=flux, name='coronal'
total = sim->get( /total_cube )
movie_map, total, /log
