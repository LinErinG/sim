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
sim->add_map, h04, intensity=0.001

; See what the map looks like so far...
test = sim->get( /total )
plot_map, test, grid=5

; Check out the ellipse parameters, for reference
test = sim->get( /ellipse )
help, test, /str

; Look at the individual maps, if desired.
test = sim->get( /maps )

