import os
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy import units as u
from astropy import constants as co

gAB = 3631. * u.Jy

class FilterSet ( object ):
    def __init__ ( self ):
        self.filter_names = []

    def __getitem__(self,key):
        assert key in self.filter_names, f"The {key} filter is not loaded!"
        return getattr ( self, key )


class FilterResponse ( object ):
    def __init__ ( self, wavelength, transmission ):
        self.wavelength = wavelength
        self.transmission = transmission
        
    def getmag ( self, wv, flux, flux_error=None ):
        ''' Get apparent magnitude from flux-calibrated spectrum '''
        minwv = np.min ( self.wavelength[self.transmission > 0.05] )*u.AA
        maxwv = np.max ( self.wavelength[self.transmission > 0.05] )*u.AA
        if (minwv < wv.min()) or (maxwv > wv.max() ):
            if flux_error is not None:
                return np.NaN, np.NaN
            else:
                return np.NaN
            
        flux_cast = np.interp ( self.wavelength, wv[np.isfinite(flux)].to(u.AA).value,
                                flux[np.isfinite(flux)].to(u.erg/u.s/u.cm**2/u.AA).value)
        
        gAB_lambda = gAB * co.c / (self.wavelength*u.AA)**2
        gAB_lambda = gAB_lambda.to(u.erg/u.s/u.cm**2/u.AA).value

        dlambda = np.diff(self.wavelength)
        dlambda = np.concatenate([[dlambda[0]], dlambda])

        numerator = np.sum ( self.wavelength * flux_cast * self.transmission * dlambda )
        denom = np.sum ( self.wavelength * gAB_lambda * self.transmission * dlambda )

        if flux_error is not None:
            fluxerr_cast = np.interp ( self.wavelength, wv[np.isfinite(flux)].to(u.AA).value,
                                       flux_error[np.isfinite(flux)].to(u.erg/u.s/u.cm**2/u.AA).value)
            num_err = np.sqrt(np.sum ((self.wavelength * fluxerr_cast * self.transmission * dlambda)**2))

            mag_err = 2.5/(np.log(10)*numerator) * num_err

            return -2.5*np.log10(numerator/denom), mag_err
            
        return -2.5*np.log10(numerator/denom)


class SDSSFilters ( FilterSet ):
    def __init__ ( self, fitsname ):
        super ( SDSSFilters, self).__init__ ()
        self.load ( fitsname )
        
    def load ( self, fitsname ):
        sdss = fits.open( fitsname )

        for i in range (1,6):
            self.filter_names.append(sdss[i].name)
            filt = FilterResponse ( sdss[i].data['wavelength'], sdss[i].data['respt'] )
            # b.c. Merian sources are small (esp. for SDSS seeing/depth), using respt
            # which is response through 1.3 airmass for a point source
            setattr ( self, sdss[i].name, filt )
            
class MerianFilters ( FilterSet ):
    def load ( self, parname ):
        df = pd.read_csv(parname, skiprows=4, delim_whitespace=True, usecols=[1,2],
                         names='wavelength transmission'.split())
        filt = FilterResponse ( df.wavelength.values, df.transmission.values )
        name = os.path.basename(parname).split('.par')[0]

        self.filter_names.append(name)
        setattr ( self, name, filt )

class Spectrum ( object ):
    pass

class GAMASpectrum ( Spectrum ):
    def __init__ ( self, fname ):
        fobj = fits.open ( fname )
        hdr = fobj[0].header
        self.flux = fobj[0].data[0] * 1e-17 * u.erg/u.s/u.cm**2 /u.AA
        wv = hdr['CRVAL1'] + (np.arange(0, fobj[0].data.shape[1]) - hdr['CRPIX1'])*hdr['CD1_1']
        self.wavelength = wv * u.AA
        self.flux_error = fobj[0].data[1] * 1e-17 * u.erg/u.s/u.cm**2 /u.AA
