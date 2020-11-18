import numpy as np
import pandas as pd

class Catalog ( object ):
    pass

class SpectralCatalog ( Catalog ):
    def __init__ ( self, fname, idname='CATAID' ):
        self.idname = idname
        self.catalog = self.read ( fname )
        self.source = fname
        
    def read ( self, fname ):
        cat = self._spec_read ( fname )
        cat = cat.set_index ( self.idname )
        cat = cat.loc[~cat.index.duplicated () ]
        return cat

class GAMACatalog ( SpectralCatalog ):
    def _spec_read ( self, fname ):
        return pd.read_csv ( fname )

class SyntheticPhotCatalog  ( Catalog ):
    def __init__ ( self, specdf, filter_set ):
        pass
