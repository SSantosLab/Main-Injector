import json
import astropy_healpix as ah
import numpy as np
from typing import Tuple
from abc import ABC, abstractmethod
from collections import namedtuple
from io import BytesIO
from astropy.table import Table
from base64 import b64decode
from io import BytesIO
from astropy.table import Table
from astropy.io import fits
from ligo.skymap.io.fits import read_sky_map, write_sky_map
from ligo.skymap.bayestar import rasterize

class Alert(ABC):
    @classmethod
    @abstractmethod
    def from_alert(self):
        pass

class GW(Alert):

    # Reference: 
    # https://stackoverflow.com/questions/1305532/how-to-convert-a-nested-python-dict-to-object
    def __init__(self, data: dict) -> None:
        for name, value in data.items():
            setattr(self, name, self._wrap(value))

    def __str__(self) -> str:
        return f'<class stages.gw.GW of event {self.superevent_id}.'
    
    @classmethod
    def from_alert(cls, alertpath: str) -> 'GW':
        instance = cls.__new__(cls)
        with open(alertpath, 'r') as f:
            alert = json.load(f)

        instance.__init__(alert)

        return instance

    def get_fields(self, prefix=''):
        """Extract all fields labels of the GW."""
        pass
    
        # attribute_names = []
        # attribute_values = []

        # for attr_name in dir(self):
        #     if not attr_name.startswith('_'):
        #         full_attr_name = f'{prefix}{attr_name}'
        #         attr_value = getattr(self, attr_name)

        #         if isinstance(attr_value, GW):
        #             # Recursively handle nested attributes
        #             nested_attributes = attr_value.get_fields(prefix=f'{full_attr_name}.')
        #             attribute_names.extend(nested_attributes._fields)
        #             attribute_values.extend(nested_attributes)

        #         else:
        #             attribute_names.append(full_attr_name)
        #             attribute_values.append(attr_value)

        # # Create a namedtuple class dynamically
        # all_attributes_namedtuple = namedtuple('AllAttributes', attribute_names)

        # # Create an instance of the namedtuple
        # attributes_namedtuple = all_attributes_namedtuple(*attribute_values)

        # return attributes_namedtuple
    
    def retrieve_skymaps(self,
                         moc_path: str,
                         flatten_path: str,
                         nside: int = 1024) -> None:
        """Retrives multiorder and flatten skymap from encoded ByteIO string.
            
        Parameters:
        -----------
        moc_path (str):
            Output path for retrieved multiorder resolution skymap.

        flatten_path (str):
            Output path for flatten HEALPix skymap.

        Returns:
        --------
            None.
        """

        if self.event.skymap is not None:

            skymap_bytes = b64decode(self.event.skymap)
            skymap = Table.read(BytesIO(skymap_bytes))
            skymap.write(moc_path, overwrite=True)

            self._flatten_skymap(moc_path, flatten_path, nside=nside)
        else:
            raise SkymapNotFound("Couldn\'t find encoded string for skymap.")
        
    def find_source_and_prob(self) -> Tuple[str,float]:
        """Returns the most probable Source and Prob."""

        sources = {
            self.event.classification.BNS: 'BNS',
            self.event.classification.NSBH: 'NSBH',
            self.event.classification.BBH: 'BBH',
            self.event.classification.Terrestrial: 'Terrestrial'
        }

        source_probs = [p[0] for p in sources.items()]
        event_prob = source_probs[np.argmax(source_probs)]
        source = sources[event_prob]

        return source, event_prob        
    
    def retrieve_coinc_map(self):
        pass

    def _wrap(self, value):
        if isinstance(value, (tuple, list, set, frozenset)): 
            return type(value)([self._wrap(v) for v in value])
        else:
            return GW(value) if isinstance(value, dict) else value
        
    def _flatten_skymap(self,
                        moc_skymap: str,
                        flatten_path: str,
                        nside: int = 1024) -> None:
        """Flatten multiorder skymap
        
        Parameters:
        -----------
            moc_skymap (str):
                Path to multiorder resolution skymap location.
            flatten_path
                Output path to flatten HEALPix skymap.
            
            Returns:
            --------
                None.
        """
        hdu = fits.open(moc_skymap)
        order = ah.nside_to_level(nside)
        table = read_sky_map(hdu, moc=True)
        table = rasterize(table, order=nside)
        write_sky_map(flatten_path, table, nest=True)
        hdu.close()

        return


class SkymapNotFound(Exception):
    pass
