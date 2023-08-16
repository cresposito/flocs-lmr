from dataclasses import dataclass
import pandas as pd
import numpy as np

@dataclass
class GrainSizeProfile():
    station: str = None
    QBelleChasse: float = None
    source: str = None
    susp_filepath: str = None
    datetime: str = None
    lat_lon: tuple = None
    bed_depth_m: float = None
    depth_frac: float = None
    depth: float = None
    conc_sand: np.ndarray = None
    conc_mud: np.ndarray = None
    conc_total: np.ndarray = None
    meas: np.ndarray = None
    ADCP_filename: str = None,
    ADCP_ens_start: float = None,
    ADCP_ens_end: float = None,


    @classmethod
    def from_dict(cls, dict):
        cls = cls()
        cls.__dict__.update(dict)
        return cls