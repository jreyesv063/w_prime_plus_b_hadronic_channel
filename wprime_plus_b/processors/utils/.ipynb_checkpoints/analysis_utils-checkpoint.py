import os
import json
import numpy as np
import pandas as pd
import awkward as ak
import pyarrow as pa
import pyarrow.parquet as pq
from datetime import datetime
from typing import List, Union
from coffea.nanoevents.methods import candidate, vector


def normalize(var: ak.Array, cut: ak.Array = None) -> ak.Array:
    """
    normalize arrays after a cut or selection

    params:
    -------
    var:
        variable array
    cut:
        mask array to filter variable array
    """
    if var.ndim == 2:
        var = ak.firsts(var)
    if cut is None:
        ar = ak.to_numpy(ak.fill_none(var, np.nan))
        return ar
    else:
        ar = ak.to_numpy(ak.fill_none(var[cut], np.nan))
        return ar


def pad_val(
    arr: ak.Array,
    value: float,
    target: int = None,
    axis: int = 0,
    to_numpy: bool = False,
    clip: bool = True,
) -> Union[ak.Array, np.ndarray]:
    """
    pads awkward array up to ``target`` index along axis ``axis`` with value ``value``,
    optionally converts to numpy array
    """
    if target:
        ret = ak.fill_none(
            ak.pad_none(arr, target, axis=axis, clip=clip), value, axis=None
        )
    else:
        ret = ak.fill_none(arr, value, axis=None)
    return ret.to_numpy() if to_numpy else ret


def build_p4(cand: ak.Array) -> ak.Array:
    """
    builds a 4-vector

    params:
    -------
    cand:
        candidate array
    """
    return ak.zip(
        {
            "pt": cand.pt,
            "eta": cand.eta,
            "phi": cand.phi,
            "mass": cand.mass,
            "charge": cand.charge,
        },
        with_name="PtEtaPhiMCandidate",
        behavior=candidate.behavior,
    )


def save_dfs_parquet(fname: str, dfs_dict: dict) -> None:
    """
    save dataframes as parquet files
    """
    table = pa.Table.from_pandas(dfs_dict)
    if len(table) != 0:  # skip dataframes with empty entries
        pq.write_table(table, fname + ".parquet")


def ak_to_pandas(output_collection: dict) -> pd.DataFrame:
    """
    cast awkward array into a pandas dataframe
    """
    output = pd.DataFrame()
    for field in ak.fields(output_collection):
        output[field] = ak.to_numpy(output_collection[field])
    return output


def save_output(
    events: ak.Array,
    dataset: str,
    output: pd.DataFrame,
    year: str,
    channel: str,
    output_location: str,
    dir_name: str,
) -> None:
    """
    creates output folders and save dfs to parquet files
    """
    with open("wprime_plus_b/data/simplified_samples.json", "r") as f:
        simplified_samples = json.load(f)
    sample = simplified_samples[year][dataset]
    partition_key = events.behavior["__events_factory__"]._partition_key.replace(
        "/", "_"
    )
    date = datetime.today().strftime("%Y-%m-%d")

    # creating directories for each channel and sample
    if not os.path.exists(
        output_location + date + "/" + dir_name + "/" + year + "/" + channel
    ):
        os.makedirs(
            output_location + date + "/" + dir_name + "/" + year + "/" + channel
        )
    if not os.path.exists(
        output_location
        + date
        + "/"
        + dir_name
        + "/"
        + year
        + "/"
        + channel
        + "/"
        + sample
    ):
        os.makedirs(
            output_location
            + date
            + "/"
            + dir_name
            + "/"
            + year
            + "/"
            + channel
            + "/"
            + sample
        )
    fname = (
        output_location
        + date
        + "/"
        + dir_name
        + "/"
        + year
        + "/"
        + channel
        + "/"
        + sample
        + "/"
        + partition_key
    )
    save_dfs_parquet(fname, output)


def prod_unflatten(array: ak.Array, n: ak.Array):
    """
    Unflattens an array and takes the product through the axis 1

    Parameters:
    -----------
        array: array to unflat
        n: array with the number of objects per event. Used to perform the unflatten operation
    """
    return ak.prod(ak.unflatten(array, n), axis=1)


def delta_r_mask(first: ak.Array, second: ak.Array, threshold: float) -> ak.Array:
    """
    Select objects from 'first' which are at least threshold away from all objects in 'second'.
    The result is a mask (i.e., a boolean array) of the same shape as first.
    
    Parameters:
    -----------
    first: 
        objects which are required to be at least threshold away from all objects in second
    second: 
        objects which are all objects in first must be at leats threshold away from
    threshold: 
        minimum delta R between objects

    Return:
    -------
        boolean array of objects in objects1 which pass delta_R requirement
    """
    mval = first.metric_table(second)
    return ak.all(mval > threshold, axis=-1)


def chi2_test(topJet, wJet, top_sigma, w_sigma, top_mass_pdg,  w_mass_pdg):
        
    t = (topJet.mass - top_mass_pdg) / top_sigma
    w = (wJet.mass - w_mass_pdg) / w_sigma
    

    chi2 = t**2 + w**2
    

    return chi2


def topTagger_mask(first, second):
        
    mask = np.logical_or(first, second)
    

    return mask



def pdg_masses():
    
    with open("wprime_plus_b/jsons/wAndtop_masses.json", "r") as f:
        pdg = json.load(f)
                    
    top_mass_pdg = pdg['pdg']['top_mass']           
    w_mass_pdg = pdg['pdg']['w_mass'] 
    
    return top_mass_pdg, w_mass_pdg


def tagger_constants(self, type: str = "hadronic"):

    # W, top and chi2
    with open("wprime_plus_b/jsons/wAndtop_masses.json", "r") as f:
        masses = json.load(f)


    top_sigma = masses[type]['top_sigma']
    top_low_mass = masses[type]['top_low_mass']
    top_up_mass = masses[type]['top_up_mass']


    w_sigma = masses[type]['w_sigma']
    w_low_mass = masses[type]['w_low_mass']
    w_up_mass = masses[type]['w_up_mass']  

    chi2 = masses[type]['chi2']

    return top_sigma, top_low_mass, top_up_mass, w_sigma, w_low_mass, w_up_mass, chi2


def entries(input_array: ak.Array):
    """
    
    With this function we obtain the number of outputs of an akward array.
    
    """
    
    data_list = ak.to_list(input_array)
    
    # Values are filtered
    values_not_none = [value for value in data_list if value is not None]
    
    return print("Valores no None:", values_not_none) 


def true_sum(input_array: ak.Array):
    
    """
    The following function gives us the number of trues of a binary awkard array.
    """

    return ak.count_nonzero(ak.num(input_array) == 1)