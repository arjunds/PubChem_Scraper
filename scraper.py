from pyteomics import mgf
import pubchempy as pcp
import numpy as np
import pandas as pd
import gc
from tqdm import tqdm

reader = mgf.MGF("ALL_GNPS.mgf")

output_csv = "active_assays.csv"

number_tracker = "index_num"

reader = list(enumerate(reader))
length = len(reader)

with open(number_tracker, "r") as f:
    curr_index = int(f.read())

pbar = tqdm(total=length, unit='compound', smoothing=0.1)

for spectrum_index, spectrum in reader:
    if spectrum_index <= curr_index:
        continue

    with open(number_tracker, 'w') as f:
        f.write(str(spectrum_index))

    pbar.update()

    smiles = spectrum['params']['smiles']
    inchi = spectrum['params']['inchi']
    name = spectrum['params']['name']
    spectrum_id = spectrum['params']['spectrumid']

    del spectrum
    gc.collect()

    using_inchi = False

    if not smiles or smiles == 'N/A':
        if not inchi or inchi == 'N/A':
           continue

        using_inchi = True

    compound = None
    try:
        if using_inchi:
            compound = pcp.get_compounds(inchi, namespace="inchi")
        else:
            compound = pcp.get_compounds(smiles, namespace="smiles")

    except:
        continue

    if compound is None:
        continue

    compound = compound[0]
    cid = compound.cid
    
    del compound
    gc.collect()

    assays = None
    try:
        assays = pd.read_csv("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + str(cid) + "/assaysummary/CSV")
    except:
        continue

    if assays is None:
        continue
    

    assays = assays[assays['Bioactivity Outcome'] == 'Active']
    assays["Name"] = np.repeat(name, len(assays))
    assays["Spectrum ID"] = np.repeat(spectrum_id, len(assays))
    
    if len(assays) > 0:
        assays.to_csv(output_csv, mode='a', header=False)
    
    del assays
    gc.collect()

pbar.close()
