import numpy as np
import pandas
import torch
import rinalmo.pretrained


direct_repeat = 'AACACCGTAATTTCTACTCTTGTAGAT'

# Read Kim et al. data
df = pandas.read_excel('https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.4104/MediaObjects/41592_2017_BFnmeth4104_MOESM79_ESM.xlsx', header=1)

grnas = df["Guide sequence (5' to 3')"].apply(lambda x: direct_repeat + x)

# Compute features
model, alphabet = rinalmo.pretrained.get_pretrained_model(model_name='giga-v1')
model = model.to(device='cuda')
model.eval()

tokens = torch.tensor(alphabet.batch_tokenize(grnas.values), dtype=torch.int64, device='cuda')
with torch.no_grad(), torch.cuda.amp.autocast():
  outputs = model(tokens)

# Features shape: samples x nucleotides x embedding dimension
np.save('rinalmo-giga-v1_features.npy', outputs['representation'].cpu().numpy())