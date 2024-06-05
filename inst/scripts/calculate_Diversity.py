import numpy as np
import pandas as pd
from tcrdist.diversity import generalized_simpsons_entropy

def calculate_Diversity(conga_clones_file, outputFile, order1, order2):
  df = pd.read_csv(conga_clones_file, sep="\t")
  sampleId = df.sampleId.unique()[0]
  order1 = int(order1)
  order2 = int(order2)
  temp = generalized_simpsons_entropy(df.loc[df['sampleId'] == sampleId, 'clone_size'], orders=np.arange(order1, order2))
  temp = temp.add_suffix("_" + sampleId)
  for sampleId in df.sampleId.unique()[1:]:
    div = generalized_simpsons_entropy(df.loc[df['sampleId'] == sampleId, 'clone_size'], orders=np.arange(order1, order2))
    div = div.add_suffix("_" + sampleId)
    temp = temp.join(div)
  temp.to_csv(outputFile, na_rep="NA")
