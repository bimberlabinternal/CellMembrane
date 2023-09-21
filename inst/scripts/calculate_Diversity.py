import numpy as np
import pandas as pd
from tcrdist.diversity import generalized_simpsons_entropy

def calculate_Diversity(conga_clones_file, OutputFile, order1, order2):
  df = pd.read_csv(conga_clones_file, sep="\t")
  df = df.rename(columns={"va_gene": "v_a_gene", "vb_gene": "v_b_gene", "cdr3a": "cdr3_a_aa", "cdr3b": "cdr3_b_aa"})
  df["newcol"] = df.clone_id
  df[["libID","clonenum"]] = df.newcol.str.split("_", expand=True)
  libID = df.libID.unique()[0]
  order1 = int(order1)
  order2 = int(order2)
  temp = generalized_simpsons_entropy(df.loc[df['libID'] == libID, 'clone_size'], orders=np.arange(order1, order2))
  temp = temp.add_suffix("_" + libID)
  for libID in df.libID.unique()[1:]:
    div = generalized_simpsons_entropy(df.loc[df['libID'] == libID, 'clone_size'], orders=np.arange(order1, order2))
    div = div.add_suffix("_" + libID)
    temp = temp.join(div)
  temp.to_csv(OutputFile, na_rep="NA")
