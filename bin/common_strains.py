    #!/usr/bin/python
import pandas as pd
import sys

pheno1 = pd.read_csv(sys.argv[1], sep='\t')
pheno2 = pd.read_csv(sys.argv[2], sep='\t')

common_strains = list(set(pheno1['strain']).intersection(set(pheno2['strain'])))

with open(sys.argv[3], 'w+') as f:
    for strain in common_strains:
        f.write(strain + '\n')