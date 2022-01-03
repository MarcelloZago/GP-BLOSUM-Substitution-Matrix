import pandas as pd
import numpy as np
import math

aa_dict = ["A", "R", "N"]

table = pd.DataFrame(np.zeros((len(aa_dict), len(aa_dict))), columns=aa_dict, index=aa_dict)

table['A']['A'] = 800
table['R']['A'] = 837
table['N']['A'] = 120
table['A']['R'] = 837
table['R']['R'] = 560
table['N']['R'] = 240
table['A']['N'] = 120
table['R']['N'] = 240
table['N']['N'] = 730

lamb = 0.347
numOfEntries = np.triu(table).sum()

table["sum"] = table.sum(axis=1)/numOfEntries
print(table)

for i in aa_dict:
    table[i] = table[i].apply(lambda x: x/numOfEntries)

for i in aa_dict:
    for j in aa_dict:
        table[i][j] = int(1/lamb * math.log((table[i][j]/(table['sum'][i]*table['sum'][j])) ,2))
table = table.drop(['sum'], axis=1)

print(table)