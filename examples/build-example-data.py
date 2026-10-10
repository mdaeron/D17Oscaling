import numpy as np
import pandas as pd

rng = np.random.default_rng(seed=42)
sigmad18O = 0.10
sigmaD17O = 0.01


samples = {
	'FOO':   {'d18O': rng.uniform(-55.5, 0, 1).item(), 'D17O': rng.normal(0, 0.1, 1).item()},
	'BAR':   {'d18O': rng.uniform(-55.5, 0, 1).item(), 'D17O': rng.normal(0, 0.1, 1).item()},
	'BAZ':   {'d18O': rng.uniform(-55.5, 0, 1).item(), 'D17O': rng.normal(0, 0.1, 1).item()},
	'VSMOW': {'d18O':   0.0, 'D17O': 0.0},
	'SLAP':  {'d18O': -55.5, 'D17O': 0.0},
}
samplekeys = list(samples.keys())

n = 20
e18 = rng.normal(0, sigmad18O, n)
e17 = rng.normal(0, sigmaD17O, n)
data = []
for k in range(n):
	s = samplekeys[k % len(samples)]
	data.append({
		'UID': f'X{k+1:02.0f}',
		'Analyte': s,
		'd18O': samples[s]['d18O'] + e18[k],
		'D17O': samples[s]['D17O'] + e17[k],
	})

df = pd.DataFrame(data)
print(df)

df.to_csv("examples/data.csv", index=False)
