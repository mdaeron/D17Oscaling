from D17Oscaling import DataSet

RM = {
	'VSMOW': {'d18O': 0., 'D17O': 0.},
	'SLAP': {'d18O': -55.5, 'D17O': 0.},
}

data = DataSet()
data.import_analyses(
	'examples/data.csv',
	columns = {
		'd18O': 'd18O_raw',
		'D17O': 'D17O_raw'
	},
)

print(data)

data.normalize(RM)
print(data)

data = DataSet()
data.import_analytes('examples/analytes.csv')

print(data)

data.normalize(RM)
print(data)
