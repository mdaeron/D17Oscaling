import numpy as np
import pandas as pd
import uncertainties as uc
from uncertainties import unumpy as unp

class DataSet():
	def __init__(self):
		pass

	def import_analyses(self, filepath, columns = {}):
		self.analyses = pd.read_csv(filepath).rename(columns = columns)
		self.compile()

	def import_analytes(self, filepath, columns = {}):
		self.analytes = pd.read_csv(
			filepath,
			converters = {
				"d18O": uc.ufloat_fromstr,
				"d18O_raw": uc.ufloat_fromstr,
				"D17O": uc.ufloat_fromstr,
				"D17O_raw": uc.ufloat_fromstr,
			},
		).rename(columns = columns).set_index('Analyte')

	def export_analyses(self, filepath):
		self.analyses.to_csv(filepath, index = False)

	def export_analytes(self, filepath):
		self.analytes.to_csv(filepath)

	def __repr__(self):
		v = vars(self)
		_ = 'DataSet('
		for a in v:
			_ += f'\n  {a} = {getattr(self, a)}'
		_ += '\n)'
		return _

	def compile(self):

		# compute pooled SD of d18O
		g = self.analyses.groupby('Analyte')['d18O_raw']
		chi2 = ((self.analyses['d18O_raw'] - g.transform('mean')) ** 2).sum()
		dof = len(self.analyses) - g.ngroups
		self.SD_d18O_raw = np.sqrt(chi2 / dof)

		# compute pooled SD of D17O
		g = self.analyses.groupby('Analyte')['D17O_raw']
		chi2 = ((self.analyses['D17O_raw'] - g.transform('mean')) ** 2).sum()
		dof = len(self.analyses) - g.ngroups
		self.SD_D17O_raw = np.sqrt(chi2 / dof)

		# create dataframe of analytes
		self.analytes = self.analyses.groupby('Analyte').agg(
			N = ('d18O_raw', 'count'),
			d18O_raw = ('d18O_raw', 'mean'),
			D17O_raw = ('D17O_raw', 'mean'),
		)

		# assign raw uncertainties
		for k, sd in (
			('d18O_raw', self.SD_d18O_raw),
			('D17O_raw', self.SD_D17O_raw),
		):
			self.analytes[k] = [
				uc.ufloat(v, sd / n**0.5)
				for v, n in zip(
					self.analytes[k],
					self.analytes['N'],
				)
			]

	def normalize(self, refmaterials, λ = 0.528):

		# populate self.refmaterials
		if hasattr(self, 'analytes'):
			loa = self.analytes.index.tolist()
		else:
			loa = self.analyses["Analyte"].unique().tolist()
		self.refmaterials = {k:v for k,v in refmaterials.items() if k in loa}

		# compute d17O of refmaterials
		for m in refmaterials.values():
			m['d17O'] = 1000 * (
				np.exp(m['D17O']/1000)
				* (1 + m['d18O']/1000) ** λ
				- 1
			)

		# compute d17O of analytes
		self.analytes['d17O_raw'] = 1000 * (
			unp.exp(self.analytes['D17O_raw']/1000)
			* (1 + self.analytes['d18O_raw']/1000) ** λ
			- 1
		)

		# compute affine transformations affecting d18O and d17O
		rm = self.analytes[self.analytes.index.isin(self.refmaterials)]
		for m in [18, 17]:

			# weighed regression of d{m}O_raw over d{m}O
			y = rm[f'd{m}O_raw'].to_numpy()
			yn, sy = unp.nominal_values(y), unp.std_devs(y)
			x = np.array([self.refmaterials[a][f'd{m}O'] for a in rm.index])

			X = np.column_stack([x, np.ones_like(x)])
			W = np.diag(sy ** -2)

			cov = np.linalg.inv(X.T @ W @ X)
			beta = cov @ X.T @ W @ yn

			# invert the affine transformation to go from d{m}O_raw to d{m}O
			inva, invb = uc.correlated_values(beta, cov)

			match m:
				case 18:
					self.a18 = 1/inva
					self.b18 = -invb/inva
				case 17:
					self.a17 = 1/inva
					self.b17 = -invb/inva

		# apply affine transformations to d18O and d17O
		self.analytes['d18O'] = self.a18 * self.analytes['d18O_raw'] + self.b18
		self.analytes['d17O'] = self.a17 * self.analytes['d17O_raw'] + self.b17

		# compute Δ17O
		self.analytes['D17O'] = 1000 * (
			unp.log(1 + self.analytes['d17O']/1000)
			- λ * unp.log(1 + self.analytes['d18O']/1000)
		)

if __name__ == '__main__':

	RM = {
		'VSMOW': {'d18O': 0., 'D17O': 0.},
		'SLAP': {'d18O': -55.5, 'D17O': 0.},
		# 'FOO': {'d18O': -12.5, 'D17O': -0.1},
	}
	data = DataSet()
	data.import_analyses('examples/data.csv', columns = dict(d18O = 'd18O_raw', D17O = 'D17O_raw'))
	data.normalize(RM)
	print(data)
	data.export_analyses('examples/analyses.csv')
	data.export_analytes('examples/analytes.csv')
	data = DataSet()
	data.import_analytes('examples/analytes.csv')
	data.normalize(RM)
	print(data)
