"""
Perform two-point normalization of triple oxygen isotope measurements
based on two or more reference materials with sufficiently different δ18O values.

The normalization applies two affine transformations:

  * δ18O_new = a_18 ∙ δ18O_old + b_18
  * δ17O_new = a_17 ∙ δ17O_old + b_17

Provided with:
  (1) a set of measurements of one or more unknown samples, with uncertainties;
  (2) a set of measurement of two or more known reference materials, with uncertaines;
  (3) nominal (δ18O, Δ17O) values for the reference materials,
`D17Oscaling` computes the normalized (δ18O, Δ17O) values for all unknown samples
and the variance-covariance matrix describing the (correlated) uncertainties of these
normalized values.

`D17Oscaling` is also capable of processing a list of “raw” individual (δ18O, Δ17O)
analyses, compute the average composition of each sample and assign uncertainties
to these average compositions based on the pooled variances of all δ18O and all Δ17O
measurements, respectively.
"""

import numpy as np
import pandas as pd
import uncertainties as uc

from pprint import pprint
from uncertainties import unumpy as unp

class DataSet:
	"""
	Store a table of individual analyses with starting ("raw") (δ18O, Δ17O) measurements
	and/or a table of analytes with starting (δ18O, Δ17O) compositions, with uncertainties.
	Provides methods to compile analytes from analyses, and to normalize the analyte
	compositions based on two or more reference materials.
	"""

	def __init__(self):
		pass

	def __repr__(self):
		s = vars(self)
		_ = []
		for a, v in s.items():
			if '\n' in str(v):
				_ += ['  '+f'[{a}] =\n{v}'.replace('\n', '\n    ')]
			else:
				_ += [f'  [{a}] = {v}']
		return 'DataSet(\n' + '\n'.join(_) + '\n)'

	def import_analyses(self, filepath, columns = {}):
		"""
		Read a CSV file with the following required fields:
			- `Analyte` (string)
			- `d18O_raw` (float-interpretable string, e.g., '-12.34')
			- `D17O_raw` (float-interpretable string, e.g., '0.1234')

		Arguments:
			* `filepath`: path of CSV file to import
			* `columns`: a dict of {new: old} column names
			(used to convert columns from the CSV to different names in the DataSet)

		Example:
			```
			data = DataSet()
			data.import_analyses(
				'data.csv',
				columns = {'d18O': 'd18O_raw', 'D17O': 'D17O_raw'},
			)
			```
		"""
		self.analyses = pd.read_csv(filepath).rename(columns = columns)
		self.compile()

	def import_analytes(self, filepath, columns = {}):
		"""
		Read a CSV file with the following required fields:
			- `Analyte` (string)
			- `d18O_raw` (UFloat-interpretable string, e.g., '-12.34±0.56')
			- `D17O_raw` (UFloat-interpretable string, e.g., '0.1234±0.0056')

		Arguments:
			* `filepath`: path of CSV file to import
			* `columns`: a dict of {new: old} column names
			(used to convert columns from the CSV to different names in the DataSet)

		Example:
			```
			data = DataSet()
			data.import_analytes(
				'data.csv',
				columns = {'d18O': 'd18O_raw', 'D17O': 'D17O_raw'},
			)
			```
		"""
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
		"Save `analyses` dataframe to CSV file"
		self.analyses.to_csv(filepath, index = False)

	def export_analytes(self, filepath):
		"Save `analytes` dataframe to CSV file"
		self.analytes.to_csv(filepath)

	def compile(self):
		"""
		Generate `analytes` dataframe based on the `analyses` dataframe.
		The δ18O and Δ17O values of each analyte are assigned based on the average over all analyses,
		with uncertainties equal to σ/sqrt(N), with σ being the pooled standard deviation of that field
		over all analytes and N being the number of analyses being averaged for this particular analyte.
		"""

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
		"""
		Nornalize the `analytes` dataframe based on two or more reference materials.

		Arguments:
			* refmaterials: a dict of nominal values of the form:
				{
					'VSMOW': {'d18O': 0., 'D17O': 0.},
					'SLAP': {'d18O': -55.5, 'D17O': 0.},
				}
			* λ (default = 0.528): the triple oxygen exponent used to compute Δ17O from (δ17O, δ18O)

		The normalized δ18O and Δ17O values are saved as UFloats in columns `d18O` and `D17O` of the
		`analytes` dataframe, and these UFloats are correlated as a result of the normalization corrections.
		Also saves `a18`, `b18`, `a17` and `b17` (all floats) corresponding to the affine transformations
		applied to `d18O_raw` and  `d17O_raw`, respectively.
		"""

		# populate self.refmaterials
		if hasattr(self, 'analytes'):
			loa = self.analytes.index.tolist()
		else:
			loa = self.analyses["Analyte"].unique().tolist()
		self.refmaterials = {k:v for k,v in refmaterials.items() if k in loa}

		# compute d17O of refmaterials
		for m in refmaterials.values():
			m['d17O'] = 1000 * float(
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
