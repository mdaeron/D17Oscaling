# D17Oscaling

Perform two-point normalization of triple oxygen isotope measurements based on two or more reference materials with sufficiently different δ18O values.

The normalization applies two affine transformations:

- δ<sup>18</sup>O<sub>new</sub> = *a*<sub>18</sub> ∙ δ<sup>18</sup>O<sub>old</sub> + *b*<sub>18</sub>
- δ<sup>17</sup>O<sub>new</sub> = *a*<sub>17</sub> ∙ δ<sup>17</sup>O<sub>old</sub> + *b*<sub>17</sub>

Provided with (1) a set of measurements of one or more unknown samples, with uncertainties; (2) a set of measurement of two or more known reference materials, with uncertaines; (3) nominal (δ<sup>18</sup>O, Δ<sup>17</sup>O) values for the reference materials, `D17Oscaling` computes the normalized (δ<sup>18</sup>O, Δ<sup>17</sup>O) values for all unknown samples and the variance-covariance matrix describing the (correlated) uncertainties of these normalized values.

`D17Oscaling` is also capable of processing a list of “raw” individual (δ<sup>18</sup>O, Δ<sup>17</sup>O) analyses, compute the average composition of each sample and assign uncertainties to these average compositions based on the pooled variances of all δ<sup>18</sup>O and all Δ<sup>17</sup>O measurements, respectively.

## Usage

First import `D17Oscaling` and define some reference materials:

```py
from D17Oscaling import DataSet

RM = {
	'VSMOW': {'d18O': 0., 'D17O': 0.},
	'SLAP': {'d18O': -55.5, 'D17O': 0.},
}
```

Then create a `DataSet()` object, import some analytical data and print out what the DataSet looks like now:

```py
# create DataSet() instance
data = DataSet()

# import some analytical data;
data.import_analyses(
	'examples/data.csv',
	columns = {
		'd18O': 'd18O_raw',
		'D17O': 'D17O_raw'
	},
)
# `columns` renames columns to 'd18O_raw' and 'D17O_raw',
# meaning that these values are not yet normalized.

# take a look
print(data)
```

<p align="center">
  <img src="https://raw.githubusercontent.com/mdaeron/D17Oscaling/refs/heads/master/warm-arrow-down.png" width="80" alt="down arrow">
</p>

This should look like this.
- The raw data we just imported is stored as a `pandas.DataFrame` in `data.analyses`.
- The pooled standard deviantion of `d18O_raw` and `D17O_raw` for this dataset is stored as `data.SD_d18O_raw` and `data.SD_D17O_raw`.
- The mean raw values for each analyte were automatically computed, with uncertainties derived from `SD_d18O_raw`, `SD_D17O_raw`, and the number of replicate analyses for each analyte. These mean values are saved as `uncertainties.UFloat` values in a `pandas.DataFrame` stored as `data.analytes`.

```txt
DataSet(
  [analyses] =
        UID Analyte   d18O_raw  D17O_raw
    0   X01     FOO -12.532655 -0.098675
    1   X02     BAR  -7.879440  0.097711
    2   X03     BAZ -50.274837 -0.126091
    3   X04   VSMOW  -0.085304  0.004308
    4   X05    SLAP -55.412060  0.021416
    5   X06     FOO -12.467660 -0.108063
    6   X07     BAR  -7.841212  0.088934
    7   X08     BAZ -50.160433 -0.138356
    8   X09   VSMOW   0.046751  0.006160
    9   X10    SLAP -55.585929  0.011290
    10  X11     FOO -12.508564 -0.105138
    11  X12     BAR  -7.943704  0.085655
    12  X13     BAZ -50.185312 -0.138463
    13  X14   VSMOW  -0.004993  0.006506
    14  X15    SLAP -55.518486  0.007433
    15  X16     FOO -12.613532 -0.098567
    16  X17     BAR  -7.725561  0.087401
    17  X18     BAZ -50.288610 -0.127896
    18  X19   VSMOW  -0.042833  0.001167
    19  X20    SLAP -55.535213  0.002187
  [SD_d18O_raw] = 0.07035790333040727
  [SD_D17O_raw] = 0.005786843071251322
  [analytes] =
             N         d18O_raw          D17O_raw
    Analyte
    BAR      4   -7.847+/-0.035   0.0899+/-0.0029
    BAZ      4  -50.227+/-0.035  -0.1327+/-0.0029
    FOO      4  -12.531+/-0.035  -0.1026+/-0.0029
    SLAP     4  -55.513+/-0.035   0.0106+/-0.0029
    VSMOW    4   -0.022+/-0.035   0.0045+/-0.0029
)
```

We can then normalize `data.analytes` based on the nominal compositions of the reference materials.

```py
data.normalize(RM)
print(data)
```

<p align="center">
  <img src="https://raw.githubusercontent.com/mdaeron/D17Oscaling/refs/heads/master/warm-arrow-down.png" width="80" alt="down arrow">
</p>

As seen below, this adds two columns named `d18O` and `D17O` to `data.analytes`, corresponding to the normalized compositions of each analyte. These values are all `uncertainties.UFloat` values, so that their correlated uncertainties are all explictly computed and accessible using the methods described in the [uncertainties documentation](https://pythonhosted.org/uncertainties). The parameters *a*<sub>18</sub>, *b*<sub>18</sub>, *a*<sub>17</sub>, and *b*<sub>17</sub> (see equations above) are also stored as `uncertainties.UFloat` values in `data.a18`, `data.b18`, `data.a17`, and `data.b17`.

```txt
DataSet(
  [analyses] =
        UID Analyte   d18O_raw  D17O_raw
    0   X01     FOO -12.532655 -0.098675
    1   X02     BAR  -7.879440  0.097711
    2   X03     BAZ -50.274837 -0.126091
    3   X04   VSMOW  -0.085304  0.004308
    4   X05    SLAP -55.412060  0.021416
    5   X06     FOO -12.467660 -0.108063
    6   X07     BAR  -7.841212  0.088934
    7   X08     BAZ -50.160433 -0.138356
    8   X09   VSMOW   0.046751  0.006160
    9   X10    SLAP -55.585929  0.011290
    10  X11     FOO -12.508564 -0.105138
    11  X12     BAR  -7.943704  0.085655
    12  X13     BAZ -50.185312 -0.138463
    13  X14   VSMOW  -0.004993  0.006506
    14  X15    SLAP -55.518486  0.007433
    15  X16     FOO -12.613532 -0.098567
    16  X17     BAR  -7.725561  0.087401
    17  X18     BAZ -50.288610 -0.127896
    18  X19   VSMOW  -0.042833  0.001167
    19  X20    SLAP -55.535213  0.002187
  [SD_d18O_raw] = 0.07035790333040727
  [SD_D17O_raw] = 0.005786843071251322
  [analytes] =
             N         d18O_raw          D17O_raw         d17O_raw           d18O             d17O            D17O
    Analyte
    BAR      4   -7.847+/-0.035   0.0899+/-0.0029   -4.062+/-0.019   -7.83+/-0.05   -4.056+/-0.025   0.085+/-0.023
    BAZ      4  -50.227+/-0.035  -0.1327+/-0.0029  -26.971+/-0.019  -50.21+/-0.05  -26.974+/-0.026  -0.143+/-0.026
    FOO      4  -12.531+/-0.035  -0.1026+/-0.0029   -6.738+/-0.019  -12.51+/-0.05   -6.733+/-0.024  -0.109+/-0.022
    SLAP     4  -55.513+/-0.035   0.0106+/-0.0029  -29.695+/-0.019  -55.50+/-0.05  -29.699+/-0.027   0.000+/-0.028
    VSMOW    4   -0.022+/-0.035   0.0045+/-0.0029   -0.007+/-0.019    0.00+/-0.05    0.000+/-0.027   0.000+/-0.027
  [refmaterials] = {'VSMOW': {'d18O': 0.0, 'D17O': 0.0, 'd17O': 0.0}, 'SLAP': {'d18O': -55.5, 'D17O': 0.0, 'd17O': -29.69864899849639}}
  [a18] = 1.0002+/-0.0009
  [b18] = 0.022+/-0.035
  [a17] = 1.0003+/-0.0009
  [b17] = 0.007+/-0.019
)
```

Alternatively, one may perform the same normalization steps based directly on the analytes table, skipping the computation of `d18O_raw` and `D17O_raw` uncertainties (assigning these uncertainties is then solely the user's responsibility).

```py
data = DataSet()
data.import_analytes('examples/analytes.csv')

print(data)
```

<p align="center">
  <img src="https://raw.githubusercontent.com/mdaeron/D17Oscaling/refs/heads/master/warm-arrow-down.png" width="80" alt="down arrow">
</p>

```txt
DataSet(
  [analytes] =
             N         d18O_raw          D17O_raw
    Analyte
    BAR      4   -7.847+/-0.035   0.0899+/-0.0029
    BAZ      4  -50.227+/-0.035  -0.1327+/-0.0029
    FOO      4  -12.531+/-0.035  -0.1026+/-0.0029
    SLAP     4  -55.513+/-0.035   0.0106+/-0.0029
    VSMOW    4   -0.022+/-0.035   0.0045+/-0.0029
)
```

And the normalized results are identical to those in the first example:

```py
data.normalize(RM)
print(data)
```

<p align="center">
  <img src="https://raw.githubusercontent.com/mdaeron/D17Oscaling/refs/heads/master/warm-arrow-down.png" width="80" alt="down arrow">
</p>

```txt
DataSet(
  [analytes] =
             N         d18O_raw          D17O_raw         d17O_raw           d18O             d17O            D17O
    Analyte
    BAR      4   -7.847+/-0.035   0.0899+/-0.0029   -4.061+/-0.019   -7.83+/-0.05   -4.056+/-0.025   0.085+/-0.023
    BAZ      4  -50.227+/-0.035  -0.1327+/-0.0029  -26.971+/-0.019  -50.21+/-0.05  -26.974+/-0.026  -0.143+/-0.025
    FOO      4  -12.531+/-0.035  -0.1026+/-0.0029   -6.738+/-0.019  -12.51+/-0.04   -6.733+/-0.024  -0.108+/-0.022
    SLAP     4  -55.513+/-0.035   0.0106+/-0.0029  -29.695+/-0.019  -55.50+/-0.05  -29.699+/-0.027   0.000+/-0.028
    VSMOW    4   -0.022+/-0.035   0.0045+/-0.0029   -0.007+/-0.019    0.00+/-0.05    0.000+/-0.026   0.000+/-0.026
  [refmaterials] = {'VSMOW': {'d18O': 0.0, 'D17O': 0.0, 'd17O': 0.0}, 'SLAP': {'d18O': -55.5, 'D17O': 0.0, 'd17O': -29.69864899849639}}
  [a18] = 1.0002+/-0.0009
  [b18] = 0.022+/-0.035
  [a17] = 1.0003+/-0.0009
  [b17] = 0.007+/-0.019
)
```
