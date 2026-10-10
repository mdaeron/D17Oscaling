# D17Oscaling

Perform two-point normalization of triple oxygen isotope measurements based on two or more reference materials with sufficiently different δ18O values.

The normalization applies two affine transformations:

- δ<sup>18</sup>O<sub>new</sub> = _a_<sub>18</sub> ∙ δ<sup>18</sup>O<sub>old</sub> + _b_<sub>18</sub>
- δ<sup>17</sup>O<sub>new</sub> = _a_<sub>17</sub> ∙ δ<sup>17</sup>O<sub>old</sub> + _b_<sub>17</sub>

Provided with (1) a set of measurements of one or more unknown samples, with uncertainties; (2) a set of measurement of two or more known reference materials, with uncertaines; (3) nominal (δ<sup>18</sup>O, Δ<sup>17</sup>O) values for the reference materials, `D17Oscaling` computes the normalized (δ<sup>18</sup>O, Δ<sup>17</sup>O) values for all unknown samples and the variance-covariance matrix describing the (correlated) uncertainties of these normalized values.

`D17Oscaling` is also capable of processing a list of “raw” individual (δ<sup>18</sup>O, Δ<sup>17</sup>O) analyses, compute the average composition of each sample and assign uncertainties to these average compositions based on the pooled variances of all δ<sup>18</sup>O and all Δ<sup>17</sup>O measurements, respectively.

## Example

```py
import D17Oscaling

```
