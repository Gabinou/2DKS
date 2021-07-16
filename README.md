# KS2D
## 2 Dimensional Kolmogorov-Smirnov test for goodness-of-fit.

KS2D is a two-dimensional extension to the Kolmogorov-Smyrnov test for goodness-of-fit.
It is used to compare datasets of points to a distributions, or two datasets of points.

In this case, we check if two-dimensional data fits a particular distribution. The extension to higher dimensions is non-trivial and requires on the order of *O*(*n*<sup>2</sup>) operations i.e. slow for large datasets.

Mainly intended to be interacted with using functions ``ks2d1s`` and ``ks2d2s``, which take as inputs one 2-column matrix and one 2D function, or two 2-column matrices respectively. 
These algorithms compute the relative probabilities of finding data in orthonormal quadrants that surround each point in the data set, then uses those to compute the K-S statistic with its distribution function (`Qks`). Look around 14.3.7 and 14.7.1 in [3] for detailed arcane mathemagic explanations. 

``ks2d1s`` and ``ks2d2s`` outputs the KS statistic *D*, and the significance level *prob* of that observed value of *D*, respectively.

## Issues

Float number representation and rounding: probabilities expected to sum to 1.0 return 0.99999999 instead, etc... No plans to implement any kind of solution to this: it sounds much more trouble than it is worth. This test is theoretically just an approximation, rounding to a couple digits seem reasonable.

``ks2d1s``: Input: no 2D density matrices. Integration method: only numerical.

## MISC

Prerequisites: *scipy*, *numpy*.

Keywords: Kolmogorov-smirnov test, multi-dimensional, two-dimensional, KS test, hypothesis testing, statistics, probability theory.

Based on the idea by Peacock (1983), an upgrade by Fasano and Franceschin (1987) with
much guidance from Numerical recipes in C by Press and Teukolsky (1996).

# References
[1] Peacock, J. A. (1983). Two-dimensional goodness-of-fit testing in astronomy. *Monthly Notices of the Royal Astronomical Society*, 202(3), 615-627.

[2] Fasano, G., & Franceschini, A. (1987). A multidimensional version of the Kolmogorov–Smirnov test. *Monthly Notices of the Royal Astronomical Society*, 225(1), 155-170.

[3] Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (1996). *Numerical recipes in C (Vol. 2).* Cambridge: Cambridge university press.
