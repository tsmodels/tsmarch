# tsmarch 1.0.1

* Small fixes to RADICAL algorithm.

* Added the FASTICA algorithm as an option to the GOGARCH model since
the RADICAL algorithm still needs work for high dimensional systems
and may be slow in those cases. This is a custom implementation which
more closely follows the original Matlab code 
(https://research.ics.aalto.fi/ica/fastica/) rather than the other 
alternatives in R.

* Added the constant correlation test of Engle and Sheppard and a flextable
method for pretty printing.
