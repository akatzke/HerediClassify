.. HerediClassify documentation master file, created by
   sphinx-quickstart on Fri Nov  8 16:25:56 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

HerediClassify documentation
============================


Welcome to the HerediClassify documentation!

HerediClassify is a tool for automated variant classification in hereditary breast and ovarian cancer.

HerediClassify follows the ACMG/AMP guidelines and implements 19/28 criteria.
Additional gene specific recommendations for *ATM* (v.1.3.0), *BRCA1* (v.1.1.0), *BRCA2* (v.1.1.0), *CDH1* (v.3.1.0), "PALB2" (v.1.1.0), "PTEN" (v.3.1.0), and "TP53" (v.1.4.0).
Uniquely, HerediClassify introduced a clear separation between evidence on protein and splicing level.
For this purpose, many rules have more than one implementation to separate between the evidence types.

..
    Insert image of implemented

For some example data of variants classified by HerediClassify please see the download section of the HerediVar website.


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   execution
   moduls
