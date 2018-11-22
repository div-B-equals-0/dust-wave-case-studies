# ---
# jupyter:
#   jupytext_format_version: '1.2'
#   kernelspec:
#     display_name: Python [conda env:astroconda]
#     language: python
#     name: conda-env-astroconda-py
#   language_info:
#     codemirror_mode:
#       name: ipython
#       version: 3
#     file_extension: .py
#     mimetype: text/x-python
#     name: python
#     nbconvert_exporter: python
#     pygments_lexer: ipython3
#     version: 3.5.1
# ---

import numpy as np
import pysynphot as S

bp1 = S.ObsBandpass('acs,wfc1,f555w')

bp1.unit_response()

bp1.obsmode.modes

S.ObsBandpass?
