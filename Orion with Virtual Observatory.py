# ---
# jupyter:
#   jupytext_format_version: '1.2'
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
#   language_info:
#     codemirror_mode:
#       name: ipython
#       version: 3
#     file_extension: .py
#     mimetype: text/x-python
#     name: python
#     nbconvert_exporter: python
#     pygments_lexer: ipython3
#     version: 3.6.3
# ---

# ## Test of using VO to download data
#

# Set up the imports.

# + {"tags": ["boilerplate"], "slideshow": {"slide_type": "skip"}}
import matplotlib
import matplotlib.pyplot as plt
# %matplotlib inline  
import requests, io, astropy
from IPython.display import Image, display, HTML
import html
import seaborn as sns

## For handling ordinary astropy Tables
from astropy.table import Table

## For reading FITS files
from astropy.io import fits

## There are a number of relatively unimportant warnings that 
## show up, so for now, suppress them:
import warnings
warnings.filterwarnings("ignore")
# -

# Faster way of downloading files (hopefully)

from urllib.request import urlretrieve

# This seems to be a simpler approach than using `requests`, where we would have to mess about with chunking and the like.  See [SO answer](https://stackoverflow.com/a/44378512/353062)
#
# Example usage:
# ``` python
# url = 'http://mirror.pnl.gov/releases/16.04.2/ubuntu-16.04.2-desktop-amd64.iso'
# dst = 'ubuntu-16.04.2-desktop-amd64.iso'
# urlretrieve(url, dst)
# ```

# ### Getting images from SkyView
#
# This based on an [example from the AAS Workshop 2018](https://github.com/NASA-NAVO/aas_workshop_2018/blob/master/heasarc/heasarc_Image_Access.ipynb)

# We will start with LP Ori and try and get images from SkyView

# #### First, use Table Access Protocol to find resources that support Simple Image Access

# Warning:  'query' is bizarrely case-sensitive:
tap_params = {
    "request":"doQuery",
    "lang":"ADQL",
    "query":"""
        select b.short_name,b.res_description, c.access_url 
        from rr.capability a natural join rr.resource b natural join rr.interface c
        where a.cap_type='simpleimageaccess' and a.ivoid like 'ivo://nasa.heasarc%' 
        order by short_name;
    """
    }
r = requests.post('http://vao.stsci.edu/RegTAP/TapService.aspx/sync', data=tap_params)
## The astropy.table  module will read this VO Table formatted result into an astropy Table:
cat_table=Table.read(io.BytesIO(r.content))
#print(r.content)
## Just look at some non-so-randomly selected rows out of the 108 returned:
cat_table[85:96]
## Note that some of the other rows have special characters that confuse the print functions and cause errors.

# Get a list of all of the resources.  Note that the values in the table are bytes, not strings, so they need to be decoded before we can do much with them.

' '.join(_.decode() for _ in cat_table['short_name'])

# Same thing, but using `map` instead of a comprehension:

' '.join(map(bytes.decode, cat_table['short_name']))

# #### Go for the 2MASS JHK images first 

row = cat_table[cat_table['short_name'] == '2MASS'.encode()]
description = row[0]['res_description'].decode()
url = row[0]['access_url'].decode()
url, description

display(row[0]['access_url'].decode())

# +
# This doesn't work for LP Ori - need to find out what Vizier calls it
#import astropy.coordinates as coord
#coord = coord.SkyCoord.from_name("tet01 Ori C")
#pos = '{},{}'.format(coord.ra.deg,coord.dec.deg)
#coord
# -

pos = '83.791114,-5.4649843'
pixel_size = 1.0/3600 # 1 arcsec
size = 0.5
N = int(size/pixel_size)
params = {'POS': pos, 'SIZE': f'{size}', "NAXIS": f"{N},{N}"}
r = requests.get(url, params=params)
twomass_table=Table.read(io.BytesIO(r.content))
twomass_table

rowf = 0
rowj = rowf+1
jpgfile = 'lp_ori_test.jpg'
fitsfile = jpgfile.replace('.jpg', '.fits')

r=requests.get(twomass_table['URL'][rowj].decode("utf-8"), stream=True)
with open(jpgfile,'wb') as f:
    f.write(r.content)

display(Image(jpgfile))

for rowf, band in [0, 'H'], [2, 'J'], [4, 'K']:
    fitsfile = f'lp_ori_2mass_{band}.fits'
    r=requests.get(twomass_table['URL'][rowf].decode("utf-8"), stream=True)
    with open(fitsfile,'wb') as f:
        f.write(r.content)
        fits.open(fitsfile).info()

# #### Now try WISE

# The description seems to have had HTML entities escaped multiple times, so we use triple-nested `html.unescape` to unwrap all that.

row = cat_table[cat_table['short_name'] == 'WISE'.encode()]
description = row[0]['res_description'].decode('utf-8')
url = row[0]['access_url'].decode()
#display(HTML(description))
description = html.unescape(html.unescape(html.unescape(description)))
text = f'<div style="background-color:#cde;color:#733;padding:20px;">{description}</div>'
display(HTML(text))

pos = '83.791114,-5.4649843'
pixel_size = 2.0/3600 # 2 arcsec
size = 1.0
N = int(size/pixel_size)
params = {'POS': pos, 'SIZE': f'{size}', "NAXIS": f"{N},{N}"}
r = requests.get(url, params=params)
wise_table=Table.read(io.BytesIO(r.content))
wise_table

r = requests.get(wise_table['URL'][1].decode("utf-8"), stream=True)
with open(jpgfile,'wb') as f:
    f.write(r.content)

display(Image(jpgfile))

for rowf, band in [0, 'w1'], [2, 'w2'], [4, 'w3'], [6, 'w4']:
    fitsfile = f'data/Orion-VO/lp_ori_wise_{band}.fits'
    r=requests.get(wise_table['URL'][rowf].decode("utf-8"), stream=True)
    with open(fitsfile, 'wb') as f:
        f.write(r.content)
    fits.open(fitsfile).info()

# Possibly, there is new bow shock over in the east. At `5:36:38.099 -5:28:36.42`

# #### And thirdly, try Herschel SPIRE
#
# This is the longer wavelength instrument.  I think that the Herschel PACS data, which includes 70 $\mu$m has to come from IRSA, not from SkyView. 

row = cat_table[cat_table['short_name'] == 'HERSCHEL-SPIRE'.encode()]
description = row[0]['res_description'].decode('utf-8')
url = row[0]['access_url'].decode()
description = html.unescape(html.unescape(html.unescape(description)))
text = f'<div style="background-color:#cde;color:#733;padding:20px;">{description}</div>'
display(HTML(text))

# That was short and sweet!

pos = '83.791114,-5.4649843'
pixel_size = 2.0/3600 # 2 arcsec
size = 1.0
N = int(size/pixel_size)
params = {'POS': pos, 'SIZE': f'{size}', "NAXIS": f"{N},{N}"}
r = requests.get(url, params=params)
spire_table=Table.read(io.BytesIO(r.content))
spire_table

# *Nothing!* Oh well.

row = cat_table[cat_table['short_name'] == 'HAlpha'.encode()]
description = row[0]['res_description'].decode('utf-8')
url = row[0]['access_url'].decode()
description = html.unescape(html.unescape(html.unescape(description)))
text = f'<div style="background-color:#cde;color:#733;padding:20px;">{description}</div>'
display(HTML(text))

pos = '83.791114,-5.4649843'
pixel_size = 180.0/3600 # low res to start with
size = 10.0 # but big field
N = int(size/pixel_size)
params = {'POS': pos, 'SIZE': f'{size}', "NAXIS": f"{N},{N}"}
r = requests.get(url, params=params)
ha_table=Table.read(io.BytesIO(r.content))
ha_table

urlretrieve(ha_table['URL'][1].decode("utf-8"), jpgfile)
display(Image(jpgfile))

# So that worked, but it is very low resolution.  It shows the sigma orionis region better than M42.



# ### Getting images from IRSA
#
# This is necessary for the Spitzer and Herschel PACS data

# First, find all the catalogs.  ALthough this is not so useful, since I want images not catalogs.

# +
#Start with a Base URL.
tap_base_url = 'https://irsa.ipac.caltech.edu/TAP/sync?'

#Attach search parameters to this.
tap_params = {
    "request": "doQuery",
    "lang": "ADQL",
    "QUERY":"SELECT * FROM TAP_SCHEMA.tables"
    }

#Send the request.
r = requests.post(tap_base_url, data = tap_params)

#Put the results in an astropy table.
table=Table.read(io.BytesIO(r.content))

#Show the table.
table.show_in_notebook()
# -

# Now, try images.  But how to get a list of possible image types?  List is here: https://irsa.ipac.caltech.edu/ibe/sia.html

# +
#IRSA Simple Image Access Base URL
sia_base_url = 'https://irsa.ipac.caltech.edu/SIA?'
#Parameters of our search
sia_params = {
    "COLLECTION":"wise_allwise",
#    "POS":"circle -164.7 -5.8 0.1",
#    "BAND":"3e-6 4e-6",
#    "COLLECTION":"spitzer_sha",
    'POS': f'circle {pos[0]} {pos[1]} 0.5', 
    }

#Send the Request
r = requests.post(sia_base_url, data = sia_params)

#Put the results in an astropy table.
table=Table.read(io.BytesIO(r.content))

#Show the table.
table.show_in_notebook()
# -

pos









# ### Looking for my B star in catalogs

# +
scs_base_url = 'https://irsa.ipac.caltech.edu/SCS?'

scs_params = {'table': 'fp_psc', 'RA': 83.881677, 'DEC': -5.0040292, 'SR':30.0/3600.0}

r = requests.get(scs_base_url, params=scs_params)

table=Table.read(io.BytesIO(r.content))

table.show_in_notebook()
# -


