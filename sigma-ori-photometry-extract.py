intab=[["IRAC-3", 5.731, 1.0], ["IRAC-4", 8, 1.0], ["WISE-4", 22.194, 1.176], ["WISE-3", 12.082, 0.04123], ["WISE-2", 4.618, 0.06087], ["PACS-R", 150, 4154.8], ["PACS-B", 70, 16619.2], ["MIPS-1", 23.68, 1.0], ["MIPS-2", 71.42, 1.0], ["MIPS-3", 155.9, 1.0], ["SPIRE-500", 500, 23.58], ["SPIRE-350", 350, 51.18], ["SPIRE-250", 250, 90.647], ["AKARI-N160", 160, 1.0], ["AKARI-WideL", 140, 1.0], ["AKARI-WideS", 90, 1.0], ["AKARI-N60", 65, 1.0], ["MSX-A", 8.28, "7.133e6"], ["MSX-C", 12.13, "2.863e7"], ["MSX-D", 14.65, "3.216e7"], ["MSX-E", 21.34, "2.476e7"]]
import json
import numpy as np
from astropy.table import Table

plotfile = "sigma-ori-seds.pdf"
wav_band = {k: float(x) for k, x, y in intab}
conv_band = {k: float(y) for k, x, y in intab}

data = json.load(open("sigma-ori-photometry.json"))

tabnames = []
for aperture, apdata in data.items():
    outdata = {
        "wav": [],
        "band": [],
        "I_nu": [],
        "d I_nu": []
    }
    for band, bdata in apdata.items():
        outdata["wav"].append(wav_band[band])
        outdata["band"].append(band)
        sb = np.mean(bdata["Source"])
        bg = np.mean(bdata["BG"])
        sig = np.std(bdata["BG"])
        print(aperture, band, sb-bg, conv_band[band], type(sb - bg))
        outdata["I_nu"].append(float(sb - bg)*conv_band[band])
        outdata["d I_nu"].append(sig*conv_band[band])

    tabname = f"sigma-ori-aperture-{aperture}.tab"
    outtab = Table(outdata)
    outtab.sort("wav")
    outtab.write(tabname, format="ascii.tab")
    tabnames.append(f"[[file:{tabname}]]")
