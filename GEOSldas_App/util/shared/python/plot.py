import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import basemap

def plotMap(
    data,
    *,
    ax=None,
    lat=None,
    lon=None,
    title=None,
    cRange=None,
    figsize=(8, 4),
    clbar=True,
    cRangeint=False,
    cmap=plt.cm.jet,
    bounding=None,
    prj="cyl",
):

    if cRange is not None:
        vmin = cRange[0]
        vmax = cRange[1]
    else:
        temp = flatData(data)
        vmin = np.percentile(temp, 5)
        vmax = np.percentile(temp, 95)
        if cRangeint is True:
            vmin = int(round(vmin))
            vmax = int(round(vmax))
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.subplots()

    if bounding is None:
        bounding = [
            np.min(lat) - 0.5,
            np.max(lat) + 0.5,
            np.min(lon) - 0.5,
            np.max(lon) + 0.5,
        ]

    mm = basemap.Basemap(
        llcrnrlat=bounding[0],
        urcrnrlat=bounding[1],
        llcrnrlon=bounding[2],
        urcrnrlon=bounding[3],
        projection=prj,
        resolution="c",
        ax=ax,
    )
    mm.drawcoastlines()
    #mm.drawstates(linestyle="dashed")
    #mm.drawcountries(linewidth=1.0, linestyle="-.")
    cs = mm.pcolormesh(lon, lat, data, cmap=cmap, vmin=vmin, vmax=vmax)
    if clbar is True:
        cb = mm.colorbar(cs, pad="5%", location="bottom")
        if 'normalized' in title:
            cb.set_ticks(np.linspace(vmin,vmax,6))
            cb.set_ticklabels([f'{10**x:.2f}' for x in np.linspace(vmin, vmax, 6)])

    if title is not None:
        ax.set_title(title)
    if ax is None:
        return fig, ax, mm
    else:
        return mm, cs

