import numpy as np
import os
from equalArea import mcbryde
import matplotlib.pyplot as plt

def getDesFootprint () :
    import insideDesFootprint
    ra,dec = insideDesFootprint.getFootprintRaDec() 
    return ra,dec
def plotDesFootprint(alpha, beta, xmin, xmax, ymin, ymax, ax) :
    import matplotlib.pyplot as plt
    import matplotlib.path
    import matplotlib.patches
    from equalArea import mcbryde
    desRa, desDec = getDesFootprint()
    x,y = mcbryde.mcbryde(desRa, desDec, alpha=alpha, beta=beta)
    ix = (x > xmin) & (x < xmax) & (y > ymin) & (y < ymax)
    #plt.plot(x[ix],y[ix],c="k",alpha=0.5)
    footprint = matplotlib.path.Path(zip(x,y))
    patch = matplotlib.patches.PathPatch(footprint, facecolor='gold', lw=1, alpha=0.066, fill=True)
    ax.add_patch(patch)
    patch = matplotlib.patches.PathPatch(footprint, edgecolor='gold', lw=1, alpha=1, fill=False)
    #patch = matplotlib.patches.PathPatch(footprint, edgecolor='gold', lw=1, alpha=0.33, fill=False)
    ax.add_patch(patch)
    
#  raHex, decHex = np.genfromtxt("512/lmc-ra-dec-prob-mjd-slot.txt", unpack=True, usecols=(0,1))
#  raHex, decHex = plotMapAndHex.getG184098_iband_hexes() 
#  alpha,beta= plotMapAndHex.mapAndHex(figure,"lmc",10,"512/",10,raHex, decHex); 
#  plotMapAndHex.lmcLabels()
#
# exp,raHex,decHex=np.genfromtxt("ra-dec.txt",unpack=True,comments="#")
# reload(plotMapAndHex);plotMapAndHex.mapAndHex(figure, "G211117",0,"/data/des30.a/data/annis/des-gw/Christmas16-event/maps/", 8, raHex, decHex, "" )
#
def mapAndHex(figure, simNumber, slot, data_dir, nslots, hexRa, hexDec, camera,
        title="", colorbar=True, slots=np.zeros(0), doHexes=True, allSky=False) :
    import healpy as hp
    import hp2np

    # get the planned observations
    ra,dec = hexRa, hexDec
    
    # get the maps for a reasonable slot
    raMap, decMap, ligoMap, maglimMap, probMap, \
        haMap, xMap, yMap, hxMap, hyMap = readMaps(data_dir, simNumber, slot)

    resolution=256
    resolution=512
    doStars = False
    image = False
    image = True
    scale = 3.
    scale = 1.
    redRa = 90.

    raMid = -1000
    #raBoxSize = 0.
    #decBoxSize = 16.
    #mod_ra = -12
    #mod_dec = 5
    mod_ra = 0; mod_dec=3; raBoxSize=15; decBoxSize=15
    raBoxSize=60

    low_limit = 21.; high_limit=23.8
    if doStars :
        low_limit = 0.0; high_limit=3.2
        low_limit = -1000.0; high_limit=750.
        low_limit = -100.0; high_limit=150.
        low_limit = -1.3; high_limit=1.8

        stars = hp.read_map("/home/s1/annis/daedalean/desgw-map/data/2MASS-J-lt-16-stardensity-hp.fits")
        #stars = hp.read_map("/home/s1/annis/daedalean/desgw-map/data/2MASS-J-lt-16-starcounts-hp.fits")
        if resolution != 512 :
            junk, junk, stars = hp2np.map2np(stars, resolution)
        ix = stars <= 0
        stars[ix]= 1.
        stars = np.log10(stars)

    if doStars:
        map = stars
        badData = False
    else :
        map = maglimMap
        badData = True
        badDataVal = -11.0

    doOrigLigoMap = False
    if doOrigLigoMap :
        ligoResolution = hp.get_nside(ligoMap)
        mapName = data_dir+"skyprobcc_cWB_complete.fits"
        origLigoRa, origLigodDec, origLigoMap = getOriginalLigoMap (mapName, ligoResolution) 
        
    else :
        #raise Exception("why are you doing this?")
        origLigoMap = ligoMap
    alpha, beta = coreMapAndHex(figure, ra, dec, raMap, decMap, camera, map, \
        low_limit, high_limit, ligoMap, origLigoMap, doLigoMap=True, doOrigLigoMap=doOrigLigoMap, \
        resolution=resolution, image=image, scale=scale, badData=badData, badDataVal=badDataVal, \
        redRa = redRa, title=title, raMid=raMid, raBoxSize=raBoxSize, decBoxSize = decBoxSize, \
        mod_ra=mod_ra, mod_dec= mod_dec , colorbar=colorbar, slots=slots, thisSlot=slot, 
        doHexes=doHexes, allSky = allSky)

    return alpha,beta

def coreMapAndHex(figure, hexRa, hexDec, raMap, decMap, camera, map, 
        low_limit, high_limit, ligoMap, origLigoMap, doLigoMap=True, doOrigLigoMap=False, 
        resolution=512, image=False, scale=1., badData=False, badDataVal=-11.0,
        redRa = 90., title="", raMid=-1000, raBoxSize=5., decBoxSize=5., mod_ra = 0, mod_dec=0.,
        doHexes = True, gradRedHiDec = -80, raGratDelRa=30., decGratDelDec=10. , colorbar=True,
        contourLabels=True , slots=np.zeros(0), thisSlot=0, allSky = False) :
    from equalArea import mcbryde
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import insideDesFootprint
    import matplotlib
    from scipy.ndimage.filters import gaussian_filter

    cmap = "cubehelix_r"
    cmap = "YlGnBu"
    cmap = matplotlib.cm.get_cmap("YlGnBu")
    cmap.set_under("w")
    cmap.set_bad("0.30")  ;# a dark gray

    # compute the image limits and midpoints (alpha, beta)
    raMin, raMax, decMin, decMax, xmin, xmax, ymin, ymax, alpha, beta = \
        computeLimits (hexRa, hexDec, raMid=raMid, raBoxSize=raBoxSize, 
        decBoxSize=decBoxSize, mod_ra=mod_ra, mod_dec = mod_dec, allSky = allSky) 
    #j,j,s_hexRa, s_hexDec = lmcHexes()
    #raMin, raMax, decMin, decMax, xmin, xmax, ymin, ymax, alpha, beta = \
    #    computeLimits (s_hexRa,s_hexDec, raMid=raMid, raBoxSize=raBoxSize, 
    #    decBoxSize=decBoxSize, mod_ra=mod_ra, mod_dec = mod_dec) 

    # project into mcbryde 
    xMap,yMap = mcbryde.mcbryde(raMap, decMap, alpha=alpha, beta=beta)
    x,y = mcbryde.mcbryde(hexRa, hexDec, alpha=alpha, beta=beta)

    if xmin==xmax and ymin==ymax:
        ix = np.ones(xMap.size).astype(bool)
        xmin=xMap.min(); xmax=xMap.max()
        ymin=yMap.min(); ymax=yMap.max()
        doHexes=False
    else :
        ix=np.nonzero((xMap > xmin) & (xMap  <= xmax) & (yMap > ymin) & (yMap <= ymax) )

    # plot the image, either as an image or as a hexbin
    plt.clf()
    if image :
        data = makeImage (xMap[ix], yMap[ix], map[ix], xmin, xmax, ymin, ymax, scale, 
            badData=badData, badDataVal=badDataVal)
        # hack to make this 5 sigma, not 10 sigma limitin mag
        print "\t\t 10sigma -> 5 sigma hack",
        data = data +0.75257 
        low_limit = low_limit+ 0.75257
        high_limit = high_limit+ 0.75257
        plt.imshow(data, cmap=cmap, vmin=low_limit, vmax=high_limit, 
            extent=[xmin, xmax, ymin, ymax])
    else :
        gridsize = 150
        gridsize = 66
        plt.hexbin(xMap[ix],yMap[ix],map[ix],vmin=low_limit, vmax=high_limit, gridsize=gridsize,
            cmap=cmap)

    # put on a graticule
    graticule (alpha, beta, xmin, xmax, ymin, ymax,  redRa = redRa, redRaDec2 = gradRedHiDec,
        raGratDelRa= raGratDelRa, decGratDelDec= decGratDelDec)

    # fig 1
    ax = figure.add_subplot(1,1,1)
    #Plot Des Footprint in the observing plots
    plotDesFootprint(alpha, beta, xmin, xmax, ymin, ymax, ax)

    if colorbar and not allSky:
        cb = plt.colorbar(shrink=0.5,pad=0.03); 
        cb.set_label("5$\sigma$ point source limiting magnitude")

    # put on the ligo contours
    gaussSigma = 7.
    gaussSigma = 3.
    if doLigoMap :
        cl_ligoMap = confidenceLevels(ligoMap)
        if gaussSigma > 0 :
            cl_ligoMap = gaussian_filter(cl_ligoMap, gaussSigma)
        plotLigoContours(xMap[ix],yMap[ix], cl_ligoMap[ix], color="w", lw=1.2, labels=contourLabels) 
    if doOrigLigoMap :
        cl_origLigoMap = confidenceLevels(origLigoMap)
        if gaussSigma > 0 :
            cl_origLigoMap = gaussian_filter(cl_origLigoMap, gaussSigma)
        plotLigoContours(xMap[ix],yMap[ix], cl_origLigoMap[ix], color="cyan", alpha=1.0, lw=1.3,
            ls="dashed",labels=contourLabels) 
        #plotLigoContours(xMap[ix],yMap[ix], cl_origLigoMap[ix], color="cyan", alpha=1.0, lw=1.3,
            #labels=contourLabels) 

    # put on the hexes
    if doHexes :
        linewidth=0.5
        linewidth=1.0
        #ax = figure.add_subplot(1,1,1)
        # this is needed for search fig 1
        ax=plotDecamHexen(ax, hexRa, hexDec, alpha, camera, beta, color="r", lw=linewidth, allSky=allSky) 
        #ix =np.invert( insideDesFootprint.insideFootprint(hexRa, hexDec))
        #ax=plotDecamHexen(ax, hexRa[ix],hexDec[ix],alpha, beta, color="orange", lw=linewidth, allSky=allSky) 
        if slots.size > 0 :
            # plot the already observed hexes as maroon
            ix = slots<thisSlot
            ax=plotDecamHexen(ax, hexRa[ix],hexDec[ix],alpha, camera, beta, color="maroon", lw=linewidth, allSky=allSky) 
            # plot the current slots hexes as yellow
            ix = slots==thisSlot
            ax=plotDecamHexen(ax, hexRa[ix],hexDec[ix],alpha, camera, beta, color="yellow", lw=linewidth, allSky=allSky) 
        

        # fig1 and fig2, lmc paper, 
        #plotLmcHexes(alpha,beta,ax)
        # and not the fig1 marcelle paper, which is change lmcHexes2 to lmcHexes


    # deal with titles, axes, etc
    plt.title(title)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.axes().set_aspect('equal'); 
    #plt.axes().set_frame_on(False); 
    plt.axes().set_xticks([]); 
    plt.axes().set_yticks([])

    plt.show()

    return alpha, beta

# compute image limits and midpoints (alpha, beta)
def computeLimits (raHex, decHex, raMid = -1000, raBoxSize=5., decBoxSize=5, mod_ra=0, mod_dec=0, allSky = False) :
    from equalArea import mcbryde
    verbose = False
    mod_alpha = 0
    mod_beta = 0
    if allSky :
        decMin = -89.9999; decMax = 89.9999
        raMin = -179.9999; raMax = 179.9999
        raBoxSize = 0.; decBoxSize = 0.
        mod_dec = 0; mod_ra = 0
    else :
        # autoscaling
        decMin = decHex.min(); decMax = decHex.max()
        raMin = raHex.min(); raMax = raHex.max()
    decMin = decMin-decBoxSize+mod_dec
    decMax = decMax+decBoxSize+mod_dec
    decMid = decMin+(decMax-decMin)/2.
    beta=-1*(decMid+mod_beta)
    boxRa = raBoxSize/np.cos(decMid*2*np.pi/360.)
    raMin = raMin-boxRa+mod_ra
    raMax = raMax+boxRa+mod_ra
    if raMid == -1000 :
        raMid = raMin+(raMax-raMin)/2.
    alpha= -1*(raMid+mod_alpha)

    print "\t raMin, raMax, decMin, decMax:", raMin, raMax, decMin, decMax
    #raise Exception
    if allSky :
        #x,y = mcbryde.mcbryde(np.array([-179.999,179.999]), np.array([-60,60]), alpha=alpha, beta=beta)
        #decMin = -60; decMax = 60.0
        x,y = mcbryde.mcbryde(np.array([raMin,raMax]), np.array([-1,1]), alpha=alpha, beta=beta)
        xmin = x.min(); xmax = x.max()
        x,y = mcbryde.mcbryde(np.array([-1,1]), np.array([decMin,decMax]), alpha=alpha, beta=beta)
        ymin = y.min(); ymax = y.max()
        xbox = 0.01*(xmax-xmin)
        ybox = 0.01*(ymax-ymin)
    else :
        # autoscaling
        #x,y = mcbryde.mcbryde(raHex, decHex, alpha=alpha, beta=beta)
        x,y = mcbryde.mcbryde(np.array([raMin,raMax]), np.array([1,1]), alpha=alpha, beta=beta)
        xmin = x.min(); xmax = x.max()
        x,y = mcbryde.mcbryde(np.array([-1,1]), np.array([decMin,decMax]), alpha=alpha, beta=beta)
        ymin = y.min(); ymax = y.max()
        xbox = 0.1*(xmax-xmin)
        ybox = 0.1*(ymax-ymin)
    xmin = xmin-xbox; xmax = xmax+xbox
    ymin = ymin-ybox; ymax= ymax+ybox

    if verbose:
        print "ra box, dec box",raMin, raMax, decMin, decMax
        print "x box, y box", xmin, xmax, ymin, ymax 
        print "alpha, beta",alpha, beta
    #raw_input("wait for human to press")
    return raMin, raMax, decMin, decMax, xmin, xmax, ymin, ymax, alpha, beta

# @profile
def makeImage (xMap, yMap, vals, xmin, xmax, ymin, ymax, scale, 
        badData=False, badDataVal=-11.0, verbose=False, too_far_away_scale=1.5) :
    import scipy.spatial
    tree = scipy.spatial.KDTree(zip(xMap, yMap))

    xsize = int(xmax)+1 - int(xmin)-1 
    ysize = int(ymax)+1 - int(ymin)-1 
    nsteps_x = int(xsize*scale)
    nsteps_y = int(ysize*scale)
    step_size = 1.0/scale

    data = np.zeros( (nsteps_y, nsteps_x) )
    
    if verbose: print xmin, xmax, ymin, ymax, "    ", xsize, ysize
    
    # this is "supposed" to be slow
    for i in range(0, nsteps_x) :
        for j in range (0, nsteps_y) :
            xpix= i*step_size  + xmin
            ypix= j*step_size  + ymin 
            
            d= tree.query([xpix,ypix])
            dist = d[0]; idx=d[1]
            val = vals[idx]
            if dist > too_far_away_scale: val=badDataVal
            data[nsteps_y-j-1,i] = val
    
    if badData:
        data = np.ma.masked_equal(data, badDataVal)
    if verbose: print data.mean()
    return data

def getOriginalLigoMap (mapName, resolution) :
    import healpy as hp
    import hp2np
    secondLigo = hp.read_map(mapName)
    secondRa,secondDec,secondLigo = hp2np.map2np(secondLigo, resolution)
    return secondRa, secondDec, secondLigo

# one could put raGratDec2 = 50 and decGratDec2 = 50 for a smaller plot
def graticule (alpha, beta, xmin, xmax, ymin, ymax, 
         raGratRa1=-179.99, raGratRa2=180+30, raGratDelRa=30., 
        raGratDec1=-89.99, raGratDec2=91, raGratDelDec= 0.1 ,
        decGratRa1=-179.99, decGratRa2=180+30, decGratDelRa=0.1, 
        decGratDec1=-89.99, decGratDec2=91, decGratDelDec=10,
        redRa = 90., redRaDec1=-90, redRaDec2=-80, redRaDelDec = 0.1) :
    import matplotlib.pyplot as plt
    from equalArea import mcbryde
    for i in np.arange( raGratRa1, raGratRa2, raGratDelRa ) :
        raLine, decLine = np.array([]), np.array([])
        for j in np.arange( raGratDec1, raGratDec2, raGratDelDec ) :
            raLine = np.append(raLine, i)
            decLine = np.append(decLine, j)
        xLine,yLine = mcbryde.mcbryde(raLine, decLine, alpha=alpha, beta=beta)
        ixg = (xLine > xmin) & (xLine < xmax) & (yLine > ymin) & (yLine < ymax)
        plt.plot(xLine[ixg],yLine[ixg],c="k",alpha=0.5, linewidth=0.5)
    doRedLine = True
    doRedLine = False
    if doRedLine :
        i = redRa
        raLine, decLine = np.array([]), np.array([])
        for j in np.arange( redRaDec1, redRaDec2, redRaDelDec ) :
            raLine = np.append(raLine, i)
            decLine = np.append(decLine, j)
        xLine,yLine = mcbryde.mcbryde(raLine, decLine, alpha=alpha, beta=beta)
        ixg, = np.where((xLine > xmin) & (xLine < xmax) & (yLine > ymin) & (yLine < ymax))
        if ixg.size > 0:
            plt.plot(xLine[ixg],yLine[ixg],c="r",alpha=1.0, linewidth=0.5)

    for i in np.arange( decGratDec1, decGratDec2, decGratDelDec ) :
        raLine, decLine = np.array([]), np.array([])
        for j in np.arange( decGratRa1, decGratRa2, decGratDelRa ) :
            raLine = np.append(raLine, j)
            decLine = np.append(decLine, i)
        xLine,yLine = mcbryde.mcbryde(raLine, decLine, alpha=alpha, beta=beta, isLine=True)
        ixg, = np.where((xLine > xmin) & (xLine < xmax) & (yLine > ymin) & (yLine < ymax))
        if ixg.size > 0:
            plt.plot(xLine[ixg],yLine[ixg],c="k",alpha=0.5, linewidth=0.5)

def plotDecamHexen(ax, ra,dec,alpha, camera, beta=0, color="r", lw=1, plateCaree=False, allSky=False) :
    if camera == 'decam':
        import decam2hp
        import matplotlib.patches 
        import matplotlib.path 
        from equalArea import mcbryde
        import matplotlib.pyplot as plt
        nHex = ra.size
        for i in range(0,nHex) :
            hexRa,hexDec = decam2hp.cameraOutline(ra[i], dec[i], camera)
            hexX,hexY = mcbryde.mcbryde(hexRa,hexDec, alpha=alpha, beta=beta, )
            if plateCaree:
                hexX,hexY = hexRa, hexDec
            if not allSky :
                hex_path = matplotlib.path.Path(zip(hexX,hexY))
                hex_patch = matplotlib.patches.PathPatch(hex_path, edgecolor=color, lw=lw, fill=False)
                ax.add_patch(hex_patch)
            else :
                # long section dealing with hexes that get split across the singularity
                # and cause lines from one side of map to another
                # split them into separate entities
                # ugly, but it seems to work.
                ix_pos = np.nonzero(hexRa>180)[0]
                ix_neg = np.nonzero(hexRa<-180)[0]
                ix, = np.where(np.nonzero((hexRa>=-180)&(hexRa<=180))[0])
                if ix.size > 0 :
                    hex_path = matplotlib.path.Path(zip(hexX[ix],hexY[ix]))
                    hex_patch = matplotlib.patches.PathPatch(hex_path, edgecolor=color, lw=lw, fill=False)
                    ax.add_patch(hex_patch)
                if ix_pos.size > 0:
                    hex_path = matplotlib.path.Path(zip(hexX[ix_pos],hexY[ix_pos]))
                    hex_patch = matplotlib.patches.PathPatch(hex_path, edgecolor=color, lw=lw, fill=False)
                    ax.add_patch(hex_patch)
                if ix_neg.size > 0:
                    hex_path = matplotlib.path.Path(zip(hexX[ix_neg],hexY[ix_neg]))
                    hex_patch = matplotlib.patches.PathPatch(hex_path, edgecolor=color, lw=lw, fill=False)
                    ax.add_patch(hex_patch)
                
            #x,y=mcbryde.mcbryde(tra[i],tdec[i], alpha=alpha, beta=beta)
            #plt.text(x,y,"{}".format(i), ha="center", va="center", color="w")
    if camera == 'hsc':
        import decam2hp
        import matplotlib.patches
        import matplotlib.path
        from equalArea import mcbryde
        import matplotlib.pyplot as plt
        nHex = ra.size
        radius = 1.5 / 2
        for i in range(0,nHex) :
            hexRa,hexDec = decam2hp.cameraOutline(ra[i], dec[i], camera)
            hexX,hexY = mcbryde.mcbryde(hexRa,hexDec, alpha=alpha, beta=beta, )
            ix_pos = np.nonzero(hexRa>180)[0]
            ix_neg = np.nonzero(hexRa<-180)[0]
            ix, = np.where(np.nonzero((hexRa>=-180)&(hexRa<=180))[0])
            if ix.size > 0 :
                hex_path = matplotlib.path.Path(zip(hexX[ix],hexY[ix]))
                hex_patch = matplotlib.patches.PathPatch(hex_path, facecolor='none', edgecolor=color, lw=lw, fill=False)
                ax.add_patch(hex_patch)

            if ix_pos.size > 0:
                hex_path = matplotlib.path.Path(zip(hexX[ix_pos],hexY[ix_pos]))
                hex_patch = matplotlib.patches.PathPatch(hex_path, facecolor='none', edgecolor=color, lw=lw, fill=False)
                ax.add_patch(hex_patch)
            
            if ix_neg.size > 0:
                hex_path = matplotlib.path.Path(zip(hexX[ix_neg],hexY[ix_neg]))
                hex_patch = matplotlib.patches.PathPatch(hex_path, facecolor='none', edgecolor=color, lw=lw, fill=False)
                ax.add_patch(hex_patch)


    return ax

def plotLigoContours(x,y, vals, color="w", alpha = 1.0, lw=0.66, ls="solid", labels=False ) :
    import matplotlib
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata

    con_levels=5
    levels=[0.50, 0.90]
    #levels1=[0.50,]
    #levels2=[0.90,]
    print "\t\t contours at confidance levels 0.5, 0.9"

    xmin = x.min(); xmax = x.max()
    ymin = y.min(); ymax = y.max()
    
    coord = np.array(zip(x,y))
    xi=np.linspace(xmin, xmax, 500)
    yi=np.linspace(ymin, ymax, 500)
    xi,yi=np.meshgrid(xi,yi)

    zi = griddata(coord,vals,(xi,yi),method="cubic")
    #zi = griddata(coord,vals,(xi,yi),method="linear")

    #print "linestyle = ",ls
    ct= plt.contour(xi,yi,zi,con_levels,linewidths=lw, linestyles=ls,
        colors=color, levels=levels, alpha=alpha)
    #ct= plt.contour(xi,yi,zi,con_levels,linewidths=lw, linestyles=ls2,
    #    colors=color, levels=levels1, alpha=alpha)
    if labels :
        inline= True
        inline= False
        fontsize =10
        fontsize =14
        #plt.clabel(ct, levels, inline=inline, fontsize=fontsize)
        #Hack

def confidenceLevels(map) :
    map = np.array(map, copy=True)
    ix = np.argsort(map)[::-1]  ;# greatest to least
    newmap = np.cumsum(map[ix]/map.sum())
    inverse_ix = np.argsort(ix)
    newmap = newmap[inverse_ix]
    return newmap

#   ra, dec, ligo, maglim, prob, ha, x,y, hx,hy = readMaps(
def readMaps(mapDir, simNumber, slot) :
    import healpy as hp
    # get the maps for a reasonable slot
    name = os.path.join(mapDir, str(simNumber) + "-"+str(slot))
    print "\t reading ",name+"-ra.hp  & etc"
    raMap     =hp.read_map(name+"-ra.hp", verbose=False);
    decMap    =hp.read_map(name+"-dec.hp", verbose=False);
    haMap     =hp.read_map(name+"-ha.hp", verbose=False);
    xMap      =hp.read_map(name+"-x.hp", verbose=False);
    yMap      =hp.read_map(name+"-y.hp", verbose=False);
    hxMap     =hp.read_map(name+"-hx.hp", verbose=False);
    hyMap     =hp.read_map(name+"-hy.hp", verbose=False);
    ligoMap   =hp.read_map(name+"-map.hp", verbose=False);
    maglimMap =hp.read_map(name+"-maglim.hp", verbose=False);
    probMap   =hp.read_map(name+"-probMap.hp", verbose=False);
    haMap=haMap/(2*np.pi/360.)
    raMap=raMap/(2*np.pi/360.)
    decMap=decMap/(2*np.pi/360.)
    return raMap, decMap, ligoMap, maglimMap, probMap, \
        haMap, xMap, yMap, hxMap, hyMap


