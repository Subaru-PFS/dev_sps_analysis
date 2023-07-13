import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec

from pfs.lam.imageAnalysis import selectRoi, getRois

def plotOnePeak(image, cx,cy, roi_size=30, doBck=False, nRows=5, vmin=None, vmax=None, verbose=False,title=None):
    indx = cy
    indy = cx

    if type(image) is str:
        hdulist = fits.open(image, "readonly")
        image = hdulist[1].data
    

    data = np.copy(image)
    outer_data, inner_data = getRois(data, cx, cy, inner_size=5, outer_size=roi_size, doBck=doBck, nRows=nRows)    
    m, s = np.mean(outer_data), np.std(outer_data)
    if verbose:
        print(f"mean: {m}")
        print(f"std: {s}")
    vmin = vmin if vmin is not None else m-s
    vmax = vmax if vmax is not None else m+s
    
    fig, (ax, ax2) = plt.subplots(ncols=2, constrained_layout=True, figsize=(18,9))
    #ax = plt.subplot(111)
    im = ax.imshow(outer_data,interpolation="none", origin="lower", vmin=vmin, vmax=vmax)
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im, ax=ax)

    #m = np.log10(np.where(outer_data>0, outer_data, -1))

    #ax2 = plt.subplot(121)
    im2 = ax2.imshow(outer_data,interpolation="none", origin="lower",cmap="gray", norm=LogNorm())
    
    fig.colorbar(im2, ax=ax2)
#    plt.show()
    if title is not None:
        plt.title(title)

def plotRoiPeak(image, peak_list, roi_size=20, raw=False, scale=True, verbose=False, savePlotFile=None, doSave=False):
    if type(image) is str:     
        hdulist = fits.open(image, "readonly")
        image = hdulist[1].data
    
    plist = pd.read_csv(peak_list) if type(peak_list) is str else peak_list
    # 2023/02/03
    # drop duplicates
    plist = plist.drop_duplicates(subset=["fiber", "wavelength"])
#    plist = plist.sort_values(["X","Y"], ascending=[True,False])
    # Use X,Y when it the list of peak, px,py otherwise
    ind_x = "px"
    ind_y = "py"
    if raw :
        ind_x = "X"
        ind_y = "Y"        
    
    nbfiber = len(plist.fiber.unique())
    nbwave = len(plist.wavelength.unique())
    # have a list of fiber from 0 to nbfiber
    listwavestoplot = pd.DataFrame(plist.wavelength.unique(), columns=["wavelength"])

    listfiberstoplot = pd.DataFrame(plist.fiber.unique(), columns=["fiber"])
    
    print("#Fiber= %d and #wavelength= %d"%(nbfiber, nbwave))
    f, axarr = plt.subplots(nbwave, nbfiber,  sharex='col', sharey='row',figsize=(12,8))
#    print(axarr.shape)
    vmin=None
    vmax=None
    
    for (wave, fiber), group in plist.groupby(['wavelength','fiber']):
        k = listwavestoplot[listwavestoplot.wavelength == wave].index.tolist()[0]
        i = listfiberstoplot[listfiberstoplot.fiber == fiber].index.tolist()[0]
        if verbose:
            print(k,i)
            print(f"px {group[ind_x]}    py: {group[ind_y]}")
        #cut_data = image[int(indx-roi_size/2):int(indx+roi_size/2), int(indy-roi_size/2):int(indy+roi_size/2)]
        cut_data = selectRoi(image, group[ind_x], group[ind_y], roi_size=roi_size)
        if nbwave == 1 and nbfiber == 1:
            axarr.set_title(f"{str(fiber)}, {str(wave)}")
            axarr.imshow(cut_data,interpolation="none", origin="lower", vmin=vmin, vmax=vmax)
        else:
            #axarr[nbwave -1 -k, nbfiber - i -1].set_title(f"{fiber}, {wave:.2f}")
            #axarr[nbwave -1 -k, nbfiber - i -1].label_outer()
            axarr[nbwave -1 -k, nbfiber - i -1].imshow(cut_data,interpolation="none", origin="lower", vmin=vmin, vmax=vmax)
            #axarr[nbwave -1 -k, nbfiber - i -1].set_ylabel(f"{wave:.2f}")
            #axarr[nbwave -1 -k, nbfiber - i -1].set_xlabel(fiber)
            axarr[nbwave -1 -k, nbfiber - i -1].grid(False)

                
    f.subplots_adjust(hspace=0.5,wspace=0.5)

    for ax, wave in zip(axarr[:,0], listwavestoplot.sort_index(ascending=False).wavelength.values) :
            ax.set_ylabel(f"{wave:.2f}", rotation='horizontal', ha='right', fontsize=20)
#            ax.set_xlabel(ax.get_xlabel(), rotation='vertical', fontsize=20)
            ax.set_yticklabels('')
            ax.set_xticklabels('')
            ax.set_frame_on(False)
    for ax, fiber in zip(axarr[-1,:], listfiberstoplot.sort_index(ascending=False).fiber.values):
#            ax.set_ylabel(ax.get_ylabel(), rotation='horizontal', ha='right', fontsize=20)
            ax.set_xlabel(fiber, rotation='vertical', fontsize=20)
            ax.set_yticklabels('')
            ax.set_xticklabels('')
            ax.set_frame_on(False)

    plt.gcf().set_facecolor('w')
    if doSave:
        f.patch.set_alpha(0.5)
        plt.savefig(savePlotFile+f"_roi_all.png")
                  
#    plt.show()
    

def plotImageQuality(dframe, vmin=-1,vmax=-1, par="EE3", hist=None, filelist=None, com=False, doSave=False, imgPath="/home/pfs/shared/Pictures/" ):
    
    # select peak center 
    # default x , y are objx and objy, but if it is the center of Mass it is oid_x and oid_y
#    x = dframe["objx"]
#    y = dframe["objy"]
#    if com :
#        x = dframe["oid_x"]
#        y = dframe["oid_y"]

    ## should now be px, py which are affected during calculation according com or not
    x = dframe["px"]
    y = dframe["py"]    
    
    z = dframe[par]
#    stat = 'ECE5' if par == 'ECE5' else 'EE5'
    stat = par
    xs = dframe[dframe[stat]>0.90].px
    ys = dframe[dframe[stat]>0.90].py
    zs = dframe[dframe[stat]>0.90][stat]

    print("%.f %% %s peaks >0.90"%(100*len(zs)/len(z),stat))
    # Set up a regular grid of interpolation points
    xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
    xi, yi = np.meshgrid(xi, yi)

    # Interpolate
    rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
    zi = rbf(xi, yi)

    if vmin == -1 :
        vmin = z.min()
    if vmax == -1 :
        vmax = z.max()
    if hist is not None:
        #fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
        fig = plt.figure(figsize=(16, 8))
        gs = gridspec.GridSpec(1, 3,
 #                      width_ratios=[3,1],
 #                      height_ratios=[1,1]
                       )
        ax1 = plt.subplot(gs[0,:2])
        ax2 = plt.subplot(gs[0,2])
        im = ax1.imshow(zi, vmin=vmin, vmax=vmax, origin='lower',
                   extent=[x.min(), x.max(), y.min(), y.max()])
        ax1.scatter(x, y, c=z, vmin=vmin, vmax=vmax)
        ax1.set_title(par)
        ax1.set_xlim([0,4095])
        ax1.set_ylim([0,4175])
        fig.colorbar(im, ax=ax1, shrink=1)
        dframe[hist].plot.hist(ax=ax2, bins=20)
        ax2.set_title("Histogram")
        if filelist is not(None):
            fig.suptitle(getFileInfo(filelist))
        
    else : 
        fig = plt.figure(figsize=(8, 8))

        plt.imshow(zi, vmin=vmin, vmax=vmax, origin='lower',
                   extent=[x.min(), x.max(), y.min(), y.max()])
        plt.scatter(x, y, c=z, vmin=vmin, vmax=vmax)
        
        plt.colorbar(shrink=0.65)
        plt.scatter(xs, ys, c="r", marker="x")

        plt.xlim(0,4095)
        plt.ylim(0,4175)
        if filelist is not(None):
            plt.title(getFileInfo(filelist))
          
    if doSave:
        fig.patch.set_alpha(0.5)
        plt.savefig(os.path.join(imgPath,\
            f"SM1_ImQual_{filelist.cam[0]}_EXP{filelist.experimentId[0]}_{par}_{filelist.obsdate[0]}.png"))  
#        plt.show()


def plotImageQualityScatter(dframe, par="EE3", vmin=-1,vmax=-1, hist=None, savePlotFile=None, com=False, doSave=False, title=None ):
    # select peak center 
    # default x , y are objx and objy, but if it is the center of Mass it is oid_x and oid_y
#    x = dframe["objx"]
#    y = dframe["objy"]
#    if com :
#        x = dframe["oid_x"]
#        y = dframe["oid_y"]

    ## should now be px, py which are affected during calculation according com or not
    x = dframe["px"]
    y = dframe["py"]   
    
    z = dframe[par]
    
    #    stat = 'ECE5' if par == 'ECE5' else 'EE5'
    val = 0.5 if par == "EE3" else 0.9

    stat = par
    xs = dframe[dframe[stat]>val].px
    ys = dframe[dframe[stat]>val].py
    zs = dframe[dframe[stat]>val][stat]

    statEE = f"{100*len(zs)/len(z):.1f}% of peak have a {par} > {val}"
        
    if vmin == -1 :
        vmin = z.min()
    if vmax == -1 :
        vmax = z.max()
    fact = 100
    if (par == "brightness") |(par == "sep_brightness"):
        fact = 1./100
    
    if hist is not None:
        #fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
        fig = plt.figure(figsize=(16, 8))
        gs = gridspec.GridSpec(1, 3,
    #                      width_ratios=[3,1],
    #                      height_ratios=[1,1]
                       )
        ax1 = plt.subplot(gs[0,:2])
        ax2 = plt.subplot(gs[0,2])
        im = ax1.scatter(x, y, c=z, s= z*fact, vmin=vmin, vmax=vmax)
        ax1.set_title(par)
        ax1.set_xlim([0,4095])
        ax1.set_ylim([0,4175])
        plt.colorbar(im,ax=ax1,shrink=1)

        dframe[hist].plot.hist(ax=ax2, bins=20)
        val = 0.5 if par == "EE3" else 0.9
        ax2.axvline(x=val, color='k')
        ax2.set_title("Histogram")
        ax2.set_xlim(vmin,vmax)
        if title is not(None):
            fig.suptitle(title+"\n"+statEE)
    else:
        fig = plt.figure(figsize=(10, 8))

        plt.scatter(x, y, c= z, s= z*fact, vmin=vmin, vmax=vmax)
        plt.colorbar(shrink=1)
        plt.xlim(0,4095)
        plt.ylim(0,4175)
        plt.xlabel('X')
        plt.ylabel('Y')

        plt.show()
        if title is not(None):
            fig.suptitle(title+"\n"+statEE)
          
    if doSave:
        fig.patch.set_alpha(0.5)
        plt.savefig(savePlotFile+f"_{par}.png", bbox_inches = "tight" )
#        plt.show()

def plotImageQualityScatterFiberWave(dframe, par="EE3", vmin=-1,vmax=-1, hist=None, savePlotFile=None, com=False, doSave=False, title=None, waveband=None ):
    '''
    
    '''
    markers = {"Ne": "^", "Ar": "o", "HgAr" : "D", "Kr": "s", "Xe": "h" }

    x = dframe["fiber"]
    y = dframe["wavelength"]
    
    z = dframe[par]
    
    #    stat = 'ECE5' if par == 'ECE5' else 'EE5'
    val = 0.5 if par == "EE3" else 0.9

    stat = par
    xs = dframe[dframe[stat]>val].px
    ys = dframe[dframe[stat]>val].py
    zs = dframe[dframe[stat]>val][stat]

    statEE = f"{100*len(zs)/len(z):.1f}% of peak have a {par} > {val}"
        
    if vmin == -1 :
        vmin = z.min()
    if vmax == -1 :
        vmax = z.max()
    fact = 100
    if (par == "brightness") |(par == "sep_brightness"):
        fact = 1./100
    
    if hist is not None:
        #fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
        fig = plt.figure(figsize=(16, 8))
        gs = gridspec.GridSpec(1, 3,
    #                      width_ratios=[3,1],
    #                      height_ratios=[1,1]
                       )
        ax1 = plt.subplot(gs[0,:2])
        ax2 = plt.subplot(gs[0,2])
        for name, group in dframe.groupby("lamp"):
            group = group.copy()
            m = markers.get(name)
#            im = ax1.scatter(x, y, c=z, s= z*fact, vmin=vmin, vmax=vmax)
            im = ax1.scatter(
                x=group["fiber"],
                y=group["wavelength"],
                c=group[par],
                s= group[par]*fact,
                vmin=vmin,
                vmax=vmax,
                marker=m,
                label=name,
            )
        ax1.invert_xaxis()
        ax1.set_title(par)
        ax1.set_xlabel('Fiber number')
        ax1.set_ylabel('wavelength (nm)')
        ax1.legend()
#        ax1.set_xlim([0,4095])
#        ax1.set_ylim([0,4175])
#        ax1.set_xlim([0,651])
#        ax1.set_ylim([630,970])
        plt.colorbar(im,ax=ax1,shrink=1)
        if waveband is not None:
            ax1.axhline(y=waveband[0], color='k', ls="--")
            ax1.axhline(y=waveband[1], color='k', ls="--")
            ax1.annotate(f"{waveband[0]}nm", (339, waveband[0]+1), xycoords="data")
            ax1.annotate(f"{waveband[1]}nm", (339, waveband[1]+1), xycoords="data")

        dframe[hist].plot.hist(ax=ax2, bins=20)
        val = 0.5 if par == "EE3" else 0.9
        ax2.axvline(x=val, color='k')
        ax2.set_title("Histogram")
        ax2.set_xlim(vmin,vmax)
        ax2.set_ylabel('count')
        ax2.set_xlabel(f'{par}')


        if title is not(None):
            if "EE" in par:
                fig.suptitle(title+"\n"+statEE)
            else:
                fig.suptitle(title)

    else:
        fig = plt.figure(figsize=(10, 8))

#        plt.scatter(x, y, c= z, s= z*fact, vmin=vmin, vmax=vmax)
        for name, group in dframe.groupby("lamp"):
            group = group.copy()
            m = markers.get(name)
#            im = ax1.scatter(x, y, c=z, s= z*fact, vmin=vmin, vmax=vmax)
            im = plt.scatter(
                x=group["fiber"],
                y=group["wavelength"],
                c=group[par],
                s= group[par]*fact,
                vmin=vmin,
                vmax=vmax,
                marker=m,
                label=name,
            )
        plt.legend()
        plt.gca().invert_xaxis()

        plt.colorbar(shrink=1)
        #plt.xlim(0,4095)
        #plt.ylim(0,4175)
        plt.xlabel('Fiber number')
        plt.ylabel('wavelength (nm)')

        plt.show()
        if title is not(None):
            if "EE" in par:
                fig.suptitle(title+"\n"+statEE)
            else:
                fig.suptitle(title)
          
    if doSave:
        fig.patch.set_alpha(0.5)
        plt.savefig(savePlotFile+f"_{par}.png", bbox_inches = "tight" )
#        plt.show()


def plotCumulativeImageQuality(dframe, par="EE5", savePlotFile=None, title=None, doSave=False, vmin=-1,vmax=-1 ):
    fig = plt.figure(figsize=(8, 8))

    criteria_values = dframe[par].values

    values, base = np.histogram(criteria_values, bins=40)
    #evaluate the cumulative
    cumulative = np.cumsum(values)
    # plot the cumulative function
    #plt.plot(base[:-1], 100*(cumulative/len(best.EE5.values)), c='blue')
    #plot the survival function
    plt.plot(base[:-1], 100*(len(criteria_values)-cumulative)/len(criteria_values), c='green')
    
    if vmin == -1 :
        vmin = base[:-1].min()
    if vmax == -1 :
        vmax = base[:-1].max()
    vl = 0.9
    if par == "EE3":
        vl = 0.5
    

    plt.xlim(vmin,vmax)
  
    plt.ylabel("% of peak")
    plt.xlabel("EE5")
    plt.axhline(y=95, c="r")
    plt.axvline(x=vl, c="r")
    plt.title(title)
    plt.gcf().set_facecolor('w')
    plt.annotate("req",(vl,95), color ="r" )
    if doSave:
        fig.patch.set_alpha(0.5)
        plt.savefig(savePlotFile+".png", bbox_inches = "tight" )
#    plt.show()
    
    
    


def plot_one_group(piston_imdata, wave, fiber, experimentId, plot_path, criteria="EE5", doSave=False) :
    group = piston_imdata.groupby(['wavelength','fiber']).get_group((wave,fiber))
    
    ax = group.plot.scatter("motor1", "EE5",\
                        title=f"Exp{experimentId} - fiber {fiber} - wave {wave}",\
                       label="EE5")
#    ax2 = ax.twinx()
#    group.plot.scatter("motor1", "sep_x2", ax=ax,color="r", label="sep_x2" )
#    group.plot.scatter("motor1", "sep_y2", ax=ax2,color="g", label="sep_y2" )
#    group.plot.scatter("motor1", "sep_ECE5", ax=ax2,color="y", label="sep ECE5" )

    #fig.patch.set_alpha(0.5)
    dat = datetime.now().isoformat(timespec='minutes') 
    if doSave : 
        plt.savefig(plot_path+f"Focus_fit_fiber{fiber}_wave{wave}_Exp{experimentId}_{dat}.png")
#    plt.show()
    

    
def plot_groups(piston_imdata, experimentId, plot_path, plot_prefix="Focus_Piston_plots", title_suffix=None, col="fiber", hue="wavelength", criteria="EE5", index="motor1", doSave=False, verbose=False) :
    grid = sns.FacetGrid(piston_imdata, col=col, hue=hue,
#                         col_wrap=4, height=3, legend_out=True)
                         col_wrap=3, height=3, legend_out=True)

    grid.map(plt.plot, index, criteria, marker="+")
    grid.fig.tight_layout(w_pad=1)
    if criteria == "EE5":
        grid.set(ylim=(0, 1))
    grid.add_legend()
    plt.subplots_adjust(top=0.80)
    if title_suffix is None:
        title = f"{plot_prefix} - Exp{experimentId}"
    else: 
        title = f"{plot_prefix} - Exp{experimentId}\n{title_suffix}"
    grid.fig.suptitle(title)
#    dat = datetime.now().isoformat(timespec='minutes') 
    dat = datetime.now().strftime("%Y-%m-%dT%Hh%M")
    if doSave:
        if verbose:
            print(f"write png file: \n {plot_path}{plot_prefix}_{col}_{hue}_Exp{experimentId}_{dat}.png")
        plt.savefig(plot_path+f"{plot_prefix}_{col}_{hue}_Exp{experimentId}_{dat}.png", transparent=True)
#    plt.show()


def plotThoughFocusData(piston_data, index="relPos", criterias=["EE5", "EE3", "2ndM"], head=0, tail=0, \
                   savePlot=False, plot_path="", plot_title=None, ylog=False, ylim=(None,None) ):
    """
    Plot Through focus data 
    piston_data: dataframe
    index: "x axis for the plot"
    criterias : "y axis"
    """
    piston = piston_data[head:piston_data.count()[0]-tail]

    grouped = piston.groupby(['wavelength','fiber'])


    ncols = len(piston.fiber.unique())
    nrows = len(piston.wavelength.unique())

    newx = np.linspace(np.min(piston[index].values), np.max(piston[index].values), 100)
    for criteria in criterias:
        if criteria in piston.columns:
            fig, axs = plt.subplots(nrows,ncols, figsize=(12,8), sharey='row', sharex=True)
            visit_info = f"{piston.visit.values.min()} to {piston.visit.values.max()}"
            if "detBoxTemp" in piston.columns:
                detBoxTemp_mean = piston.detBoxTemp.mean()
            else:
                detBoxTemp_mean = np.nan
            if "ccdTemp" in piston.columns:
                ccdTemp_mean = piston.ccdTemp.mean()
            else:
                ccdTemp_mean = np.nan                
            temp_info = f"detBox: {detBoxTemp_mean:.1f}K  ccd: {ccdTemp_mean:.1f}K"
            cam_info = piston.cam.unique()[0]
            #date_info = piston.obsdate[0].split('T')[0]
            date_info =""
            fca_focus = piston.fcaFocus.mean()
            fig.suptitle(f"{cam_info.upper()} ExpId {str(int(piston.experimentId.unique()[0]))} - {visit_info} - {criteria} - {temp_info} - FCA Focus {fca_focus:.1f}mm - {date_info}")
            plt.subplots_adjust(top=0.93)
            for (name, df), ax in zip(grouped, axs.flat):
                ax.set_title(f"{name[0]:.2f}, {name[1]}")
                df.plot.scatter(x=index,y=criteria, ax=ax)
                if criteria == "EE5" or criteria == "EE3":
                    ax.set_ylim(0,1)
                else :
                    bot, top = ylim                    
                    if ylog:
                        ax.set_yscale('log')
                    else :
                        bot = 0 if bot is None else bot
                    ax.set_ylim(bottom=bot, top=top)

            if savePlot:
                if plot_title is None:
                    plot_title = f"{cam_info.upper()}_ExpId_{str(int(piston.experimentId.unique()[0]))}_{criteria}_thFocusPlot{date_info}.png"
                plt.savefig(plot_path+plot_title)


def plotPeaksBrightness(df, doSave=False, plot_title=None,savePlotFile=None):
    """
    plot peak brightness by waves and fibers
    df: dataframe coming from ImageQualityToCsv (need to have sep_brightness wavelength and fiber)
    """
    fig = plt.figure(figsize=(10,12))
    plt.subplot(2, 1, 1)
    df.set_index("fiber").groupby(["wavelength"]).sep_brightness.plot(subplots=False, stacked=True,legend=True, logy=True, style="--*",  sharex=True)
    plt.legend(bbox_to_anchor=(1.0, 1.0))
    plt.ylabel("ADU")

    plt.subplot(2, 1, 2)
    df.set_index("fiber").groupby(["wavelength"]).sep_brightness.plot(subplots=False, stacked=True,legend=True, style="--*",  sharex=True)
    plt.legend(bbox_to_anchor=(1.0, 1.0))
    plt.ylabel("ADU")

    plt.suptitle(plot_title)
    plt.gcf().set_facecolor('w')

    if doSave:
        fig.patch.set_alpha(0.5)
        plt.savefig(savePlotFile+".png", bbox_inches='tight')
