import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime


def plotOnePeak(image, cx,cy, roi_size=30, doBck=False, nRows=5, vmin=None, vmax=None, verbose=False):
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
    plt.show()



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
    plt.show()
    

    
def plot_groups(piston_imdata, experimentId, plot_path, plot_prefix="Focus_Piston_plots", title_suffix=None, col="fiber", hue="wavelength", criteria="EE5", doSave=False, verbose=False) :
    grid = sns.FacetGrid(piston_imdata, col=col, hue=hue,
#                         col_wrap=4, height=3, legend_out=True)
                         col_wrap=3, height=3, legend_out=True)

    grid.map(plt.plot, "motor1", criteria, marker="+")
    grid.fig.tight_layout(w_pad=1)
    if criteria == "EE5":
        grid.set(ylim=(0, 1))
    grid.add_legend()
    plt.subplots_adjust(top=0.75)
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
    plt.show()