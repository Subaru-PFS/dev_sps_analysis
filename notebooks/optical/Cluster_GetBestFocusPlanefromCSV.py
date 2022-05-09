from pfs.lam.imageAnalysis import neighbor_outlier_filter

from pfs.lam.detAnalysis import *
from pfs.lam.detFocusAnalysis import *
from pfs.lam.fileHandling import *
from pfs.lam.style import *
from pfs.lam.analysisPlot import *
import glob
import lsst.daf.persistence as dafPersist

from pfs.lam.misc import find_nearest
import time, sys
import getopt

from datetime import datetime
import matplotlib
matplotlib.use('Agg')


datetime.now().strftime("%Y-%m-%dT%Hh%M")


def main(argv):
    visit = ''
    peak = ''
    basePath = os.path.join('/data/drp/analysis/sm3/')
    drpPath = "/data/drp"
    repo = "sps"
    repoRoot = f"{drpPath}/{repo}"
    experimentId = None
    rerun = "ginga"
    extra = ""

    try:
        opts, args = getopt.getopt(argv,"hc:b:e:r:",["cam=", "basePath=", "experimentId=","rerun="])
    except getopt.GetoptError:
        print('Cluster_GetBestFocusPlanefromCSV.py -b <basePath> -c <cam> -e <experimentId>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('GetIMqual2csv-standalone.py -v <visit> -p <peakfile> -c <cam>')
            sys.exit()
        elif opt in ("-b", "--basePath"):
            basePath = arg
        elif opt in ("-c", "--cam"):
            cam = arg
        elif opt in ("-e", "--experimentId"):
            experimentId = arg
        elif opt in ("-r", "--rerun"):
            rerun = arg

           


    #outpath = "output\\" if outpath is None else outpath

    # define defaut parameters
    roi_size = 30
    seek_size = 60

    com = True  # Center Of Mass
    doBck = True
    head = 0
    tail = 0
    criteria = 'EE5'
    criterias = [criteria]
    piston_index = 'motor1'

    verbose = False
    doPrint = False
    
    doNeighborFilter = True
    thres_neighbor = -0.1 # default -0.1   
    
    
    doSave = True
    doSaveCsv = True
    
    
    files = []
    

    csvPath = basePath+"Exp"+str(experimentId)+"/"+rerun+"_roi"+str(roi_size)+"/doBck"+str(doBck)+"/"+extra
    if not os.path.exists(csvPath):
        csvPath = basePath+"Exp"+str(experimentId)+"/"+rerun+"/"+"roi"+str(roi_size)+"/doBck"+str(doBck)+"/"+extra

    dataPath = csvPath
    print(dataPath)

    searchFile = f"{csvPath}Imquality_{cam}_Exp{experimentId}*"
    print("\n"+searchFile+"\n")
    files.extend(glob.glob(searchFile))

    if verbose:
        print(*files, sep="\n")

    piston_imdata = pd.concat(map(pd.read_csv,files)).reset_index().sort_values(by="motor1")

    
    # hack debut sm2
    print(f"Wavelength used: {piston_imdata.wavelength.unique()}")
    minPos = piston_imdata[['motor1','motor2','motor3']]
    minPos = minPos - minPos.min()
    piston_imdata['relPos'] = minPos['motor1']

    
    waves = piston_imdata.wavelength.sort_values().unique()
    fibers = piston_imdata.fiber.sort_values().unique()
    print(fibers, waves)
    
    
    # filtering the data
    
    
    if doNeighborFilter :
        tot = []
        for  group_name, series in piston_imdata.groupby(['wavelength','fiber']):
            series = neighbor_outlier_filter(series, "EE5", thres_neighbor, absolute=True)
            tot.append(series)
            piston_filtered = pd.concat(tot)
    
    # 2nd Moment calculation from SEP
    if "sep_x2" in piston_filtered.columns:
        piston_filtered["sep_2ndM"]= piston_filtered.apply(lambda x: np.mean([x["sep_x2"],x["sep_y2"]]) , axis=1)
        piston_filtered["sep_radius_x"]= np.sqrt(piston_filtered["sep_x2"])
        piston_filtered["sep_radius_y"]= np.sqrt(piston_filtered["sep_y2"])
        piston_filtered["sep_radius_xy"]= piston_filtered.apply(lambda x: np.mean([x["sep_radius_x"],x["sep_radius_y"]]) , axis=1)
    
    piston_filtered = piston_filtered[piston_filtered.EE5_nbh_flag]
    
    # Remove a wavelenght if need 
    # faudrait faire ça mieux...
#    if cam[0] == "r":
#        filter_wave = waves[-1]
#        print(f"{filter_wave:.2f} filtered")
#        piston_filtered = piston_filtered[piston_filtered.wavelength != filter_wave]
    
    #
    # Fit
    #
    width = None
    bounds = [(-200, 0.86, 150), (500, 1, 163)]
    extra_title = f"_bounds_086_163_err_max_{piston_index}"
    
    thfoc_data = getAllBestFocus(piston_filtered.dropna(subset=['EE5']), criterias=['EE5'], index=piston_index, doPlot=True, savePlot=doSave, plot_path=csvPath, plot_title=f"Exp{experimentId}_thfocus_fit{extra_title}",\
                             bounds=bounds)
    # save in CSV
    thfoc_data.to_csv(f"{csvPath}imquality_{cam}_Exp{experimentId}_doBck{str(doBck)}_piston_thfocfitdata{extra_title}.csv")

    # Best Plane
    piston_plane = getBestPlane(thfoc_data[thfoc_data.criteria == "EE5"].dropna(), coords=["px", "py", "focus"], order=1, doPlot=True, exp=experimentId, plot_path=dataPath, plot_title=f"Focus3DPlane_{cam}_Exp{experimentId}", savePlot=doSave)

    txtfile = f"SM1_{cam.upper()}_BestFocusPlane_Exp{experimentId}_doBck{str(doBck)}_{criteria}{extra_title}.dat"

    invMat, invMatName = getFocusInvMat(cam)
    foc = findMotorPos(piston_plane, cam=cam)

    (tip, tilt, tip_mic, tilt_mic), defocus = getPlaneInfo(piston_plane, doPrint=True)
    fca_focus = piston_filtered.fcaFocus.mean()
    txt = f"{datetime.now()} \n"
    txt += f"{csvPath}\n"
    txt += f"ExpId {experimentId}\n"
    txt += f"{invMatName}\n"
    txt += f"{foc}\n"
    txt += f"xcu_{cam} motors moveCcd a={foc[0]:.2f} b={foc[1]:.2f} c={foc[2]:.2f} microns abs\n"
    txt += f"Tip {tip:.3e} rad => {tip_mic:.1f} microns\n"
    txt += f"Tilt {tilt:.3e} rad => {tilt_mic:.1f} microns\n"
    txt += f"{defocus:.1f} microns\n"

    print(txt)

    # debug B2 
    #dataPath = "/drp/analysis/sm2/debugB2"

    if doSave :
        text_file = open(csvPath+txtfile, "w")
        text_file.write(txt)
        text_file.close()

    if "detBoxTemp" in piston_filtered.columns:
        detBoxTemp_mean = piston_filtered.detBoxTemp.mean()
    else:
        detBoxTemp_mean = np.nan
    if "ccdTemp" in piston_filtered.columns:
        ccdTemp_mean = piston_filtered.ccdTemp.mean()
    else:
        ccdTemp_mean = np.nan          


    plot_prefix = f"Focus_Piston_plots_doBck{str(doBck)}_{criteria}{extra_title}"
    title_suffix = f"a={foc[0]:.2f} b={foc[1]:.2f} c={foc[2]:.2f} microns"
    title_suffix += f"- DetBox {detBoxTemp_mean:.1f}K - ccdTemp {ccdTemp_mean:.1f}K - FCA Focus {fca_focus:.1f}mm \n"
    title_suffix += f"Tip {tip:.3e} rad => {tip_mic:.1f}µm -- Tilt {tilt:.3e} rad => {tilt_mic:.1f}µm -- Center focus {defocus:.1f}µm"
    plot_groups(piston_filtered.sort_values("motor1")[piston_filtered.EE5_nbh_flag], experimentId, dataPath, plot_prefix=plot_prefix, title_suffix=title_suffix,\
                col="fiber", hue="wavelength", criteria=criteria, doSave=doSave, verbose=True)
    plot_groups(piston_filtered.sort_values("motor1")[piston_filtered.EE5_nbh_flag], experimentId, dataPath, plot_prefix=plot_prefix, title_suffix=title_suffix,\
                col="wavelength", hue="fiber", criteria=criteria, doSave=doSave)

    
    

if __name__ == "__main__":
    start = time.time()
    main(sys.argv[1:])
    finish = time.time()
    elapsed = finish - start
    print(f"Time elapsed: {elapsed}")

    
