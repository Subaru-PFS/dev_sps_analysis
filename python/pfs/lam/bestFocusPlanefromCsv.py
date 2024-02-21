from pfs.lam.imageAnalysis import neighbor_outlier_filter
from pfs.lam.detAnalysis import *
from pfs.lam.detFocusAnalysis import *
from pfs.lam.fileHandling import *
from pfs.lam.style import *
from pfs.lam.analysisPlot import *
from pfs.lam.misc import find_nearest

import lsst.daf.persistence as dafPersist

import glob
import time, sys, os
import argparse
from datetime import datetime

import matplotlib
matplotlib.use('Agg')




def main(experimentId, cam, rerun, criteria, basePath, drpPath, repo, roi_size, doBck, piston_index, filtered_waves):

    dat = datetime.now().strftime("%Y-%m-%dT%Hh%M")
    print(dat)

    filtered_waves = filtered_waves if filtered_waves is None else np.array(filtered_waves)

    
    repoRoot = f"{drpPath}/{repo}"
    extra = ""

    # define defaut parameters
    com = True  # Center Of Mass
    head = 0
    tail = 0
    criterias = [criteria]

    verbose = False
    doPrint = False
    
    doNeighborFilter = True
    thres_neighbor = -0.1 # default -0.1   
    
    
    doSave = True
    doSaveCsv = True
    
    
    files = []
    

    csvPath = os.path.join(basePath,"Exp"+str(experimentId)+"/"+rerun+"_roi"+str(roi_size)+"/doBck"+str(doBck)+"/"+extra)
    if not os.path.exists(csvPath):
        csvPath = os.path.join(basePath,f"sm{cam[1]}","Exp"+str(experimentId)+"/"+rerun+"/"+"roi"+str(roi_size)+"/doBck"+str(doBck)+"/"+extra)

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
    
    # Remove a wavelenght if needed 
    #
    if filtered_waves is not None:
        print(filtered_waves)
        if filtered_waves.size != 0:
            for wave in filtered_waves:
                filter_wave = waves[np.searchsorted(waves, wave)]
                print(f"{filter_wave:.2f} filtered")
                piston_filtered = piston_filtered[piston_filtered.wavelength != filter_wave]
    
    #
    # Fit
    #
    width = None
    #bounds = [(-200, 0.86, 150), (500, 1, 163)]
    if criteria == "EE5":
        bounds = [(-200, 0.86, 150), (500, 1, 163)]
        vmin = 0.75
        vmax = 1
        extra_title = f"_bounds_086_163_err_max_{piston_index}"
    else :
        bounds = None
        vmin = -1
        vmax = -1
        extra_title = f"_{piston_index}"

    
    thfoc_data = getAllBestFocus(piston_filtered.dropna(subset=['EE5']), criterias=['EE5'], index=piston_index, doPlot=True, savePlot=doSave, plot_path=csvPath, plot_title=f"Exp{experimentId}_thfocus_fit{extra_title}_{dat}",\
                             bounds=bounds)
    # save in CSV
    thfoc_data.to_csv(f"{csvPath}imquality_{cam}_Exp{experimentId}_doBck{str(doBck)}_piston_thfocfitdata{extra_title}_{dat}.csv")
    
    # plot Estimation of criteria EE5
    criteria_max = f"{criteria}_max"
    plot_title = f"{cam} Exp{experimentId} {criteria} thFocus estimation doBck{str(doBck)}"
    plot_file = f"{csvPath}{cam}_Exp{experimentId}_{criteria}_thFocusEstimation_doBck{str(doBck)}_{extra_title}"
    plotImageQualityScatter(thfoc_data, par=criteria_max, hist=criteria_max, savePlotFile=plot_file+f"_{dat}.csv", title=plot_title, doSave=doSave,\
                       vmin=vmin, vmax=vmax)
    plotImageQualityScatterFiberWave(thfoc_data, par=criteria_max, hist=criteria_max, savePlotFile=plot_file+f"wavefiber_{dat}.csv", title=plot_title, doSave=doSave,\
                       vmin=vmin, vmax=vmax)
    plotCumulativeImageQuality(thfoc_data, par=criteria_max, savePlotFile=plot_file+f"cumulative_{dat}.csv", title=plot_title+" cumulative", doSave=doSave,vmin=vmin, vmax=vmax )
    

    # Best Plane
    piston_plane = getBestPlane(thfoc_data[thfoc_data.criteria == "EE5"].dropna(), coords=["px", "py", "focus"], order=1, doPlot=True, exp=experimentId, plot_path=dataPath, plot_title=f"Focus3DPlane_{cam}_Exp{experimentId}_{dat}", savePlot=doSave)

    txtfile = f"SM1_{cam.upper()}_BestFocusPlane_Exp{experimentId}_doBck{str(doBck)}_{criteria}{extra_title}_{dat}.dat"

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
    
    parser = argparse.ArgumentParser(description="Cluster_GetImqual2csv.py argument parser")
    
    parser.add_argument("-c","--cam", type=str, help="camera to process like 'r3'", required=True)
    parser.add_argument("-r","--rerun", type=str, help="data rerun. could be 'ginga/detrend' or your rerun folder ", required=True)
    parser.add_argument("-e","--experimentId", type=int ,help="experimentId or visit_set_id")
    parser.add_argument("-b","--basePath", type=str, default="/data/drp/analysis",help="path where to store the results")
    parser.add_argument("--criteria", type=str, default="EE5",help="thfocus criteria, EE5 by default")
    parser.add_argument("--drpPath", type=str, default="/data/drp",help="main drp folder")
    parser.add_argument("--repo", type=str, default="sps",help="drp repository")
    parser.add_argument("--roi_size", type=int, default=24,help="roi_size in px used to calculate the total flux")
    parser.add_argument("--doBck", action="store_true" ,default=True,help="local bck substraction. default=True")
    parser.add_argument("--piston_index", type=str, default="motor1",help="index used for throughfocus analysis")
    parser.add_argument("--filtered_waves", type=float, nargs="+", default=[],help="waves to filtered")
    

    
    args = parser.parse_args()
    
    cam = args.cam
    rerun = args.rerun
    experimentId = args.experimentId
    criteria = args.criteria
    basePath = args.basePath
    drpPath = args.drpPath
    repo = args.repo
    roi_size = args.roi_size
    doBck = args.doBck
    piston_index = args.piston_index
    filtered_waves = args.filtered_waves
    filtered_waves = filtered_waves if filtered_waves is None else np.array(filtered_waves)


    
    main(experimentId, cam, rerun, criteria, basePath, drpPath, repo, roi_size, doBck, piston_index, filtered_waves)
    finish = time.time()
    elapsed = finish - start
    print(f"Time elapsed: {elapsed}")

    