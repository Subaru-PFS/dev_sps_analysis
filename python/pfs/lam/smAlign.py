import numpy as np

class arm():
    def __init__(self, arm):
        self.RfrontBell=372.5 #radius in mm
        self.VPHG_Set = 0.867 # arsec / µm
        self.arm = arm
        
    def getRollCu(self, df, extremFibers=(2,650), doPrint=False):
        # This shift is calculated using position of a center wavelength for the 2 extreme fibers (2 and 650)
        # TBO: for visible camera 
        # sign is given by doing differnce between 2 and 650
        # for nir it is between 650 and 2 
        i0 = 0
        i1 = 1
        if self.arm == "n":
            i0 = 1
            i1 = 0
        cu_roll_shift = df[df.fiber == extremFibers[i0]].y.values - df[df.fiber == extremFibers[i1]].y.values
        distance =  df[df.fiber == extremFibers[0]].x.values - df[df.fiber == extremFibers[1]].x.values

        # Calculation of shift
        roll_angle = math.atan(np.mean(cu_roll_shift/distance))
        cu_roll_shift = cu_roll_shift.mean()

        roll_shim = self.RfrontBell * math.tan(roll_angle)
        annot = "roll_shift %.1f px\n"%cu_roll_shift
        annot += "roll_angle %.1f arsec \n"%(math.degrees(roll_angle)*3600)
        annot += "roll_shim %.2f mm\n"%roll_shim

        if doPrint:
            print(distance)
            print(roll_angle)
            print(annot)
            print()

        if abs(math.degrees(roll_angle)*3600) < 50:
            if doPrint:
                print(f"{abs(math.degrees(roll_angle)*3600)} < 50 arsec : alignment is ok")
        return annot, roll_shim, roll_angle
    
    def plot_CU_roll(self, df, dataId, midWaves, fname=None, doSavePlot=False, site='LAM', doPrint=False):
    
        annot, roll_shim, roll_angle = self.getRollCu(df, doPrint=doPrint)

        fig, ax = plt.subplots(figsize=(12,8))
        plt.grid()
        plt.title(f"{site} - SM{dataId['spectrograph']} - {dataId['arm'].upper()}{dataId['spectrograph']} - CU Roll adjustment \n visitId = {dataId['visit']}")

        lns1 = ax.plot("x", "y", "*-", label=f"{midWaves[0]}nm", data=df[df.wavelength == midWaves[0]])
        secax = ax.secondary_xaxis(1, functions=(px2fiber, fiber2px))
        secax.set_xlabel('Fiber')
        axy = ax.twinx()
        lns2 = axy.plot("x", "y", "*-", label=f"{midWaves[1]}nm", data=df[df.wavelength == midWaves[1]], color="r")

        # added these three lines
        lns = lns1+lns2
        labs = [l.get_label() for l in lns]
        ax.legend(lns, labs, loc=0)

        #ax = df.groupby("wavelength")[["x","y"]].plot(kind="scatter")
        #ax.plot("x", "y", data=df.groupby("wavelength")) #  "*-", label=f"{midWave}nm")
        color=None
        if abs(math.degrees(roll_angle)*3600) <= 50 :
            color = "green"
        ax.annotate(annot,  xy=(.45,.45),xycoords="figure fraction", color=color)
        if doSavePlot:
        #fig.patch.set_alpha(0.5)
            plt.savefig(fname)
    #        plt.savefig(fname, transparent=True)


    
    
    def plot_vphg_roll(self, df, dataId, fibers, fname=None, doSavePlot=False, site='LAM'):

        fig, ax = plt.subplots(figsize=(12,8))
        plt.grid()

        plt.xlabel("X pixel")
        plt.ylabel("Y pixel")
        plt.title(f"{site} - SM{dataId['spectrograph']} - {dataId['arm'].upper()}{dataId['spectrograph']} - VPHG Roll adjustment \n visitId = {dataId['visit']}")


        x0 = df.x.values[0] 

        lns1 = ax.plot( df.x, df.y, 'o', label=f"{*fibers,}", color="black")

        # fit data 
        popt = np.polyfit( df.x, df.y, 1)
        print(popt)
        t = np.arange(df.x.min(), df.x.max()+1)
        lns2 = ax.plot(t,np.polyval(popt, t), 'b--', label="linear fit")

        annot = f"Shift { df.x.max() - df.x.min():.1f} px over {df.y.max() - df.y.min():.1f} px \n" 

        Roll_VPHG_shim_fit = np.rad2deg(1/popt[0])*3600*self.VPHG_Set
        annot += f"VPHG Roll Shim: {Roll_VPHG_shim_fit:.0f} um \n"

        print(np.rad2deg(popt[0])*3600*self.VPHG_Set)

        print(annot)

        plt.annotate(annot, xy=(.58,.2),xycoords="figure fraction")

        lns = lns1+lns2
        labs = [l.get_label() for l in lns]
        ax.legend(lns, labs, loc=0)


        if doSavePlot:
        #fig.patch.set_alpha(0.5)
            plt.savefig(fname)
    #        plt.savefig(fname, transparent=True)
        return Roll_VPHG_shim_fit