def removeClosePeak(df, dist=80, doPlot=False):
    """
    Set a 'ignore' flag for points which closer to the distance <dist>
    input pandas dataframe
    output dataframe
    Use KDTree see:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query_ball_point.html#scipy.spatial.KDTree.query_ball_point
    """
    from scipy import spatial
    # h
    points = np.c_[df.px.ravel(), df.py.ravel()]
    tree = spatial.KDTree(points)
    
    badpoints = []
    for point in points:
        result = tree.query_ball_point(point, dist)
        # remove the point itself 
        if len(result) > 1 : badpoints.append(result)
    if doPlot:
        points = np.asarray(points)
        plt.plot(points[:,0], points[:,1], '.')
        for badpoint in badpoints:
            nearby_point = points[badpoint]
            plt.plot(nearby_point[:,0], nearby_point[:,1], 'x')
        plt.show()
    
    # set a flag "ignore" in the dataframe to identify the peak to be ignored

#   df["ignore"] = 0
    for result in badpoints:
        nearby_points = points[result]
        for nearby_point in nearby_points:
            df.loc[(df.px == nearby_point[0]) & (df.py == nearby_point[1]), "ignore"] = 1
    
#    return df.where(df.ignore<1).dropna()
    return df


def removeFluxPeak(df, fmin=2000, fmax=45000, doPlot=False):
    """
    Set a 'ignore' flag for points which a brigthness higher than <fmin> and lower than <fmax>
    input pandas dataframe
    output dataframe
    """
    if doPlot:
        ax= df.plot.scatter(x="px", y="py")
        df.where(df.brightness> fmax).plot.scatter(x="px", y="py", c="red", ax=ax)
    
    df.loc[(df.brightness> fmax) | (df.brightness< fmin), "ignore"] = 1
    
#    return df.where((df.brightness< fmax) & (df.brightness> fmin)).dropna()
    return df


def removeOnePeak(df, px=None, py=None, thresPix =3, peak=None, fiber=None, doPlot=False):
    """
    Set a 'ignore' flag for points which coordinate is (<px> -/+ <thresPix> ,<py> -/+ <thresPix> ) 
    input pandas dataframe
    output dataframe
    """
    df.loc[((py-thresPix<df.py) & (df.py< py+thresPix) & (px-thresPix<df.px) & (df.px< px+thresPix)) , "ignore"] = 1
    
#    return df.where((df.brightness< fmax) & (df.brightness> fmin)).dropna()
    return df

