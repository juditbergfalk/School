import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.animation as ani
from matplotlib.image import imread
import numpy as np

# Generate random diameters with given probability distribution
def SFD(x, slope, q):
    """Size-frequency distribution"""
    return slope * x**(-q)

def RandomDiameters(L1, L2, SFD=SFD, size=1, slope=1, q=3, control=1e6):
    """Generating random diameters in km with given probability distribution.
    
       Parameters
       ----------
           L1 : float
               Lower limit of the impactor diameter in km
           L2 : float
               Upper limit of the impactor diameter in km
           SFD : function
               Probability distribution function with the option to change the slope and y-intercept
               Format: slope * x^(-q)
               Default is the function SFD
           slope : float
               Slope of the probability function
               Default is 1
           q : float
               y-intercept of the probability function
               q >= 2 for a reasonable distribution
               Default is 3
           size : float
               Number of impactors
               Default is 1
           control : float
               Number of max calculations to avoid computer crash
               Default is 10^6
           
       Returns
       -------
           Random diameters in km: np.array
       
       Example
       -------
           Generate 10 impactors between 10 km and 100 km diameter with given SFD:
           >> print(RandomDiameters(10,100,size=10))
           >> [40.49 14.8  13.21 90.74 22.48 49.47 22.85 10.83 10.72 11.12]
    """
       
    N = []                                 # list of random numbers
    count = 0                              # counter
    assert q > 2  
    while len(N) < size and count < control:
        x = np.random.uniform(low=L1,high=L2)
        x = np.round(x,2)                  # round to two decimals
        function = SFD(x,slope,q)      # define function
        if np.random.uniform(0,1) <= function:
            N += [x]
        count += 1
    N = np.array(N)
    return N

# Calculate crater sizes 
def CraterCalculator(D, v_i=10000, g=9.81, rho_p=3000, rho_t=5500):
    """Calculating the size of the impactors from the size of a crater.
       
       Parameters
       ----------
           D : float
               Diameters of craters in km
           v_i : float
               Impact velocity in m/s
               Default is 10 km/s
           g : float
               Gravitational constant in m/s^2
               Default is 9.81 m/s^2
           rho_p : float
               Mass density of impactor in kg/m^3
               Default is 3000 kg/m^3
           rho_t : float
               Mass density of surface in kg/m^3
               Default is 5500 kg/m^3
           
       Returns
       -------
           Diameter of impactor in km: float
       
       Example
       -------
           What diameter of impactor makes 10 km craters?
           >> print(CraterCalculator(10))
           >> 1.22
    """
    D = D*1000           # convert km to m
    L = (g/v_i**2)**(1/3) * (rho_t/rho_p)**(1/3) * D**(4/3)
    L = L/1000           # convert m to km
    L = np.round(L,2)    # round to two decimals
    return L

# Calculate random crater sizes
def CraterCreator(low, high, SFD=SFD, size=1, 
                  v_i=10000, g=9.81,
                  rho_p=3000, rho_t=5500,
                  slope=1, q=3):
    """Creating craters from random diameter impactors, given lower and upper limit of crater sizes.
       
       Parameters
       ----------
           low : float
               Lower limit for crater diameter in km
           high : float
               Upper limit for crater diameter in km
           SFD : function
               Probability distribution function with the option to change the slope and y-intercept
               Format: slope * x^(-q)
               Default is the function SFD
           size : float
               Number of impactors
               Default is 1
           v_i : float
               Impact velocity in m/s
               Default is 10 km/s
           g : float
               Gravitational constant in m/s^2
               Default is 9.81 m/s^2
           rho_p : float
               Mass density of impactor in kg/m^3
               Default is 3000 kg/m^3
           rho_t : float
               Mass density of surface in kg/m^3
               Default is 5500 kg/m^3
           slope : float
               Slope of the probability function
               Default is 1
           q : float
               y-intercept of the probability function
               q >= 2 for a reasonable distribution
               Default is 3
           
       Returns
       -------
           List of crater diameters in km : np.array
       
       Example
       -------
           Generate 10 craters between 10 km and 100 km diameters.
           >> print(CraterCreator(10,100,size=10))
           >> [11.94 15.12 17.63 11.05 10.27 10.81 10.09 27.74 10.15 13.36]

    """
    # Finding the lower limit of impactors from given lower limit of crater size
    L1 = CraterCalculator(low,v_i=v_i,g=g,rho_p=rho_p,rho_t=rho_t)
    # Finding the upper limit of impactors from given lower limit of crater size
    L2 = CraterCalculator(high,v_i=v_i,g=g,rho_p=rho_p,rho_t=rho_t)
    
    # Generate random diameter impactors
    randomimpactors = RandomDiameters(L1,L2,SFD=SFD,size=size,slope=slope,q=q)
    # Convert km to m
    randomimpactors = randomimpactors*1000
    
    # Calculate crater diameters
    D = (v_i**2/g)**(1/4) * (rho_p/rho_t)**(1/4) * randomimpactors**(3/4)
    D = np.round(D/1000,2)             # Convert m to km and round to two decimals
    return D

# Random position and crater size generator
def RandomPositionAndArea(D1, D2, N=4000,
                          xlimlow=0.0, xlimhigh=500.0,
                          ylimlow=0.0, ylimhigh=500.0,
                          dpi=72, SFD=SFD, 
                          v_i=10000, g=9.81,
                          rho_p=3000, rho_t=5500, 
                          slope=1, q=3):
    """Random position and crater size generator.
    
       Parameters
       ----------
           D1 : float
               Lower limit for crater diameter in km
           D2 : float
               Upper limit for crater diameter in km
           N : integer
               Number of values generated
               Default is 3000 
           xlimlow : float
               Lower limit of x-positions
               Default is 0.0
           xlimhigh : float
               Upper limit of x-positions
               Default is 500.0
           ylimlow : float
               Lower limit of y-positions
               Default is 0.0
           ylimhigh : float
               Upper limit of y-positions
               Default is 500.0
           dpi : float
               Dots per inch of the plot
               Default is 72
           SFD : function
               Probability distribution function with the option to change the slope and y-intercept
               Format: slope * x^(-q)
               Default is the function SFD
           v_i : float
               Impact velocity in m/s
               Default is 10 km/s
           g : float
               Gravitational constant in m/s^2
               Default is 9.81 m/s^2
           rho_p : float
               Mass density of impactor in kg/m^3
               Default is 3000 kg/m^3
           rho_t : float
               Mass density of surface in kg/m^3
               Default is 5500 kg/m^3
           slope : float
               Slope of the probability function
               Default is 1
           q : float
               y-intercept of the probability function
               q >= 2 for a reasonable distribution
               Default is 3
                      
       Returns
       -------
           Uniform random x-positions, y-positions and crater sizes. 
       """

    # Rescale crater sizes for a 500x500 plot
    scale = dpi 
    relativecraters = scale * (CraterCreator(D1,D2,SFD=SFD,size=N,v_i=v_i,
                                                   g=g,rho_p=rho_p,rho_t=rho_t,
                                                   slope=slope,q=q))

    # Calculate random positions and craters sizes
    xposition = np.random.uniform(low=xlimlow, high=xlimhigh, size=N)
    yposition = np.random.uniform(low=ylimlow, high=ylimhigh, size=N)
    return xposition, yposition, relativecraters

# Loop 
def CraterSimulation(D1, D2, N=4000,
                     xlimlow=0.0, xlimhigh=500.0,
                     ylimlow=0.0, ylimhigh=500.0, 
                     dpi=72, SFD=SFD,
                     v_i=10000, g=9.81,
                     rho_p=3000, rho_t=5500,
                     slope=1, q=3):
    """Simulating craters on a given surface area and calculating the times of saturation. 
       Returns lists of values.
      
       Parameters
       ----------
           D1 : float
               Lower limit for crater diameter in km
           D2 : float
               Upper limit for crater diameter in km
           N : integer
               Number of values generated
               Default is 3000 
           xlimlow : float
               Lower limit of x-positions
               Default is 0.0
           xlimhigh : float
               Upper limit of x-positions
               Default is 500.0
           ylimlow : float
               Lower limit of y-positions
               Default is 0.0
           ylimhigh : float
               Upper limit of y-positions
               Default is 500.0
           dpi : float
               Dots per inch of the plot
               Default is 72
           SFD : function
               Probability distribution function with the option to change the slope and y-intercept
               Format: slope * x^(-q)
               Default is the function SFD
           size : float
               Number of impactors
               Default is 1
           v_i : float
               Impact velocity in m/s
               Default is 10 km/s
           g : float
               Gravitational constant in m/s^2
               Default is 9.81 m/s^2
           rho_p : float
               Mass density of impactor in kg/m^3
               Default is 3000 kg/m^3
           rho_t : float
               Mass density of surface in kg/m^3
               Default is 5500 kg/m^3
           slope : float
               Slope of the probability function
               Default is 1
           q : float
               y-intercept of the probability function
               q >= 2 for a reasonable distribution
               Default is 3
           
       Returns
       -------
           List of x and y values and craters sizes at 25%, 50%, 75% and 100% saturation : list
           List of times of saturation and the number of craters : list
           List of the used random x and y positions and crater sizes

    """
    # Plot
    fig = plt.figure(figsize=(6.0,6.0),dpi=dpi)
    
    # Saving the times of saturation
    nsat25, nsat50, nsat75, nsat100 = [],[],[],[]

    # Store positions at the times of saturation
    position25, position50, position75, position100 = [],[],[],[]
    
    # Store saturation
    saturationvalues = []
    
    # Extract position and crater size values
    xposition, yposition, cratersizes = RandomPositionAndArea(D1, D2, N=N,
                                                              xlimlow=xlimlow, xlimhigh=xlimhigh,
                                                              ylimlow=ylimlow, ylimhigh=ylimhigh,
                                                              dpi=dpi, SFD=SFD, 
                                                              v_i=v_i, g=g,
                                                              rho_p=rho_p, rho_t=rho_t, 
                                                              slope=slope, q=q)
    
    # Create an array (same length as the plot array below) 
    # representing the saturated surface with all values equal to 255
    satsurface = 255*np.ones(559872)

    # Loop
    saturation = 0
    time = 0
    indices = list(range(len(xposition)))
    for i,x,y,area in zip(indices,xposition,yposition,cratersizes):
        
        # Count time
        time += 1000
        
        # Plot
        plt.scatter(x, y, s=area, c='black', alpha=1, marker='o')
        plt.xlim(xlimlow,xlimhigh)
        plt.ylim(ylimlow,ylimhigh)
        fig.tight_layout(pad=0)
        plt.axis('off')

        # Plot figure to RGB array
        fig.canvas.draw()
        figarray = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)

        # Crater plot to array, 0 for white and 255 for black (reverse order than in general)
        craterfig = 255-figarray

        # Calculate percentage of saturation on surface
        saturation = np.sum(craterfig)/np.sum(satsurface)
        saturationvalues.append(saturation)

        # Save saturated surface image for 25%, 50%, 75% and 100%
        if saturation > 0.25:
            if len(nsat25)==0:
                position25.append(i)
                nsat25.append(time)
                nsat25.append(saturation)
                print(saturation)
                print(i)
                print("mark 25%")

        if saturation > 0.50:
            if len(nsat50)==0:
                position50.append(i)
                nsat50.append(time)
                nsat50.append(saturation)
                print(saturation)
                print(i)
                print("mark 50%")

        if saturation > 0.75:
            if len(nsat75)==0:
                position75.append(i)
                nsat75.append(time)
                nsat75.append(saturation)
                print(saturation)
                print(i)
                print("mark 75%")
            
        if saturation == 1.0:
            if len(nsat100)==0:
                position100.append(i)
                nsat100.append(time)
                nsat100.append(saturation)
                print(saturation)
                print(i)
                print("mark 100%")
                break
          
    return position25,position50,position75,position100,nsat25,nsat50,nsat75,nsat100,xposition,yposition,cratersizes,saturationvalues;

# Sizes of craters vs saturation changes
def SaturationChanges(cratersizes):
    """Calculating saturation changes for each individual crater sizes generated during the simulation.
    
       Parameters
       ----------
           cratersizes : list or float
               List or float of crater sizes in km generated in the simulation
           
       Returns
       -------
           A list of the changes in saturation for each individual crater sizes.
           
       Example
       -------
           How much does the saturation of a surface change when a 10 km crater is created?
           >> print(SaturationChanges([10]))
           >> [0.00011252572016460906]
           
    """

    satchanges = []
    fig = plt.figure(figsize=(6.0,6.0))  
    for craters in cratersizes:

        plt.cla()
        plt.scatter(250, 250, s=craters, c='black', alpha=1)
        plt.xlim(0,500)
        plt.ylim(0,500)
        fig.tight_layout(pad=0)
        plt.axis('off')

        # Plot figure to RGB array
        fig.canvas.draw()
        figarray = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)

        # Crater plot to array, 0 for white and 255 for black (reverse order than in general)
        craterfig = 255-figarray

        # Calculate percentage of saturation on surface
        satsurface = 255*np.ones(len(craterfig))
        saturation = np.sum(craterfig)/np.sum(satsurface)
        satchanges.append(saturation)

    plt.clf()
    return satchanges

# How to determine that a crater is recognizable on a surface
def CraterCount(saturationvalues,satchanges):
    """Counting the craters over time by measuring the change in saturation.
    
       Parameters
       ----------
           saturationvalues : list
               List of surface saturation percentages generated during simulation
           satchanges : list
               List of changes is saturation for each individual crater sizes
           
       Returns
       -------
           A list of craters counted over time during the simulation."""

    # Store crater count
    cratersovertime = []

    cratercount = 0
    indices = list(range(len(saturationvalues)))
    prev_sat = 0

    for i,saturation,changes in zip(indices,saturationvalues,satchanges):

        # Difference in saturation
        diff_sat = saturation - prev_sat

        # Difference between changes in saturation
        diff_change_sat = diff_sat / changes

        # If the changes are more than 25%
        if diff_change_sat >= 0.25:
            cratercount += 1

        prev_sat = saturation    
        cratersovertime.append(cratercount)
        
    return cratersovertime

# Save saturation 
def SaveSaturation(n,x,y,cratersizes,saturation,time,cratercount,dpi=72):
    """Creates plots of the surface at different saturation times from previously calculated data, 
        then saves it to a png file.
        
       Parameters
       ----------
           n : integer
               number of run
           x : list
               List of random x positions on a plot at different times of saturation, extracted from function
           y : list
               List of random y positions on a plot at different times of saturation, extracted from function
           cratersizes : list
               List of random crater sizes on a plot at different times of saturation, extracted from function
           saturation : float
               Percentage of saturation for the title of the plot, extracted from function
           time : float
               List of time, extracted from function
           cratercount : list
               Crater count, used for the time in the title of the plot, extracted from function
           dpi : float
               Dots per inch of the plot
               Default is 72
           
       Returns
       -------
           Saved png image of plot."""
    fig = plt.figure(figsize=(6.0,6.0),dpi=dpi)
    plt.style.use('default')
    ax = plt.gca()
    ax.set_facecolor('xkcd:slate grey')
    plt.scatter(x, y, s=cratersizes, c='#f0f0f0', alpha=0.75, linewidth=2, marker='o')
    plt.suptitle("Saturation = {:.3f}%".format(saturation*100), fontsize="19")
    plt.title("Duration = {:.2e} years, {} Visible Craters.".format(time,cratercount[int(time/1000 - 1)]), fontsize="17")
    plt.xlim(0,500)
    plt.ylim(0,500)
    fig.savefig('craters{}-{}.png'.format(n,np.round(np.int(saturation*100))))
    plt.clf();
    
# Making a movie
def ImpactMovie(n,xposition,yposition,cratersizes,saturationvalues,time,cratercount):
    """Makes a movie about the impact simulation.
        
       Parameters
       ----------
           n : integer
               number of run
           xposition : list
               List of random x positions on a plot at different times of saturation, extracted from function
           yposition : list
               List of random y positions on a plot at different times of saturation, extracted from function
           cratersizes : list
               List of random crater sizes on a plot at different times of saturation, extracted from function
           saturationvalues : list
               Percentage of saturation calculated during simulation, extracted from function
           time : float
               List of time, extracted from function
           cratercount : list
               Crater count, used for the time in the title of the plot, extracted from function
           dpi : float
               Dots per inch of the plot
               Default is 72
           
       Returns
       -------
           Saved png image of plot.
       """
    
    writer = ani.FFMpegWriter(fps=25)

    fig = plt.figure(figsize=(6.0,6.0),dpi=72)

    with writer.saving(fig, 'Craters{}.mp4'.format(n), 100):

        count = 1
        # loop 
        for x,y,craters,time,sat,cratercount in zip(xposition,yposition,cratersizes,time,saturationvalues,cratercount):
            
            ax = plt.gca()
            ax.set_facecolor('xkcd:slate grey')
            plt.scatter(x, y, s=craters, c='#f0f0f0', alpha=0.75, linewidth=2, marker='o')
            plt.suptitle("Duration = {:.2e} years".format(time))
            plt.title("Saturation = {:.3f}%, {} Visible Craters.".format(sat*100,cratercount))
            plt.xlim(0,500)
            plt.ylim(0,500)
            # save the current plot as a movie frame
            writer.grab_frame()
            count += 1
    
    