
import numpy as np
import matplotlib.pyplot as plt


#---------------------------#
# FUNCTION: readdata        #
#---------------------------#

def readdata(filename):
    """
    Reads the data of the CSV file

    Parameters
    ----------
    filename : str
        The name of the file that will be read.

    Returns
    -------
    numpy.ndarray
        The file after it has been turned into a numpy array format.
    """
    # Purpose: Read in Harvard Forest data
    return np.genfromtxt(filename, delimiter= "," , skip_header= 1)
    

#---------------------------#
# FUNCTION: summarizedata   #
#---------------------------#

def summarizedata(hf):
    """
    Returns data relating to the file and plots CO2 data relating to time.

    Parameters
    ----------
    hf : numpy.ndarray
        Array that contains the Harvard Forest data.

    Returns
    -------
    ndata : int
        Number of data points.
    hfmean : float
        Average value of CO2 flux.
    hf25 : float
        25th quartile of of CO2 flux.
    hf75 : float
        75th quartile of of CO2 flux.
    """
    # Purpose:
    #	Create a plot of the CO2 data as a time series
    #	Calculate several summary statistics about the dataset.
    # Read in data from the Harvard Forest
    # Create a time series plot of the data
    
    #defining the variables for the plot
    x = hf[:,0] + (hf[:,2] / 365)
    y = hf[:,3]
    #creating the labels of the plot
    plt.title("CO2 Flux vs. Time") 
    plt.xlabel("Time (years)") 
    plt.ylabel("CO2 Flux")
    #making the plot
    plt.scatter(x, y)
    
    # Display the plot
    plt.show()

    # Find the number of data points
    ndata = len(x)
    
    # Find the mean of the data
    total_y = sum(y)
    hfmean = round((total_y / ndata) , 3)
    
    # Find the 25th percentile of the data points
    hf25 = round(np.quantile(y , 0.25), 3)
    
    # Find the 75th percentile of the data points
    hf75 = round(np.quantile(y , 0.75),3)
    
    # Return the number of data points, the mean, the 25th percentile, and the 75th percentile
    return ndata,hfmean,hf25,hf75

#---------------------------#
# FUNCTION: missingdata     #
#---------------------------#

def missingdata(hf):
    """
    Returns the percentage of missing data from each year and plots it.

    Parameters
    ----------
    hf : numpy.ndarray
        Array that contains the Harvard Forest data.

    Returns
    -------
    missing_data : list
        List of the perentage of the missing data for each year.
    """
    # Purpose: 
    # Find the percentage of missing data points in each year
    # Read in the data from the Harvard Forest
    # Find the percentage of data in each year that are missing
    
    #finding min and max year in the data
    min_y = int(np.min(hf[:,0]))
    max_y = int(np.max(hf[:,0]))
    
    #making arrays for x and missing data that will be used for plots
    missing_data = []
    x = []
    #iterating through the years and counting how many days there are per year
    for i in range(min_y, max_y + 1):
        count_missing = np.size((np.array(np.where(hf[:,0] == i))))
        #checks if year is a leap year and then calculates % missing data
        if i % 4 == 0:
            miss_pcnt = round(((366 - count_missing) / 366 ) * 100 , 3)
        else:
            miss_pcnt = round(((365 - count_missing) / 365 ) * 100 , 3)
        #adds found info to the arrays for x and missing data
        missing_data.append(miss_pcnt)
        x.append(i)
    # Create a bar plot showing the percentage of data that are
    # missing in each year
    plt.bar(x,missing_data)
    #adding labels
    plt.title("Percent Data Missing Each Year") 
    plt.xlabel("Year") 
    plt.ylabel("Percentage of data missing") 
    
    # Display the plot
    plt.show()
    
    # Return the result
    return missing_data

#---------------------------#
# FUNCTION: seasonalcycle   #
#---------------------------#

def seasonalcycle(hf):
    """
    Finds the mean CO2 flux for each month and then plots those values from
    1995 to 2000 inclusive.

    Parameters
    ----------
    hf : numpy.ndarray
        Array that contains the Harvard Forest data.

    Returns
    -------
    mean_month : numpy.ndarray
        An array that contains the years and months and then the mean CO2 flux
        for that year/month.
    """
    # Purpose: 
    # 	Find the average flux by month for each year
    #	Plot the average monthly flux for years 1995 through (and including) 2000.
    # Read in the data from the Harvard Forest
    
    # Determine the range of years
    min_y = int(np.min(hf[:,0]))
    max_y = int(np.max(hf[:,0]))
    dif = max_y - min_y + 1
    
    # Find the average CO2 flux in each month
    
    #defining empty array to be filled
    mean_month = np.zeros((dif,12))
    #iterating through years and then months
    for i in range(min_y, max_y + 1):
        for month in range(1,13):
            #getting the index and mean for each month
            index = (hf[:,0] == i) & (hf[:,1] == month)
            co2_flux = hf[index,3]
            month_avg = np.mean(co2_flux)
            #adding that found value to the empty array
            mean_month[i - min_y , month - 1] = month_avg
    # Plot the average fluxes by month in several different years
    x = np.arange(1,13)
    #getting the y values for each year and plotting them
    for year in range(1995,2001):
        index = year - min_y
        plt.plot(x, mean_month[index], label=str(year) )
    #adding labels
    plt.title("Average CO2 Flux per Year") 
    plt.xlabel("Month") 
    plt.ylabel("Mean CO2 Flux per month") 
    plt.legend()
    
    # Display the plot
    plt.show()
    
    # Return the result
    return mean_month

#---------------------------#
# FUNCTION: HFregression    #
#---------------------------#

def HFregression(hf):
    """
    Calculates the model for what the CO2 flux value should be and then plots
    that against the actual CO2 vlaue.

    Parameters
    ----------
    hf : numpy.ndarray
        Array that contains the Harvard Forest data.

    Returns
    -------
    beta : numpy.ndarray
        Regression coefficients that are used in the regression model.
    model_est : numpy.ndarray
        Model estimate of for the CO2 flux values.
    """
    # Purpose: 
        #    Create a regression model for CO2 fluxes
        #    Visualize the outputs of the model
    # Read in the data from the Harvard Forest
    # Create the X matrix for the regression
    
    
    # Estimate the regression coefficients
    
    #getting co2 flux and environmnetal factors
    z = hf[:, 3]
    env_fact = hf[:, 4:]
    #creating the X column of ones and then adding to env_fact
    X_ones = np.ones((len(hf), 1))
    X = np.concatenate((X_ones, env_fact) , axis = 1)
    #calculating beta
    beta = np.linalg.inv( np.transpose(X) .dot(X) ) .dot(np.transpose(X)) .dot(z)
    
    # Create the model estimate
    model_est = X.dot(beta)
    
    # Calculate the correlation coefficient
    corr_coef = round(np.corrcoef(z, model_est)[0,1],3)
    
    # Plot the model estimate and add the correlation coefficient to the plot
    
    #defining the x axis if the plot
    years = hf[:,0] + (hf[:,2] / 365)
    
    #make the figuresize to fit both subplots
    plt.figure(figsize=(12, 11))
    #create the first subplot
    plt.subplot(211)
    plt.title("Time Series of Model Fit") 
    plt.xlabel("Time (years)") 
    plt.ylabel("CO2 Flux")

    #making the plot
    plt.plot(years, z, label = "Data of CO2 Flux" )
    plt.plot(years, model_est, label = "Estimated CO2 FLux")
    plt.legend()
    plt.text(0.2, 0.85, f'r = {corr_coef}', transform=plt.gca().transAxes)
    
    # Create a plot of the model components
    Xbeta = X * beta
    #create the legend names
    label_names = ["Intercept", "Net Radiation", "Air Temperature", "Water Vapor", "Wind Speed"]
    
    #create the second plot with all the components
    plt.subplot(212)
    for i in range(Xbeta.shape[1]):
        plt.plot(years, Xbeta[:, i], label= label_names[i])
    plt.xlabel('Time (years)')
    plt.ylabel('CO2 FLux')
    plt.title('Contribution of Different Model Components')
    plt.legend()

    # Display the plot
    plt.show()
    
    # Return the regression coefficients 
    return beta , model_est


#----------------------------#
# FUNCTION: averagecarbon    #
#----------------------------#

def averagecarbon(hf, modelest):
    """
    Calculates the average CO2 flux value from the model numbers it then plots
    that number for every year and shows if the forest was a net sink or not.

    Parameters
    ----------
    hf : numpy.ndarray
        Array that contains the Harvard Forest data.
    model_est : numpy.ndarray
        Model estimate of for the CO2 flux values.

    Returns
    -------
    avg_c : numpy.ndarray
        Average CO2 flux value according to the model.
    """
    # PURPOSE: calculate the average carbon flux from the 
    # model for each year of the simulations. Create a time series plot
    # showing the carbon flux by year.
    # Read in the data from the Harvard Forest
    
    # Calculate the average, modeled CO2 flux for each year
    
    #creating an array of all the years
    min_y = int(np.min(hf[:,0]))
    max_y = int(np.max(hf[:,0]))
    all_years = np.arange(min_y, max_y + 1)
    dif = max_y - min_y + 1
    #creating empty array
    avg_c = np.zeros(dif)
    #iterating through the years and finding the average carbon value
    for i in all_years:
        index = np.where(hf[:, 0] == i)[0]
        c_total = np.sum(modelest[index])
        days = len(index)
        avg = round((c_total/days) , 3)
        avg_c[i-min_y] = avg
        
    # Create a plot of points showing the average modeled CO2 flux for each year
    plt.scatter(all_years, avg_c)
    plt.title("Average Annual CO2 Flux at Harvard Forest")
    plt.xlabel("Year")
    plt.ylabel("Average CO2 Flux")
    plt.axhline(0, color='k')
    
    # Display the plot
    plt.show()
    print(type(avg_c))
    # Return the result
    return avg_c

"""
It seems as if the forest was a carbon sink until around 2000 and then there
might have been deforestation or something in 200 that caused there to be 
carbon emitted form the forest for a few years. Then recently the forest has
started being a large carbon sink. It has taken in the most carbon ever in
the more recent years.
"""

#-----------------------------------------#
# Execute the functions defined above     #
#-----------------------------------------#

if __name__ == "__main__": 
    filename               = 'harvard_forest-1.csv'
    hf                     = readdata(filename)
    #summarizedata(hf)
    #missingdata(hf)
    #seasonalcycle(hf)
    #HFregression(hf)
    
#    ndata,hfmean,hf25,hf75 = summarizedata(hf)
#    missing_data           = missingdata(hf)
#    month_means            = seasonalcycle(hf)
#    betas, modelest        = HFregression(hf)
#    means                  = averagecarbon(hf, modelest) 
