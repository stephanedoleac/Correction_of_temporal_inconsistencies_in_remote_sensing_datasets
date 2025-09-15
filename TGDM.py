import numpy as np
import math

class TGDM() :
    """
    Class implementing the Temporal Gap Detection Method (TGDM) 
    described in van Oostende et al. (2022).
    
    It creates an object containing:
     - the initial time series,
     - the summed observations time series,
     - the yearly mask, 
     - the corrected time series.
    
    Remarks:
     - Only works for an individual time series with daily frequency. 
     - The time series must begin on January 1st and end on December 31st.
     - Leap days (February 29th) must be removed before use.
    """
    
    def __init__(self, time_serie, window, n_years) :
        """
        Initialize a TGDM object.

        Parameters
        ----------
        time_serie : np.array
            Original time series to be corrected.
        window : int
            Window size (in days) for computing summed observations.
        n_years : int
            Number of years in the time series.
        """
        
        self.time_serie = time_serie
        self.window = window
        self.nbr_years = n_years        
        self.len = time_serie.shape[0]
        
        # Arrays to store intermediate and final results
        self.summed_obs = np.empty(self.len)
        self.mask = np.empty(365)
        self.corrected_time_serie = np.empty(self.len)
        
        # Build TGDM outputs immediately upon initialization
        self.build_TGDM()
        
    
    def build_TGDM(self) :
        """
        Run the TGDM algorithm on the input time series.
        Produces:
         - summed observations time series
         - minimum observations per day (mask)
         - corrected time series
        """
        
        time_serie = self.time_serie
        window = self.window
        n_years = self.nbr_years
    
        # Step 1: compute summed observations per window
        summed_obs = self.summed_observations_time_serie(time_serie, window)
        
        # Step 2: compute minimum observations per day of year
        min_obs = self.min_obs_per_day(summed_obs, n_years)
    
        # Step 3: mask out days without sufficient observations
        time_serie_TGDM = self.remove_days_without_obs(time_serie, min_obs, n_years)
    
        # Save results
        self.summed_obs = summed_obs
        self.mask = min_obs
        self.corrected_time_serie = time_serie_TGDM

    
    def summed_observations_time_serie(self, time_serie, window) :
        """
        Compute the summed observations time series (see Figure 2B in van Oostende et al. 2022). 
        This gives, for each day, the number of valid observations within a moving window.
        The first and last 'half_window' days are set to NaN.

        Parameters
        ----------
        time_serie : np.array
            Original daily time series.
        window : int
            Window size (in days).

        Returns
        -------
        np.array
            Time series of summed observations.
        """
    
        n = len(time_serie)
        half_window = math.floor(window/2)
        
        summed_obs = np.empty([n-window+1])
    
        # Slide the window across the time series
        for k in range(half_window, len(time_serie)-window+1) :
            values_in_window = time_serie[k-half_window:k+half_window+1]
            # Count number of non-NaN values in the window
            summed_obs[k] = window - np.isnan(values_in_window).sum()
        
        # Add NaN padding at the beginning and end
        return np.concatenate((np.array([np.nan]*half_window), summed_obs, np.array([np.nan]*half_window)))


    def min_obs_per_day(self, summed_obs, n_years) :
        """
        Compute the minimum number of observations ever recorded 
        for each day of the year (see Figure 2C in van Oostende et al. 2022).

        Parameters
        ----------
        summed_obs : np.array
            Output from summed_observations_time_serie.
        n_years : int
            Number of years in the time series.

        Returns
        -------
        np.array
            Array of length 365 containing the minimum number of 
            observations per day-of-year.
        
        Remarks
        -------
        Only valid if the time series starts on January 1st and 
        ends on December 31st (without leap days).
        """
        
        min_number_of_obs_per_day = np.empty([365])
        
        # For each day of the year, check across all years
        for day in range(365) :
            min_number_of_obs_per_day[day] = np.nanmin(summed_obs[day:-1:365])
            
        return min_number_of_obs_per_day
        

    def remove_days_without_obs(self, time_serie, min_obs_per_day, n_years) :
        """
        Mask out the days where the minimum number of observations is zero.
        Returns a corrected time series with NaNs on those days.

        Parameters
        ----------
        time_serie : np.array
            Original time series.
        min_obs_per_day : np.array
            Output from min_obs_per_day.
        n_years : int
            Number of years in the time series.

        Returns
        -------
        np.array
            Corrected time series with NaNs where coverage is insufficient.
        
        Remarks
        -------
        Only valid if the time series starts on January 1st and 
        ends on December 31st (without leap days).
        """
    
        # Build mask (NaN if no observations, 1 otherwise)
        yearly_mask = np.where(min_obs_per_day == 0, np.nan, 1)            
            
        # Apply the mask year by year
        time_serie_TGDM = np.empty([len(time_serie)])*np.nan
        
        for yr in range(n_years) :
            time_serie_TGDM[yr*365:(yr+1)*365] = time_serie[yr*365:(yr+1)*365]*yearly_mask
        
        return time_serie_TGDM
