using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GlycoFragworkDLL.Utils
{
    /// <summary>
    /// Class for performing operations on a spectrum/spectra
    /// Data is stored here in the form of lists.
    /// </summary>
    public class SpectralUtilities
    {

        List <List<double>> _AllBinnedSpectra ;
        List <List<double>> _Clusters;
        List <string> _SpectraNames ; 
        List <string> _ClusterNames ;
        List<List<GlypID.Peaks.clsPeak>> _AllOriginalSpectra; 

        float _minMz = 100;
        float _maxMz = 2000;
        double _binsize = 1.0;
        int _rounds = 20;
        float _Tmin = Convert.ToSingle(0.6);
        float _Tau = 1;
        float _delta; 


        public SpectralUtilities()
        {
            _AllBinnedSpectra = new List<List<double>>() ;
            _Clusters = new List<List<double>>(); 
            _SpectraNames = new List<string>();
            _ClusterNames = new List<string>() ;
            _AllOriginalSpectra = new List<List<GlypID.Peaks.clsPeak>>(); 
        }

        public void Clear()
        {
            _AllBinnedSpectra.Clear();
            _Clusters.Clear(); 
            _SpectraNames.Clear();
            _ClusterNames.Clear();
            _AllOriginalSpectra.Clear(); 
        }

        /// <summary>
        /// Add a bunch of peaks to a list of spectra, bin them and add the binned spectra to another list
        /// </summary>
        /// <param name="peaks"></param>
        /// <param name="spectraName"></param>
        public void AddPeaksToList(ref GlypID.Peaks.clsPeak[] peaks, string spectraName)
        {
            // add original spectrum
            List<GlypID.Peaks.clsPeak> thisOrigSpectra = new List<GlypID.Peaks.clsPeak>();
            for (int i = 0; i < peaks.Length; i++)
            {
                thisOrigSpectra.Add(peaks[i]);
            }
            _AllOriginalSpectra.Add(thisOrigSpectra); 


            // Bin spectrum
            List<double> thisBinnedSpectra = new List<double>() ; 
            BinNormalize(ref peaks, ref thisBinnedSpectra, _minMz, _maxMz, _binsize) ;          

            // Add binned spectrum to another list and to Cluster object
            _AllBinnedSpectra.Add(thisBinnedSpectra);
            _Clusters.Add(thisBinnedSpectra); 

            // Store names to keep track
            _ClusterNames.Add(spectraName);
            _SpectraNames.Add(spectraName); 
        }

       
        /// <summary>
        ///  Function to perform hieratchial clustering spectra based on cosine similarity.  Each spectrum must have been binned and added to cluster object
        /// </summary>
        /// <returns></returns>
        public int ClusterSpectraInList()
        {
            int num_clusters = 0;
            if (_Clusters.Count>0)            
            {
                _delta = (1 - _Tmin) / _rounds; 
                _Tau = (float)0.9+_delta;

                for (int i = 0; i < _rounds; i++)
                {
                    _Tau = _Tau - _delta;
                    for (int j = 1; j < _Clusters.Count; j++)
                    {
                        for (int k = 0; k < j; k++)
                        {
                            if (j < _Clusters.Count)
                            {
                                List<double> a = _Clusters[j];
                                List<double> b = _Clusters[k];
                                double dot = 0;
                                double val1 = 0;
                                double val2 = 0;
                                for (int n = 0; n < a.Count(); n++)
                                {
                                    dot += a[n] * b[n];
                                    val1 += Math.Pow(a[n], 2);
                                    val2 += Math.Pow(b[n], 2);
                                }
                                double CosineSimilarity = (dot / (Math.Sqrt(val1) * Math.Sqrt(val2)));

                                if (CosineSimilarity >= _Tau)
                                {
                                    for (int l = 0; l < _Clusters[k].Count(); l++)
                                    {
                                        _Clusters[k][l] = ((_Clusters[j][l] + _Clusters[k][l]) / 2);

                                    }
                                    _SpectraNames[k] = "(" + _SpectraNames[k] + " and " + _SpectraNames[j] + "-{Cosim:" + Convert.ToString(CosineSimilarity) + "})";
                                    _ClusterNames[k] = _ClusterNames[k] + "-" + _ClusterNames[j];
                                    _Clusters.RemoveAt(j);
                                    _SpectraNames.RemoveAt(j);
                                    _ClusterNames.RemoveAt(j);
                                }
                            }
                        }
                    }

                }

                num_clusters = _ClusterNames.Count; 
            }

            return num_clusters; 
        }

        public void AssignClusters(ref List<string> clusternames)
        {
            _ClusterNames.Clear();
            foreach (string s in clusternames)
            {
                _ClusterNames.Add(s); 
            }

        }

        // Get cluster names
        public void GetClusterNames(ref List<string> clusternames)
        {
            clusternames.Clear();
            foreach (string s in _ClusterNames)
            {
                clusternames.Add(s);
            }
        }

        // Given a cluster index, return the list of spectra in the form of original indices
        public void GetOriginalIDFromCluster(int clusterIndex, ref List<int> index)
        {
            string clustername = _ClusterNames[clusterIndex];
            string[] spectra = clustername.Split('-');

            if (spectra.Length > 0)
            {
                foreach (string s in spectra)
                {
                    string[] parts = s.Split('_');
                    int spectra_index = Convert.ToInt32(parts[2]);
                    index.Add(spectra_index);
                }
            }
            else
            {             
                // Should never reach here 
                index.Add(0); 
            }
        }
        
        
      

        /// <summary>
        /// Given a cluster index, return the representative of that cluster. Currently, the spectrum with maximum SNR in a given range is returned
        /// </summary>
        /// <param name="clusterIndex"></param>
        /// <param name="repPeaks"></param>
        /// <param name="minmz"></param>
        /// <param name="maxmz"></param>
        /// <param name="return_original_id">true or false. If true, return the original index(among all spectra), else returns the index within the cluster</param>
        /// <returns></returns>
        public int GetRepresentativePeaksFromCluster(int clusterIndex, ref GlypID.Peaks.clsPeak[] repPeaks, double minmz, double maxmz, bool return_original_id)
        {
            double max_snr = 0;
            int spectra_with_max_snr = -1;
            int orig_index = -1; 
            string clustername = _ClusterNames[clusterIndex];
            string[] spectra = new string[1];

            if (clustername.Contains('-'))
                spectra = clustername.Split('-');
            else
                spectra[0] = clustername; 

            if (spectra.Length > 0)
            {
                foreach (string s in spectra)
                {
                    List<GlypID.Peaks.clsPeak> this_peaks = new List<GlypID.Peaks.clsPeak>();
                    string[] parts = s.Split('_');
                    int spectra_index = Convert.ToInt32(parts[2]);
                    this_peaks = _AllOriginalSpectra[spectra_index];
                    double noise_level = CalculateAverage(ref this_peaks, 0, minmz, maxmz);
                    double signal_level = CalculateAverage(ref this_peaks, noise_level, minmz, maxmz);
                    double snr = signal_level / noise_level;
                    if (snr > max_snr)
                    {
                        max_snr = snr;
                        spectra_with_max_snr = spectra_index;
                        orig_index = Convert.ToInt32(parts[1]);
                    }
                }
                if (spectra_with_max_snr >= 0)
                {
                    repPeaks = new GlypID.Peaks.clsPeak[_AllOriginalSpectra[spectra_with_max_snr].Count];
                    for (int i = 0; i < _AllOriginalSpectra[spectra_with_max_snr].Count; i++)
                    {
                        repPeaks[i] = _AllOriginalSpectra[spectra_with_max_snr][i];
                    }
                }
            }
         
            if (return_original_id)
                return orig_index;
            else
                return spectra_with_max_snr;                 
            
        }

        /// <summary>
        /// Sinple function to calculate average of intensities within an mz range and above a threshold
        /// </summary>
        /// <param name="peaks"> List of cls peaks</param>
        /// <param name="threshold">threshold</param>
        /// <param name="min_mz">min mz</param>
        /// <param name="max_mz">max mz</param>
        /// <returns> the average</returns>
        public double CalculateAverage(ref List<GlypID.Peaks.clsPeak> peaks, double threshold, double min_mz, double max_mz)
        {
            double sum = 0;
            int count = 0; 
            double avg = 0 ; 
            for (int i = 0; i < peaks.Count; i++)
            {
                if (peaks[i].mdbl_mz >= min_mz && peaks[i].mdbl_mz <= max_mz)
                {
                    if (peaks[i].mdbl_intensity > threshold)
                    {
                        sum += peaks[i].mdbl_intensity;
                        count++;
                    }
                }
            }
            if (count > 0)
                avg = (sum/count);
            return avg; 
        }



        /// <summary>
        /// Function to normalize peaks into bins of set width within mz range
        /// </summary>
        /// <param name="?"></param>
        /// Return values in binnedIntensities
        public void BinNormalize(ref GlypID.Peaks.clsPeak[] peaks, ref List<double> binnedIntensities, float min_mz, float max_mz, double bin_size)
        {
            int num_bins = (int)((max_mz - min_mz) / bin_size)+1;
            binnedIntensities = new List<double>();
            for (int i =0; i < num_bins;i++)
            {
                binnedIntensities.Add(0) ; 
            }
            
            for (int i = 0; i < peaks.Length; i++)
            {
                int bin = (int) Math.Round ((peaks[i].mdbl_mz - min_mz) / bin_size);
                binnedIntensities[bin]+= peaks[i].mdbl_intensity; 
            }
           
        }

        public string EncodePeaks(ref GlypID.Peaks.clsPeak[] peaks)
        {
            string spectrum = "";

            if (peaks != null)
            {
                float[] mzs = new float[peaks.Length];
                float[] intensities = new float[peaks.Length];

                for (int i = 0; i < peaks.Length; i++)
                {
                    mzs[i] = (float)peaks[i].mdbl_mz;
                    intensities[i] = (float)peaks[i].mdbl_intensity;
                }

                spectrum = GlypID.Utils.EncodeData(ref mzs, ref intensities);
            }
            else
            {
                spectrum = "";
            }

            return spectrum; 
        }

    }
}
