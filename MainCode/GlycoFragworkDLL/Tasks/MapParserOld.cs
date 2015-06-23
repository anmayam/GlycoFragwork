using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using GlycoFragworkDLL.Classes;
using GlycoFragworkDLL.Utils;

namespace GlycoFragworkDLL.Tasks
{
   public  class MapParserOld
    {
        MapRecord _GlycoMap;
        Utilities _utils;
     
        GlypID.Scoring.clsScoringParameters _ScoringParameters;
        GlypID.HornTransform.clsHornTransformParameters _TransformParameters;
        double _binSize;
        double _PPMDiff;
 
        

        public MapParserOld()
        {
            _GlycoMap = new GlycoFragworkDLL.Classes.MapRecord();
            _utils = new GlycoFragworkDLL.Utils.Utilities();
            _ScoringParameters = new GlypID.Scoring.clsScoringParameters();
            _TransformParameters = new GlypID.HornTransform.clsHornTransformParameters();
           
            _binSize = 0.05;
            _PPMDiff = 20; 

        }

        public MapParserOld(MapRecord mrecord, GlypID.Scoring.clsScoringParameters score_params, GlypID.HornTransform.clsHornTransformParameters transform_params)
        {
            _GlycoMap = new MapRecord();
            _GlycoMap = mrecord; 
        }
       
        /// <summary>
        /// Runs through all mrecords in the map and searches for putative glycopeptide match based on mass
        /// </summary>
        /// <param name="_glycoMap"></param>
        /// <param name="fastaFile"></param>
        /// <param name="glycanListFile"></param>
        public void SearchMapForNGlycopeptides(ref MapRecord _glycoMap, string fastaFile, string glycanListFile)
        {
            GlypID.Sequence.clsSequence[] sequences = new GlypID.Sequence.clsSequence[1];
            GlypID.Glycan.clsGlycan[] glycans = new GlypID.Glycan.clsGlycan[1];
            GlypID.Glycopeptide.clsGlycopeptide GlycoPeptide = new GlypID.Glycopeptide.clsGlycopeptide();


            GlycoPeptide.LoadGlycansFromList(glycanListFile, ref glycans);
            GlycoPeptide.LoadNGlycopeptidesFromFasta(fastaFile, ref sequences);

            _glycoMap._AllMLNRecords.ForEach(delegate(MultiAlignRecord m)
          {
              GlypID.ETDScoring.clsETDScoringScanResults[] etdScoringResults = new GlypID.ETDScoring.clsETDScoringScanResults[1];
              etdScoringResults[0] = new GlypID.ETDScoring.clsETDScoringScanResults();
              etdScoringResults[0].mdbl_mono_mw = m.Mass;
              GlypID.enmGlycanType type = GlypID.enmGlycanType.NA; 
              GlycoPeptide.SearchForGlycopeptides(ref etdScoringResults, ref glycans, ref sequences, type, false) ;
              m._PeptideSeq = new string[etdScoringResults.Length];
              m._ProteinName = new string[etdScoringResults.Length];
              m._GlycanComposition = new string[etdScoringResults.Length];
              for (int i = 0 ; i < etdScoringResults.Length ; i++)
              {
                  m._GlycanComposition[i] = etdScoringResults[i].mstr_glycan_composition;
                  m._PeptideSeq[i] = etdScoringResults[i].mstr_pep_seq;
                  m._ProteinName[i] = etdScoringResults[i].mstr_pro_seq_name; 
              }
          });
        }

        /// <summary>
        /// Function to calculate representative fragmentation spectra for each multi record. 
        /// </summary>
        /// <param name="_glycoMap"></param>
        public void CalculateRepresentativeFragmentationSpectraThroughSum(ref Classes.MapRecord _glycoMap)
        {
            for (int i = 0; i < _glycoMap._AllMLNRecords.Count; i++)
            {
                MultiAlignRecord m = _glycoMap._AllMLNRecords[i];
                int num_records = m._AssociatedUMCRecords.Count;
                int num_cid_peaks = (int)((_TransformParameters.MaxMZ - _TransformParameters.MinMZ) / _binSize) + 1;
                int num_hcd_peaks = (int) ((_ScoringParameters.MaxHCDMz - _ScoringParameters.MinHCDMz) / _binSize) + 1;
                int num_etd_peaks = (int)((_TransformParameters.MaxMZ - _TransformParameters.MinMZ) / _binSize) + 1;
                GlypID.Peaks.clsPeak[] summedCIDPeaks = new GlypID.Peaks.clsPeak[num_cid_peaks];
                GlypID.Peaks.clsPeak[] summedHCDPeaks = new GlypID.Peaks.clsPeak[num_hcd_peaks];
                GlypID.Peaks.clsPeak[] summedETDPeaks = new GlypID.Peaks.clsPeak[num_etd_peaks];

                for (int j = 0; j < num_cid_peaks; j++)
                {
                    summedCIDPeaks[j] = new GlypID.Peaks.clsPeak();
                    summedCIDPeaks[j].mdbl_mz = _TransformParameters.MinMZ + j * _binSize;
                    summedCIDPeaks[j].mdbl_intensity = 0;
                }
                for (int j = 0; j < num_etd_peaks; j++)
                {
                    summedETDPeaks[j] = new GlypID.Peaks.clsPeak();
                    summedETDPeaks[j].mdbl_mz = _TransformParameters.MinMZ + j * _binSize;
                    summedETDPeaks[j].mdbl_intensity = 0; 
                }
                for (int j = 0; j < num_hcd_peaks; j++)
                {
                    summedHCDPeaks[j] = new GlypID.Peaks.clsPeak();
                    summedHCDPeaks[j].mdbl_mz = _TransformParameters.MinMZ + j * _binSize;
                    summedHCDPeaks[j].mdbl_intensity = 0; 
                }


                for (int j = 0; j < num_records; j++)
                {
                    if (m._AssociatedUMCRecords[j].CIDPeaks[0] != null)
                    {
                        GlypID.Peaks.clsPeak[] thisCIDPeaks = new GlypID.Peaks.clsPeak[m._AssociatedUMCRecords[j].CIDPeaks.Length];
                        Array.Copy(m._AssociatedUMCRecords[j].CIDPeaks, thisCIDPeaks, m._AssociatedUMCRecords[j].CIDPeaks.Length);
                        BinNormalize(ref thisCIDPeaks, (float)_TransformParameters.MinMZ, (float)_TransformParameters.MaxMZ, _binSize);
                        AddToPeaks(ref summedCIDPeaks, ref thisCIDPeaks);
                    }
                    if (m._AssociatedUMCRecords[j].HCDPeaks[0] !=null)
                    {
                        GlypID.Peaks.clsPeak[] thisHCDPeaks = new GlypID.Peaks.clsPeak[m._AssociatedUMCRecords[j].HCDPeaks.Length];
                        Array.Copy(m._AssociatedUMCRecords[j].HCDPeaks, thisHCDPeaks, m._AssociatedUMCRecords[j].HCDPeaks.Length);
                        BinNormalize(ref thisHCDPeaks, (float)_ScoringParameters.MinHCDMz, (float)_ScoringParameters.MaxHCDMz, _binSize);
                        AddToPeaks(ref summedHCDPeaks, ref thisHCDPeaks);
                    }
                    if (m._AssociatedUMCRecords[j].ETDPeaks[0] !=null)
                    {
                        GlypID.Peaks.clsPeak[] thisETDPeaks = new GlypID.Peaks.clsPeak[m._AssociatedUMCRecords[j].ETDPeaks.Length];
                        Array.Copy(m._AssociatedUMCRecords[j].ETDPeaks, thisETDPeaks, m._AssociatedUMCRecords[j].ETDPeaks.Length);
                        BinNormalize(ref thisETDPeaks, (float)_TransformParameters.MinMZ, (float)_TransformParameters.MaxMZ, _binSize);
                        AddToPeaks(ref summedETDPeaks, ref thisETDPeaks);
                    }
                }
                m._RepresentativeCIDPeaks = summedCIDPeaks;
                m._RepresentativeETDPeaks = summedETDPeaks;
                m._RepresentativeHCDPeaks = summedHCDPeaks; 
                
            }
        }

        public void CalculateRepresentativeParentTransform(ref Classes.MapRecord _glycoMap)
        {
            // This might be problematic
            for (int i = 0; i < _glycoMap._AllMLNRecords.Count; i++)
            {
                MultiAlignRecord m = _glycoMap._AllMLNRecords[i];
                int num_records = m._AssociatedUMCRecords.Count;
                double min_ppm = _PPMDiff;
                int min_ppm_index = -1;
                for (int j = 0; j < num_records; j++)
                {
                    double ppm = _utils.CalculateDelMassPPM(m.Mass, m._AssociatedUMCRecords[j].TransformResult.mdbl_mono_mw);
                    if (ppm < min_ppm)
                    {
                        min_ppm_index = j;
                        min_ppm = ppm;
                    }
                }
                m._RepresentativeTransformResult = m._AssociatedUMCRecords[min_ppm_index].TransformResult; 
            }
        }

        public void AddToPeaks(ref GlypID.Peaks.clsPeak[] peak, ref GlypID.Peaks.clsPeak[] peaksToAdd)
        {
            if (peak.Length != peaksToAdd.Length)
                throw new Exception("Peaks to sum must be equal length");

            for (int i = 0; i < peaksToAdd.Length; i++)
            {
                peak[i].mdbl_intensity = peak[i].mdbl_intensity + peaksToAdd[i].mdbl_intensity; 
            }
            
        }

        /// <summary>
        /// Function to normalize peaks into bins of set width 
        /// </summary>
        /// <param name="?"></param>
        public void BinNormalize(ref GlypID.Peaks.clsPeak[] peaks, float min_mz, float max_mz, double bin_size)
        {
            int num_bins = (int)((max_mz - min_mz) / bin_size)+1; 
            GlypID.Peaks.clsPeak[] baseSpectrumPeaks = new GlypID.Peaks.clsPeak[num_bins];

            for (int i = 0; i < baseSpectrumPeaks.Length; i++)
            {
                baseSpectrumPeaks[i] = new GlypID.Peaks.clsPeak(); 
                baseSpectrumPeaks[i].mdbl_mz = min_mz + i * bin_size;
                baseSpectrumPeaks[i].mdbl_intensity = 0; 
            }

            for (int i = 0; i < peaks.Length; i++)
            {
                int bin = (int) ((peaks[i].mdbl_mz - min_mz) / bin_size);
                baseSpectrumPeaks[i].mdbl_intensity += peaks[i].mdbl_intensity; 
            }

            Array.Clear(peaks, 0, peaks.Length);
            peaks = baseSpectrumPeaks; 
        }
    }
}
