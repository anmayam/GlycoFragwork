using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GlycoFragworkDLL.Classes
{
    public class FragEvents
    {
        GlypID.enmProfileType _cidProfileType;
        GlypID.enmProfileType _hcdProfileType;
        GlypID.enmProfileType _etdProfileType;

        int _id; 

        int _hcd_scan;
        int _cid_scan;
        int _etd_scan;
        int _parent_scan;

        double _parent_scan_time; 

        double _parent_mz;
       /* double _hcd_parent_mz;
        double _etd_parent_mz;
        double _cid_parent_mz; */
        GlypID.HornTransform.clsHornTransformResults _TransformResult;

        float[] _cid_mz_values;
        float[] _cid_intensity_values;
        GlypID.Peaks.clsPeak[] _CIDPeaks;

        float[] _hcd_mz_values;
        float[] _hcd_intensity_values;
        GlypID.Peaks.clsPeak[] _HCDPeaks;

        float[] _etd_mz_values;
        float[] _etd_intensity_values;
        GlypID.Peaks.clsPeak[] _ETDPeaks;

        GlypID.enmGlycanType _glycanType;
        double _hcdScore;
        double _etdScore;
        double _cidScore;

        GlycopeptideRecord _gp_record;
        float _fdr;
        double _cid_sequencing_score;

        bool _FalseHit; 

        /*string _ProteinName;
        string _PeptideSeq;
        string _NGlycoSite;
        double _PeptideMass;
        string _GlycanComposition;
        double _GlycanMass;
        double _Y1;
        double _PPM;*/


        public FragEvents()
        {
            _id = -1; 

            _TransformResult = new GlypID.HornTransform.clsHornTransformResults();

            _cid_intensity_values = new float[1];
            _cid_mz_values = new float[1];
            _CIDPeaks = new GlypID.Peaks.clsPeak[1];

            _hcd_mz_values = new float[1];
            _hcd_intensity_values = new float[1];
            _HCDPeaks = new GlypID.Peaks.clsPeak[1];

            _etd_intensity_values = new float[1];
            _etd_mz_values = new float[1];
            _ETDPeaks = new GlypID.Peaks.clsPeak[1];

            _glycanType = GlypID.enmGlycanType.NA;
            _hcdScore = 1;
            _etdScore = 0;
            _cidScore = 0;

            _hcd_scan = 0;
            _etd_scan = 0;
            _cid_scan = 0;

            _fdr = 1.0f;
            _cid_sequencing_score = 0; 



            /*_ProteinName = "";
            _PeptideSeq = "";
            _NGlycoSite = "";

            _PeptideMass = 0;
            _GlycanComposition = "";
            _GlycanMass = 0;
            _PPM = 0;
            _Y1 = 0;*/

            _gp_record = new GlycopeptideRecord();
            _parent_mz = 0;
            _parent_scan =  0;
            _parent_scan_time = 0; 

            _FalseHit = false; 
            /*_hcd_parent_mz = 0;
            _etd_parent_mz = 0;
            _cid_parent_mz = 0;*/ 
        }

        public FragEvents(FragEvents e)
        {
            _id = e._id; 
            _TransformResult = e._TransformResult;
            _parent_mz = e._parent_mz;
            _parent_scan = e._parent_scan;
            _parent_scan_time = e._parent_scan_time; 
            _gp_record = e._gp_record;
            _glycanType = e._glycanType; 

            _cid_scan = e._cid_scan; 
            _cid_intensity_values = e._cid_intensity_values;
            _cid_mz_values = e._cid_mz_values; 
            _CIDPeaks = e._CIDPeaks;
            _cidProfileType = e._cidProfileType;
            _cidScore = e._cidScore;
           // _cid_parent_mz = e._cid_parent_mz; 

            _hcd_scan = e._hcd_scan;
            _hcd_intensity_values = e._hcd_intensity_values;
            _hcd_mz_values = e._hcd_mz_values; 
            _HCDPeaks = e._HCDPeaks;
            _hcdProfileType = e._hcdProfileType;
            _hcdScore = e._hcdScore;
           // _hcd_parent_mz = e._hcd_parent_mz;

            _etd_scan = e._etd_scan;
            _etd_intensity_values = e._etd_intensity_values;
            _etd_mz_values = e._etd_mz_values;
            _ETDPeaks = e._ETDPeaks;
            _etdProfileType = e._etdProfileType;
            _etdScore = e._etdScore;
           // _etd_parent_mz = e._etd_parent_mz; 

            _fdr = e._fdr;
            _cid_sequencing_score = e._cid_sequencing_score;
            _FalseHit = e._FalseHit;

        }

     /*   public double PeptideMass
        {
            get { return _PeptideMass; }
            set { _PeptideMass = value; }
        }
        public double GlycanMass
        {
            get { return _GlycanMass; }
            set { _GlycanMass = value; }
        }

        public string ProteinName
        {
            get { return _ProteinName; }
            set { _ProteinName = value; }
        }
        public string PeptideSeq
        {
            get { return _PeptideSeq; }
            set { _PeptideSeq = value; }
        }
        public string NGlycoSite
        {
            get { return _NGlycoSite; }
            set { _NGlycoSite = value; }
        }

        public string GlycanComposition
        {
            get { return _GlycanComposition; }
            set { _GlycanComposition = value; }
        }
        public double PPM
        {
            get { return _PPM; }
            set { _PPM = value; }
        }
        public double Y1
        {
            get { return _Y1; }
            set { _Y1 = value; }
        }*/

        public GlycopeptideRecord GP_Record
        {
            get { return _gp_record; }
            set { _gp_record = value; }
        }

        public float[] CIDMzs
        {
            get { return _cid_mz_values; }
            set { _cid_mz_values = value; }
        }
        public float[] CIDIntensities
        {
            get { return _cid_intensity_values; }
            set { _cid_intensity_values = value; }
        }

        public GlypID.Peaks.clsPeak[] CIDPeaks
        {
            get { return _CIDPeaks; }
            set { _CIDPeaks = value; }
        }
        public float[] HCDMzs
        {
            get { return _hcd_mz_values; }
            set { _hcd_mz_values = value; }
        }
        public float[] HCDIntensities
        {
            get { return _hcd_intensity_values; }
            set { _hcd_intensity_values = value; }
        }
        public GlypID.Peaks.clsPeak[] HCDPeaks
        {
            get { return _HCDPeaks; }
            set { _HCDPeaks = value; }
        }

        public float[] ETDMzs
        {
            get { return _etd_mz_values; }
            set { _etd_mz_values = value; }
        }
        public float[] ETDIntensities
        {
            get { return _etd_intensity_values; }
            set { _etd_intensity_values = value; }
        }
        public GlypID.Peaks.clsPeak[] ETDPeaks
        {
            get { return _ETDPeaks; }
            set { _ETDPeaks = value; }
        }
        public GlypID.HornTransform.clsHornTransformResults TransformResult
        {
            get { return _TransformResult; }
            set { _TransformResult = value; }
        }
        public GlypID.enmGlycanType GlycanType
        {
            get { return _glycanType; }
            set { _glycanType = value; }
        }
        public double HCDScore
        {
            get { return _hcdScore; }
            set { _hcdScore = value; }
        }
        public double CIDScore
        {
            get { return _cidScore; }
            set { _cidScore = value; }
        }
        public double CIDSequencingScore
        {
            get { return _cid_sequencing_score; }
            set { _cid_sequencing_score = value; }
        }
        public double ETDScore
        {
            get { return _etdScore; }
            set { _etdScore = value; }
        }
        public GlypID.enmProfileType CIDProfileType
        {
            get { return _cidProfileType; }
            set { _cidProfileType = value; }
        }
        public GlypID.enmProfileType HCDProfileType
        {
            get { return _hcdProfileType; }
            set { _hcdProfileType = value; }
        }
        public GlypID.enmProfileType ETDProfileType
        {
            get { return _etdProfileType; }
            set { _etdProfileType = value; }
        }
        public double ParentMz
        {
            get { return _parent_mz; }
            set { _parent_mz = value; }
        }
        public double ParentScanTime
        {
            get { return _parent_scan_time; }
            set { _parent_scan_time = value; }
        }

        public void ClearTransformRecord()
        {
            _parent_mz = 0;
            _parent_scan = 0;             
            TransformResult = new GlypID.HornTransform.clsHornTransformResults();
        }
        public void ClearCID()
        {
            _cid_scan = 0;
            _cidScore = 0;
            _cid_sequencing_score = 0; 
            //_cid_parent_mz = 0; 

            _cid_intensity_values = null; 
            /*if (_cid_mz_values.Length > 0)
                Array.Clear(_cid_mz_values, 0, _cid_mz_values.Length);
            if (_cid_intensity_values.Length > 0)
                Array.Clear(_cid_intensity_values, 0, _cid_intensity_values.Length);
            if (_CIDPeaks.Length > 0)
                Array.Clear(_CIDPeaks, 0, _CIDPeaks.Length);*/
        }
        public void ClearHCD()
        {
            _hcd_scan = 0;
            _hcdScore = 1;
            //_hcd_parent_mz = 0; 
            if (_hcd_intensity_values.Length > 0)
                Array.Clear(_hcd_intensity_values, 0, _hcd_intensity_values.Length);
            if (_hcd_mz_values.Length > 0)
                Array.Clear(_hcd_mz_values, 0, _hcd_mz_values.Length);
            if (_HCDPeaks.Length > 0)
                Array.Clear(_HCDPeaks, 0, _HCDPeaks.Length);
        }
        public void ClearETD()
        {
            _etd_scan = 0;
            _etdScore = 0;
            _fdr = 1.0f; 
           // _etd_parent_mz = 0; 
            if (_etd_intensity_values.Length > 0)
                Array.Clear(_etd_intensity_values, 0, _etd_intensity_values.Length);
            if (_etd_mz_values.Length > 0)
                Array.Clear(_etd_mz_values, 0, _etd_mz_values.Length);
            if (_ETDPeaks.Length > 0)
                Array.Clear(_ETDPeaks, 0, _ETDPeaks.Length);
        }
        public void ClearGPInfo()
        {
            _glycanType = GlypID.enmGlycanType.NA;
            GP_Record = new GlycopeptideRecord(); 
        }
        public bool FalseHit
        {
           get { return _FalseHit ; }
           set {_FalseHit  = value;}
        }
        public int ParentScan
        {
            get { return _parent_scan; }
            set { _parent_scan = value; }
        }
        public int HCDScan
        {
            get { return _hcd_scan; }
            set { _hcd_scan = value; }
        }
        public int ETDScan
        {
            get { return _etd_scan; }
            set { _etd_scan = value; }
        }
        public int CIDScan
        {
            get { return _cid_scan; }
            set { _cid_scan = value; }
        }
        public float FDR
        {
            get { return _fdr; }
            set { _fdr = value; }
        }
        public int ID
        {
            get { return _id; }
            set { _id = value; }
        }

        



    }
}
