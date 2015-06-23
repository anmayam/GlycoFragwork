using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GlycoFragworkDLL.Classes
{
    public class GlycopeptideRecord
    {
        GlypID.Sequence.clsSequence _sequence ;
        GlypID.Glycan.clsGlycan _glycan;

        bool _Is_Decoy; // either peptide or glycan
        string _Glycan_Sequence; 

        double _GP_Mono_Mass;
        double _GP_Average_Mass; 
        double _Glycan_Mono_Mass;
        double _Glycan_Average_Mass; 
        double _Pep_Mono_Mass;
        double _Pep_Average_Mass; 
       
        Utils.Utilities _utils;

        public GlycopeptideRecord()
        {
            _sequence = new GlypID.Sequence.clsSequence();
            _glycan = new GlypID.Glycan.clsGlycan();
            _utils = new GlycoFragworkDLL.Utils.Utilities(); 
        }

        public GlypID.Sequence.clsSequence Sequence
        {
            get { return _sequence; }
            set { _sequence = value; }
        }

        public GlypID.Glycan.clsGlycan Glycan
        {
            get { return _glycan; }
            set { _glycan = value; }
        }

        public string GlycanSequence
        {
            get { return _Glycan_Sequence; }
            set { _Glycan_Sequence = value; }
        }

        public double GP_Mono_Mass
        {
            get { return _GP_Mono_Mass; }
            set { _GP_Mono_Mass = value; }            
        }

        public double GP_Average_Mass
        {
            get { return _GP_Average_Mass; }
            set { _GP_Average_Mass = value; }
        }

        public double GlycanMonoMass
        {
            get { return _Glycan_Mono_Mass; }
            set { _Glycan_Mono_Mass = value; }
        }

        public double GlycanAverageMass
        {
            get { return _Glycan_Average_Mass; }
            set { _Glycan_Average_Mass = value; }
        }

        public double SequenceMonoMass
        {
            get { return _Pep_Mono_Mass; }
            set { _Pep_Mono_Mass = value; }
        }
        public double SequenceAverageMass
        {
            get { return _Pep_Average_Mass; }
            set { _Pep_Average_Mass = value; }
        }
        public bool IsDecoy
        {
            get { return _Is_Decoy; }
            set { _Is_Decoy = value; }
        }
            

       
    }
}
