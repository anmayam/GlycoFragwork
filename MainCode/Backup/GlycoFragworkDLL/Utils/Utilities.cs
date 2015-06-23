using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GlycoFragworkDLL.Utils
{
    public class Utilities
    {
       public  double _CC_MASS = 1.00727638;// SHOULD CHECK THIS WITH CHuan-yih
       public double _H2O_mass = 18.01056;

        public Utilities()
        {
        }
        public double CalculateDelMassPPM(double mass1, double mass2)
        {
            double ppm = (Math.Abs(mass2 - mass1) / mass1) * 1000000;
            return ppm; 
        }

        public double CalculateMass(double mz, short cs)
        {
            double mass = (mz - _CC_MASS) * cs ; 
            return mass ; 
        }
        
        public double CalculateMz(double mass, short cs)
        {
            double mz = mass/cs + _CC_MASS ; 
            return mz ; 
        }

        public double CalculateDelMassDa(double mass1, double mass2)
        {
            return (Math.Abs(mass2-mass1)) ; 
        }

        public double CalculateGPMass(double seq_mass, double glycan_mass)
        {
            double gpmass = seq_mass + glycan_mass - _H2O_mass;
            return gpmass; 
        }

        
    }
}
