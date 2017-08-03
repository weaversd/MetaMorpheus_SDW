﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using Proteomics;

namespace EngineLayer
{
    [Serializable]
    public abstract class CompactPeptideBase
    {

        #region Protected Fields

        protected static readonly double nitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;
        protected static readonly double oxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        protected static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;
        protected static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

        #endregion Protected Fields

        #region Public Fields

        public double[] CTerminalMasses { get; protected set; }
        public double[] NTerminalMasses { get; protected set; }
        public double MonoisotopicMassIncludingFixedMods { get; protected set; }

        #endregion Public Fields

        public override bool Equals(object obj)
        {
            var cp = obj as CompactPeptideBase;
            if (cp == null)
                return false;
            return (
                ((double.IsNaN(MonoisotopicMassIncludingFixedMods) && double.IsNaN(cp.MonoisotopicMassIncludingFixedMods)) || Math.Abs(MonoisotopicMassIncludingFixedMods - cp.MonoisotopicMassIncludingFixedMods) < 1e-7)
                && CTerminalMasses.SequenceEqual(cp.CTerminalMasses)
                && NTerminalMasses.SequenceEqual(cp.NTerminalMasses)
                );
        }

        public override int GetHashCode()
        {
            unchecked
            {
                var result = 0;
                foreach (double b in CTerminalMasses)
                    result = (result * 31) ^ b.GetHashCode();
                return result;
            }
        }

        #region Public Methods

        public double[] ProductMassesMightHaveDuplicatesAndNaNs(List<ProductType> productTypes)
        {
            int massLen = 0;
            bool containsAdot = productTypes.Contains(ProductType.Adot);
            bool containsB = productTypes.Contains(ProductType.B);
            bool containsBnoB1 = productTypes.Contains(ProductType.BnoB1ions);
            bool containsC = productTypes.Contains(ProductType.C);
            bool containsX = productTypes.Contains(ProductType.X);
            bool containsY = productTypes.Contains(ProductType.Y);
            bool containsZdot = productTypes.Contains(ProductType.Zdot);

            if (containsAdot)
                throw new NotImplementedException();
            if (containsBnoB1)
                massLen += NTerminalMasses.Length - 1;
            if (containsB)
                massLen += NTerminalMasses.Length;
            if (containsC)
                massLen += NTerminalMasses.Length;
            if (containsX)
                throw new NotImplementedException();
            if (containsY)
                massLen += CTerminalMasses.Length;
            if (containsZdot)
                massLen += CTerminalMasses.Length;

            double[] massesToReturn = new double[massLen];

            int i = 0;
            for (int j = 0; j < NTerminalMasses.Length; j++)
            {
                var hm = NTerminalMasses[j];
                if (containsBnoB1)
                {
                    if (j > 0)
                    {
                        massesToReturn[i] = hm;
                        i++;
                    }
                }
                if (containsB)
                {
                    massesToReturn[i] = hm;
                    i++;
                }
                if (containsC)
                {
                    massesToReturn[i] = hm + nitrogenAtomMonoisotopicMass + 3 * hydrogenAtomMonoisotopicMass;
                    i++;
                }
            }
            for (int j = 0; j < CTerminalMasses.Length; j++)
            {
                var hm = CTerminalMasses[j];
                if (containsY)
                {
                    massesToReturn[i] = hm + waterMonoisotopicMass;
                    i++;
                }
                if (containsZdot)
                {
                    massesToReturn[i] = hm + oxygenAtomMonoisotopicMass - nitrogenAtomMonoisotopicMass;
                    i++;
                }
            }
            return massesToReturn;
        }

        #endregion Public Methods

        #region Protected Methods

        protected IEnumerable<double> ComputeFollowingFragmentMasses(PeptideWithSetModifications yyy, double prevMass, int oneBasedIndexToLookAt, int direction)
        {
            ModificationWithMass residue_variable_mod = null;
            do
            {
                prevMass += Residue.ResidueMonoisotopicMass[yyy[oneBasedIndexToLookAt - 1]];

                yyy.allModsOneIsNterminus.TryGetValue(oneBasedIndexToLookAt + 1, out residue_variable_mod);
                if (residue_variable_mod == null)
                {
                    yield return prevMass;
                }
                else if (residue_variable_mod.neutralLosses.Count == 1)
                {
                    prevMass += residue_variable_mod.monoisotopicMass - residue_variable_mod.neutralLosses.First();
                    yield return prevMass;
                }
                else
                {
                    foreach (double nl in residue_variable_mod.neutralLosses)
                    {
                        var theMass = prevMass + residue_variable_mod.monoisotopicMass - nl;
                        yield return theMass;
                        if ((direction == 1 && oneBasedIndexToLookAt + direction < yyy.Length) ||
                            (direction == -1 && oneBasedIndexToLookAt + direction > 1))
                            foreach (var nextMass in ComputeFollowingFragmentMasses(yyy, theMass, oneBasedIndexToLookAt + direction, direction))
                                yield return nextMass;
                    }
                    break;
                }
                oneBasedIndexToLookAt += direction;
            } while ((oneBasedIndexToLookAt > 1 && direction == -1) || (oneBasedIndexToLookAt < yyy.Length && direction == 1));
        }

        #endregion Protected Methods
    }
}
