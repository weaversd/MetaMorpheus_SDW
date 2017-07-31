﻿using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public class BinGenerationTest
    {

        #region Public Methods

        [Test]
        public static void TestBinGeneration()
        {
            List<MassDiffAcceptor> massDiffAcceptors = new List<MassDiffAcceptor> { new OpenSearchMode() };
            SearchTask st = new SearchTask()
            {
                DoHistogramAnalysis = true,
                MassDiffAcceptors = massDiffAcceptors,
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain,
                SearchDecoy = false
            };

            string proteinDbFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "BinGenerationTest.xml");
            string mzmlFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "BinGenerationTest.mzml");

            Protein prot1 = new Protein("MEDEEK", "prot1");
            Protein prot2 = new Protein("MENEEK", "prot2");

            ModificationMotif motif;
            ModificationMotif.TryGetMotif("D", out motif);
            ModificationWithMass mod = new ModificationWithMass(null, null, motif, ModificationSites.Any, 10, null, new List<double> { 0 }, null, null);

            var possMod1 = prot1.Digest(st.Protease, st.MaxMissedCleavages, st.MinPeptideLength, st.MaxPeptideLength, st.InitiatorMethionineBehavior, new List<ModificationWithMass>()).First();
            var pep1_0 = possMod1.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4096, 3).First();
            var pep1_10 = possMod1.GetPeptidesWithSetModifications(new List<ModificationWithMass> { mod }, 4096, 3).Last();

            Protein prot3 = new Protein("MAAADAAAAAAAAAAAAAAA", "prot3");

            var possMod3 = prot3.Digest(st.Protease, st.MaxMissedCleavages, st.MinPeptideLength, st.MaxPeptideLength, st.InitiatorMethionineBehavior, new List<ModificationWithMass>()).First();
            var pep2_0 = possMod3.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4096, 3).First();
            var pep2_10 = possMod3.GetPeptidesWithSetModifications(new List<ModificationWithMass>() { mod }, 4096, 3).Last();

            Protein prot4 = new Protein("MNNDNNNN", "prot4");
            var possMod4 = prot4.Digest(st.Protease, st.MaxMissedCleavages, st.MinPeptideLength, st.MaxPeptideLength, st.InitiatorMethionineBehavior, new List<ModificationWithMass>()).First();
            var pep3_10 = possMod4.GetPeptidesWithSetModifications(new List<ModificationWithMass>() { mod }, 4096, 3).Last();

            List<PeptideWithSetModifications> pepsWithSetMods = new List<PeptideWithSetModifications> { pep1_0, pep1_10, pep2_0, pep2_10, pep3_10 };
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(pepsWithSetMods);

            List<Protein> proteinList = new List<Protein> { prot1, prot2, prot3, prot4 };

            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlFilePath, false);
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinList, proteinDbFilePath);

            st.RunTask(
                TestContext.CurrentContext.TestDirectory,
                new List<DbForTask> { new DbForTask(proteinDbFilePath, false) },
                new List<string> { mzmlFilePath },
                null);

            Assert.AreEqual(3, File.ReadLines(Path.Combine(TestContext.CurrentContext.TestDirectory, @"aggregate_OpenSearch.mytsv")).Count());
        }

        #endregion Public Methods

    }
}