using Chemistry;
using EngineLayer;
using EngineLayer.Indexing;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;


namespace GuiFunctions.MetaDraw
{
    public class ProteinCoverage : Protein
    {
        public ProteinCoverage(string sequence, string accession, string organism = null, List<Tuple<string, string>> geneNames = null, IDictionary<int, List<Modification>> oneBasedModifications = null, List<ProteolysisProduct> proteolysisProducts = null, string name = null, string fullName = null, bool isDecoy = false, bool isContaminant = false, List<DatabaseReference> databaseReferences = null, List<SequenceVariation> sequenceVariations = null, List<SequenceVariation> appliedSequenceVariations = null, string sampleNameForVariants = null, List<DisulfideBond> disulfideBonds = null, List<SpliceSite> spliceSites = null, string databaseFilePath = null, bool addTruncations = false) : base(sequence, accession, organism, geneNames, oneBasedModifications, proteolysisProducts, name, fullName, isDecoy, isContaminant, databaseReferences, sequenceVariations, appliedSequenceVariations, sampleNameForVariants, disulfideBonds, spliceSites, databaseFilePath, addTruncations)
        {
        }
    }
}




