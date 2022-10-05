using Chemistry;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text.RegularExpressions;
using EngineLayer.GlycoSearch;
using System.IO;

namespace EngineLayer
{
    public class PeptideFromTsv
    {
        public string FullSequence { get; }
        public string FileNameWithoutExtension { get; }
        public string ProteinAccession { get; }
        public string BaseSeq { get; }
        public string PeptideMonoMass { get; }
        public string ProteinName { get; }
        public string GeneName { get; }
        public string OrganismName { get; }
        
        // First amino acid in protein is amino acid number 1, which differs from internal code numbering with N-terminus as 1
        // This numbering is for the peptide location within the protein
        public string StartAndEndResiduesInProtein { get; }
        public string PreviousAminoAcid { get; }
        public string NextAminoAcid { get; }
        public string DecoyContamTarget { get; }

        public PeptideFromTsv(string line, char[] split, Dictionary<string, int> parsedHeader)
        {
            var spl = line.Split(split).Select(p => p.Trim('\"')).ToArray();

            //Required properties
            FileNameWithoutExtension = spl[parsedHeader[PeptideTsvHeader.FileName]].Trim();

            // remove file format, e.g., .raw, .mzML, .mgf
            // this is more robust but slower than Path.GetFileNameWithoutExtension
            if (FileNameWithoutExtension.Contains('.'))
            {
                foreach (var knownSpectraFileExtension in GlobalVariables.AcceptedSpectraFormats)
                {
                    FileNameWithoutExtension = Path.GetFileName(FileNameWithoutExtension.Replace(
                        knownSpectraFileExtension, string.Empty, StringComparison.InvariantCultureIgnoreCase));
                }
            }


            BaseSeq = RemoveParentheses(spl[parsedHeader[PeptideTsvHeader.BaseSequence]].Trim());
            FullSequence = spl[parsedHeader[PeptideTsvHeader.FullSequence]];
            PeptideMonoMass = spl[parsedHeader[PeptideTsvHeader.PeptideMonoMass]].Trim();
            DecoyContamTarget = spl[parsedHeader[PeptideTsvHeader.DecoyContaminantTarget]].Trim();

            //For general peptides
           
            ProteinAccession = (parsedHeader[PsmTPeptideTsvHeadersvHeader.ProteinAccession] < 0)
                ? null
                : spl[parsedHeader[PeptideTsvHeader.ProteinAccession]].Trim();
            ProteinName = (parsedHeader[PeptideTsvHeader.ProteinName] < 0)
                ? null
                : spl[parsedHeader[PeptideTsvHeader.ProteinName]].Trim();
            GeneName = (parsedHeader[PeptideTsvHeader.GeneName] < 0)
                ? null
                : spl[parsedHeader[PeptideTsvHeader.GeneName]].Trim();
            OrganismName = (parsedHeader[PeptideTsvHeader.OrganismName] < 0)
                ? null
                : spl[parsedHeader[PeptideTsvHeader.OrganismName]].Trim();
            StartAndEndResiduesInProtein = (parsedHeader[PeptideTsvHeader.StartAndEndResiduesInProtein] < 0)
                ? null
                : spl[parsedHeader[PeptideTsvHeader.StartAndEndResiduesInProtein]].Trim();
            PreviousAminoAcid = (parsedHeader[PeptideTsvHeader.PreviousAminoAcid] < 0)
                ? null
                : spl[parsedHeader[PeptideTsvHeader.PreviousAminoAcid]].Trim();
            NextAminoAcid = (parsedHeader[PeptideTsvHeader.NextAminoAcid] < 0)
                ? null
                : spl[parsedHeader[PeptideTsvHeader.NextAminoAcid]].Trim();

        }

        /// <summary>
        /// Constructor used to disambiguate PeptideFromTsv to a single psm object
        /// </summary>
        /// <param name="psm">psm to disambiguate</param>
        /// <param name="fullSequence">sequence of ambiguous psm to use</param>
        public PeptideFromTsv(PeptideFromTsv peptide, string fullSequence, int index = 0, string baseSequence = "")
        {
            // psm is not ambiguous
            if (!peptide.FullSequence.Contains("|"))
            {
                FullSequence = fullSequence;
                BaseSeq = baseSequence == "" ? peptide.BaseSeq : baseSequence;
                StartAndEndResiduesInProtein = peptide.StartAndEndResiduesInProtein;
                ProteinAccession = peptide.ProteinAccession;
                ProteinName = peptide.ProteinName;
                GeneName = peptide.GeneName;
                PeptideMonoMass = peptide.PeptideMonoMass;
            }
            // potentially ambiguous fields
            else
            {
                FullSequence = fullSequence;
                BaseSeq = baseSequence == "" ? peptide.BaseSeq.Split("|")[index] : baseSequence;
                StartAndEndResiduesInProtein = peptide.StartAndEndResiduesInProtein.Split("|")[index];
                ProteinAccession = peptide.ProteinAccession.Split("|")[index];
                ProteinName = peptide.ProteinName.Split("|")[index];
                GeneName = peptide.GeneName.Split("|")[index];

                if (peptide.PeptideMonoMass.Split("|").Count() == 1)
                {
                    PeptideMonoMass = peptide.PeptideMonoMass.Split("|")[0];
                }
                else
                {
                    PeptideMonoMass = peptide.PeptideMonoMass.Split("|")[index];
                }
            }

            // non ambiguous fields
            FileNameWithoutExtension = peptide.FileNameWithoutExtension;
            OrganismName = peptide.OrganismName;
            PreviousAminoAcid = peptide.PreviousAminoAcid;
            NextAminoAcid = peptide.NextAminoAcid;
            DecoyContamTarget = peptide.DecoyContamTarget;
        }

        //Used to remove Silac labels for proper annotation
        public static string RemoveParentheses(string baseSequence)
        {
            if (baseSequence.Contains("("))
            {
                string updatedBaseSequence = "";
                bool withinParentheses = false;
                foreach (char c in baseSequence)
                {
                    if (c == ')') //leaving the parentheses
                    {
                        withinParentheses = false;
                    }
                    else if (c == '(') //entering the parentheses
                    {
                        withinParentheses = true;
                    }
                    else if (!withinParentheses) //if outside the parentheses, preserve this amino acid
                    {
                        updatedBaseSequence += c;
                    }
                    //else do nothing
                }
                return updatedBaseSequence;
            }
            return baseSequence;
        }

        /// <summary>
        /// Parses the full sequence to identify mods
        /// </summary>
        /// <param name="fullSequence"> Full sequence of the peptide in question</param>
        /// <returns> Dictionary with the key being the amino acid position of the mod and the value being the string representing the mod</returns>
        public static Dictionary<int, List<string>> ParseModifications(string fullSeq)
        {
            // use a regex to get all modifications
            string pattern = @"\[(.+?)\]";
            Regex regex = new(pattern);

            // remove each match after adding to the dict. Otherwise, getting positions
            // of the modifications will be rather difficult.
            //int patternMatches = regex.Matches(fullSeq).Count;
            Dictionary<int, List<string>> modDict = new();
            
            RemoveSpecialCharacters(ref fullSeq);
            MatchCollection matches = regex.Matches(fullSeq);
            int currentPosition = 0;
            foreach (Match match in matches)
            {
                GroupCollection group = match.Groups;
                string val = group[1].Value;
                int startIndex = group[0].Index;
                int captureLength = group[0].Length;
                int position = group["(.+?)"].Index;

                List<string> modList = new List<string>();
                modList.Add(val);
                // check to see if key already exist
                // if there is a missed cleavage, then there will be a label on K and a Label on X modification.
                // And, it'll be like [label]|[label] which complicates the positional stuff a little bit.
                // if the already key exists, update the current position with the capture length + 1.
                // otherwise, add the modification to the dict.

                // int to add is startIndex - current position
                int positionToAddToDict = startIndex - currentPosition;
                if (modDict.ContainsKey(positionToAddToDict))
                {
                    modDict[positionToAddToDict].Add(val);
                }
                else
                {
                    modDict.Add(positionToAddToDict, modList);
                }
                currentPosition += startIndex + captureLength;
            }
            return modDict;
        }

        /// <summary>
        /// Fixes an issue where the | appears and throws off the numbering if there are multiple mods on a single amino acid.
        /// </summary>
        /// <param name="fullSeq"></param>
        /// <param name="replacement"></param>
        /// <param name="specialCharacter"></param>
        /// <returns></returns>
        public static void RemoveSpecialCharacters(ref string fullSeq, string replacement = @"", string specialCharacter = @"\|")
        {
            // next regex is used in the event that multiple modifications are on a missed cleavage Lysine (K)
            Regex regexSpecialChar = new(specialCharacter);
            fullSeq = regexSpecialChar.Replace(fullSeq, replacement);
        }


        public override string ToString()
        {
            return FullSequence;
        }
    }
}