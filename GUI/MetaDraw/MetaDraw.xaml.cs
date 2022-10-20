using EngineLayer;
using GuiFunctions;
using Nett;
using OxyPlot;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Data;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MetaDraw.xaml
    /// </summary>
    public partial class MetaDraw : Window
    {
        private ParentChildScanPlotsView itemsControlSampleViewModel;
        private MetaDrawLogic MetaDrawLogic;
        private readonly DataTable propertyView;
        private ObservableCollection<string> plotTypes;
        private ObservableCollection<string> PsmStatPlotFiles;
        public PtmLegendViewModel PtmLegend;
        public ChimeraLegendViewModel ChimeraLegend;
        private ObservableCollection<ModTypeForTreeViewModel> Modifications = new ObservableCollection<ModTypeForTreeViewModel>();
        private static List<string> AcceptedSpectraFormats = new List<string> { ".mzml", ".raw", ".mgf" };
        private static List<string> AcceptedResultsFormats = new List<string> { ".psmtsv", ".tsv" };
        private static List<string> AcceptedSpectralLibraryFormats = new List<string> { ".msp" };
        private SettingsViewModel SettingsView;

        public MetaDraw()
        {
            UsefulProteomicsDatabases.Loaders.LoadElements();

            InitializeComponent();

            InitializeColorSettingsView();
            MetaDrawLogic = new MetaDrawLogic();
            BindingOperations.EnableCollectionSynchronization(MetaDrawLogic.PsmResultFilePaths, MetaDrawLogic.ThreadLocker);
            BindingOperations.EnableCollectionSynchronization(MetaDrawLogic.SpectraFilePaths, MetaDrawLogic.ThreadLocker);
            BindingOperations.EnableCollectionSynchronization(MetaDrawLogic.FilteredListOfPsms, MetaDrawLogic.ThreadLocker);
            BindingOperations.EnableCollectionSynchronization(MetaDrawLogic.PsmsGroupedByFile, MetaDrawLogic.ThreadLocker);

            itemsControlSampleViewModel = new ParentChildScanPlotsView();
            ParentChildScanViewPlots.DataContext = itemsControlSampleViewModel;

            propertyView = new DataTable();
            propertyView.Columns.Add("Name", typeof(string));
            propertyView.Columns.Add("Value", typeof(string));
            dataGridProperties.DataContext = propertyView.DefaultView;

            dataGridScanNums.DataContext = MetaDrawLogic.PeptideSpectralMatchesView;

            Title = "MetaDraw: version " + GlobalVariables.MetaMorpheusVersion;
            base.Closing += this.OnClosing;

            ParentChildScanView.Visibility = Visibility.Collapsed;

            PsmStatPlotFiles = new ObservableCollection<string>();
            selectSourceFileListBox.DataContext = PsmStatPlotFiles;
            plotTypes = new ObservableCollection<string>();
            SetUpPlots();
            plotsListBox.ItemsSource = plotTypes;

            ExportButton.Content = "Export As " + MetaDrawSettings.ExportType;
        }

        private void Window_Drop(object sender, DragEventArgs e)
        {
            string[] files = ((string[])e.Data.GetData(DataFormats.FileDrop)).OrderBy(p => p).ToArray();

            if (files != null)
            {
                foreach (var draggedFilePath in files)
                {
                    if (File.Exists(draggedFilePath))
                    {
                        AddFile(draggedFilePath);
                    }
                }
            }
        }

        private void AddFile(string filePath)
        {
            var theExtension = GlobalVariables.GetFileExtension(filePath).ToLowerInvariant();

            if (AcceptedSpectraFormats.Contains(theExtension))
            {
                if (!MetaDrawLogic.SpectraFilePaths.Contains(filePath))
                {
                    MetaDrawLogic.SpectraFilePaths.Add(filePath);

                    if (MetaDrawLogic.SpectraFilePaths.Count == 1)
                    {
                        spectraFileNameLabel.Text = filePath;
                    }
                    else
                    {
                        spectraFileNameLabel.Text = "[Mouse over to view files]";
                    }

                    spectraFileNameLabel.ToolTip = string.Join("\n", MetaDrawLogic.SpectraFilePaths);
                    resetSpectraFileButton.IsEnabled = true;
                }
            }
            else if (AcceptedResultsFormats.Contains(theExtension))
            {
                if (!MetaDrawLogic.PsmResultFilePaths.Contains(filePath))
                {
                    MetaDrawLogic.PsmResultFilePaths.Add(filePath);

                    if (MetaDrawLogic.PsmResultFilePaths.Count == 1)
                    {
                        psmFileNameLabel.Text = filePath;
                        psmFileNameLabelStat.Text = filePath;
                    }
                    else
                    {
                        psmFileNameLabel.Text = "[Mouse over to view files]";
                        psmFileNameLabelStat.Text = "[Mouse over to view files]";
                    }

                    psmFileNameLabel.ToolTip = string.Join("\n", MetaDrawLogic.PsmResultFilePaths);
                    resetPsmFileButton.IsEnabled = true;

                    psmFileNameLabelStat.ToolTip = string.Join("\n", MetaDrawLogic.PsmResultFilePaths);
                    resetPsmFileButtonStat.IsEnabled = true;
                }
            }
            else if (AcceptedSpectralLibraryFormats.Contains(theExtension))
            {
                // TODO: display this somewhere in the GUI
                if (!MetaDrawLogic.SpectralLibraryPaths.Contains(filePath))
                {
                    MetaDrawLogic.SpectralLibraryPaths.Add(filePath);
                    specLibraryLabel.Text = filePath;
                    specLibraryLabel.ToolTip = string.Join("\n", MetaDrawLogic.SpectralLibraryPaths);
                    resetSpecLibraryButton.IsEnabled = true;
                }
            }
            else
            {
                MessageBox.Show("Cannot read file type: " + theExtension);
            }
        }

        /// <summary>
        /// Event triggers when a different cell is selected in the PSM data grid
        /// </summary>
        private void dataGridScanNums_SelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            if (dataGridScanNums.SelectedItem == null || sender == null)
            {
                ClearPresentationArea();
                return;
            }

            var strings = sender.ToString();

            MetaDrawLogic.CleanUpCurrentlyDisplayedPlots();
            wholeSequenceCoverageHorizontalScroll.ScrollToLeftEnd();
            plotView.Visibility = Visibility.Visible;
            PsmFromTsv psm = (PsmFromTsv)dataGridScanNums.SelectedItem;

            // Chimera plotter
            if (((Grid)MetaDrawTabControl.SelectedContent).Name == "chimeraPlotGrid")
            {
                List<PsmFromTsv> chimericPsms = MetaDrawLogic.FilteredListOfPsms
                    .Where(p => p.Ms2ScanNumber == psm.Ms2ScanNumber && p.FileNameWithoutExtension == psm.FileNameWithoutExtension).ToList();
                MetaDrawLogic.DisplayChimeraSpectra(chimeraPlot, chimericPsms, out List<string> error);
                if (error != null && error.Count > 0)
                    Debugger.Break();
                ClearPresentationArea();
                wholeSequenceCoverageHorizontalScroll.Visibility = Visibility.Collapsed;


                if (MetaDrawSettings.ShowLegend)
                {
                    ChimeraLegend = new(chimericPsms);
                    ChimeraLegendControl.DataContext = ChimeraLegend;
                }
                return;
            }
            else
            {
                wholeSequenceCoverageHorizontalScroll.Visibility = Visibility.Visible;
            }

            SetSequenceDrawingPositionSettings(true);
            // Selection of ambiguous psm => clean up the canvases and show the option box
            if (psm.FullSequence.Contains('|') && sender.ToString() != "System.Object")
            {
                // clear all drawings of the previous non-ambiguous psm
                ClearPresentationArea();

                AmbiguousWarningTextBlocks.Visibility = Visibility.Visible;
                AmbiguousSequenceOptionBox.Visibility = Visibility.Visible;

                // create a psm object for each ambiguous option and add it to the dropdown box
                var fullSeqs = psm.FullSequence.Split('|');
                foreach (var fullSeq in fullSeqs)
                {
                    PsmFromTsv oneAmbiguousPsm = new(psm, fullSeq);
                    AmbiguousSequenceOptionBox.Items.Add(oneAmbiguousPsm);
                }
                return;
            }
            // Psm selected from ambiguous dropdown => adjust the psm to be drawn
            else if (psm.FullSequence.Contains('|') && sender.ToString() == "System.Object")
            {
                psm = (PsmFromTsv)AmbiguousSequenceOptionBox.SelectedItem;
            }
            // Selection of non-ambiguous psm => clear items in the drop down
            else if (!psm.FullSequence.Contains('|'))
            {
                AmbiguousSequenceOptionBox.Items.Clear();
                AmbiguousWarningTextBlocks.Visibility = Visibility.Collapsed;
                AmbiguousSequenceOptionBox.Visibility = Visibility.Collapsed;
                wholeSequenceCoverageHorizontalScroll.Visibility = Visibility.Visible;
                if (psm.BaseSeq.Length > MetaDrawSettings.NumberOfAAOnScreen)
                {
                    GrayBox.Opacity = 0.7;
                }
                else
                {
                    GrayBox.Opacity = 0;
                }
            }

            // display the ion and elements correctly
            MetaDrawSettings.DrawMatchedIons = true;


            // define initial limits for sequence annotation
            double maxDisplayedPerRow = (int)Math.Round((UpperSequenceAnnotaiton.ActualWidth - 10) / MetaDrawSettings.AnnotatedSequenceTextSpacing, 0) + 7;
            MetaDrawSettings.SequenceAnnotationSegmentPerRow = (int)Math.Floor(maxDisplayedPerRow / (double)(MetaDrawSettings.SequenceAnnotaitonResiduesPerSegment + 1));

            // draw the annotated spectrum
            MetaDrawLogic.DisplaySequences(stationarySequenceCanvas, scrollableSequenceCanvas, sequenceAnnotationCanvas, psm);
            MetaDrawLogic.DisplaySpectrumMatch(plotView, psm, itemsControlSampleViewModel, out var errors);

            // add ptm legend if desired
            if (MetaDrawSettings.ShowLegend)
            {
                
                int descriptionLineCount = MetaDrawSettings.SpectrumDescription.Count(p => p.Value);
                descriptionLineCount += (int)Math.Floor((psm.ProteinName.Length - 20) / 26.0);
                if (psm.ProteinAccession.Length > 10)
                    descriptionLineCount++;
                double verticalOffset = descriptionLineCount * 14;
                
                PtmLegend = new PtmLegendViewModel(psm, verticalOffset);
                ChildScanPtmLegendControl.DataContext = PtmLegend;
                SequenceCoveragePtmLegendControl.DataContext = PtmLegend;
            }

            //draw the sequence coverage if not crosslinked
            if (psm.ChildScanMatchedIons == null)
            {
                MetaDrawLogic.DrawSequenceCoverageMap(psm, sequenceText, map); //TODO: figure out how to show coverage on crosslinked peptides
                ParentChildScanView.Visibility = Visibility.Collapsed;
                SequenceCoverageAnnotationView.Visibility = Visibility.Visible;
            }
            else
            {
                ParentChildScanView.Visibility = Visibility.Visible;
                SequenceCoverageAnnotationView.Visibility = Visibility.Collapsed;
            }

            mapViewer.Width = map.Width;

            if (errors != null && errors.Any())
            {
                MessageBox.Show(errors.First());
                return;
            }

            // display PSM properties
            propertyView.Clear();
            System.Reflection.PropertyInfo[] temp = psm.GetType().GetProperties();

            for (int i = 0; i < temp.Length; i++)
            {
                if (temp[i].Name == nameof(psm.MatchedIons))
                {
                    propertyView.Rows.Add(temp[i].Name, string.Join(", ", psm.MatchedIons.Select(p => p.Annotation)));
                }
                else if (temp[i].Name == nameof(psm.VariantCrossingIons))
                {
                    propertyView.Rows.Add(temp[i].Name, string.Join(", ", psm.VariantCrossingIons.Select(p => p.Annotation)));
                }
                else
                {
                    propertyView.Rows.Add(temp[i].Name, temp[i].GetValue(psm, null));
                }
            }
        }

        private void selectSpectraFileButton_Click(object sender, RoutedEventArgs e)
        {
            string filterString = string.Join(";", AcceptedSpectraFormats.Select(p => "*" + p));

            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog
            {
                Filter = "Spectra Files(" + filterString + ")|" + filterString,
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = true
            };
            if (openFileDialog1.ShowDialog() == true)
            {
                foreach (var filePath in openFileDialog1.FileNames.OrderBy(p => p))
                {
                    AddFile(filePath);
                }
            }
        }

        private void selectPsmFileButton_Click(object sender, RoutedEventArgs e)
        {
            string filterString = string.Join(";", AcceptedResultsFormats.Concat(AcceptedSpectralLibraryFormats).Select(p => "*" + p));

            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog
            {
                Filter = "Result Files(" + filterString + ")|" + filterString,
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = false
            };
            if (openFileDialog1.ShowDialog() == true)
            {
                foreach (var filePath in openFileDialog1.FileNames.OrderBy(p => p))
                {
                    AddFile(filePath);
                }
            }
        }

        private void selectSpecLibraryButton_Click(object sender, RoutedEventArgs e)
        {
            string filterString = string.Join(";", AcceptedResultsFormats.Concat(AcceptedSpectralLibraryFormats).Select(p => "*" + p));

            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog
            {
                Filter = "Result Files(" + filterString + ")|" + filterString,
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = false
            };
            if (openFileDialog1.ShowDialog() == true)
            {
                foreach (var filePath in openFileDialog1.FileNames.OrderBy(p => p))
                {
                    AddFile(filePath);
                }
            }
        }

        private void resetFilesButton_Click(object sender, RoutedEventArgs e)
        {
            if (((Button)sender).Name.Equals("resetSpectraFileButton"))
            {
                spectraFileNameLabel.Text = "None Selected";
                MetaDrawLogic.CleanUpSpectraFiles();
            }

            else if (((Button)sender).Name.Equals("resetPsmFileButton"))
            {
                psmFileNameLabel.Text = "None Selected";
                MetaDrawLogic.CleanUpPSMFiles();
            }

            else if (((Button)sender).Name.Equals("resetSpecLibraryButton"))
            {
                specLibraryLabel.Text = "None Selected";
                MetaDrawLogic.CleanUpSpectralLibraryFiles();
            }
            else
            {
                MetaDrawLogic.CleanUpResources();
            }

            // if a psm is selected
            if (MetaDrawLogic.ScrollableSequence != null)
            {
                ClearPresentationArea();
                MetaDrawLogic.FilteredListOfPsms.Clear();
            }
        }

        private void OnClosing(object sender, CancelEventArgs e)
        {
            MetaDrawLogic.CleanUpResources();
        }

        private void settings_Click(object sender, RoutedEventArgs e)
        {
            // save current selected PSM
            var selectedItem = dataGridScanNums.SelectedItem;
            var settingsWindow = new MetaDrawSettingsWindow(SettingsView);
            var result = settingsWindow.ShowDialog();

            ExportButton.Content = "Export As " + MetaDrawSettings.ExportType;
            // re-select selected PSM
            if (result == true)
            {
                // refresh chart
                dataGridScanNums_SelectedCellsChanged(null, null);

                // filter based on new settings
                MetaDrawLogic.FilterPsms();
            }

            if (selectedItem != null)
            {
                dataGridScanNums.SelectedItem = selectedItem;
            }
        }

        private async void loadFilesButton_Click(object sender, RoutedEventArgs e)
        {
            // check for validity
            propertyView.Clear();
            if (!MetaDrawLogic.SpectraFilePaths.Any())
            {
                MessageBox.Show("Please add a spectra file.");
                return;
            }

            if (!MetaDrawLogic.PsmResultFilePaths.Any())
            {
                MessageBox.Show("Please add a search result file.");
                return;
            }

            // load the spectra file
            ToggleButtonsEnabled(false);
 
            prgsFeed.IsOpen = true;
            prgsText.Content = "Loading data...";

            // Add EventHandlers for popup click-in/click-out behaviour
            Deactivated += new EventHandler(prgsFeed_Deactivator);
            Activated += new EventHandler(prgsFeed_Reactivator);

            var slowProcess = Task<List<string>>.Factory.StartNew(() => MetaDrawLogic.LoadFiles(loadSpectra: true, loadPsms: true));
            await slowProcess;
            var errors = slowProcess.Result;

            if (errors.Any())
            {
                string errorList = string.Join("\n", errors);
                MessageBox.Show(errorList);
            }

            PsmStatPlotFiles.Clear();
            foreach (var item in MetaDrawLogic.PsmsGroupedByFile)
            {
                PsmStatPlotFiles.Add(item.Key);
            }

            // done loading - restore controls
            this.prgsFeed.IsOpen = false;

            // Remove added EventHandlers
            Deactivated -= new EventHandler(prgsFeed_Deactivator);
            Activated -= new EventHandler(prgsFeed_Reactivator);

            ToggleButtonsEnabled(true);
        }

        /// <summary>
        /// Deactivates the "loading data" popup if one clicks out of the main window
        /// </summary>
        private void prgsFeed_Deactivator(object sender, EventArgs e)
        {
            prgsFeed.IsOpen = false;
        }

        /// <summary>
        /// Reactivates the "loading data" popup if one clicks into the main window
        /// </summary>
        private void prgsFeed_Reactivator(object sender, EventArgs e)
        {
            prgsFeed.IsOpen = true;
        }

        private void TextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            string txt = (sender as TextBox).Text;
            MetaDrawLogic.FilterPsmsByString(txt);
        }

        /// <summary>
        /// Exports images of the parent and child scan
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void PDFButton_Click(object sender, RoutedEventArgs e)
        {
            if (dataGridScanNums.SelectedCells.Count == 0)
            {
                MessageBox.Show("Please select at least one scan to export");
                return;
            }

            SetSequenceDrawingPositionSettings();
            List<PsmFromTsv> items = new List<PsmFromTsv>();

            foreach (var cell in dataGridScanNums.SelectedItems)
            {
                var psm = (PsmFromTsv)cell;
                items.Add(psm);
            }

            string directoryPath = Path.Combine(Path.GetDirectoryName(MetaDrawLogic.PsmResultFilePaths.First()), "MetaDrawExport",
                    DateTime.Now.ToString("yyyy-MM-dd", CultureInfo.InvariantCulture));


            Canvas legendCanvas = null;
            Vector ptmLegendLocationVector = new();
            List<string> errors = new();
            if (((Grid)MetaDrawTabControl.SelectedContent).Name == "chimeraPlotGrid")
            {
                if (MetaDrawSettings.ShowLegend)
                {
                    ChimeraLegendControl chimeraLegendCopy = new();
                    chimeraLegendCopy.DataContext = ChimeraLegendControl.DataContext;
                    legendCanvas = new();
                    legendCanvas.Children.Add(chimeraLegendCopy);
                    Size legendSize = new Size((int)ChimeraLegendControl.ActualWidth, (int)ChimeraLegendControl.ActualHeight);
                    legendCanvas.Measure(legendSize);
                    legendCanvas.Arrange(new Rect(legendSize));
                    legendCanvas.UpdateLayout();
                }
                MetaDrawLogic.ExportPlot(chimeraPlot, null, items, itemsControlSampleViewModel,
                    directoryPath, out errors, legendCanvas);
            }
            else if (((Grid)MetaDrawTabControl.SelectedContent).Name == "PsmAnnotationGrid")
            {

                if (MetaDrawSettings.ShowLegend)
                {
                    PtmLegendControl ptmLegendCopy = new();
                    ptmLegendCopy.DataContext = ChildScanPtmLegendControl.DataContext;
                    legendCanvas = new();
                    legendCanvas.Children.Add(ptmLegendCopy);
                    Size ptmLegendSize = new Size((int)ChildScanPtmLegendControl.ActualWidth, (int)ChildScanPtmLegendControl.ActualHeight);
                    legendCanvas.Measure(ptmLegendSize);
                    legendCanvas.Arrange(new Rect(ptmLegendSize));
                    legendCanvas.UpdateLayout();
                    ptmLegendLocationVector = (Vector)ChildScanPtmLegendControl.GetType().GetProperty("VisualOffset", BindingFlags.NonPublic | BindingFlags.Instance).GetValue(ChildScanPtmLegendControl);
                    ptmLegendLocationVector.X = PsmAnnotationGrid.ActualWidth - ChildScanPtmLegendControl.ActualWidth;
                }
                MetaDrawLogic.ExportPlot(plotView, stationarySequenceCanvas, items, itemsControlSampleViewModel,
                    directoryPath, out errors, legendCanvas, ptmLegendLocationVector);
            }

            if (errors != null && errors.Any())
            {
                MessageBox.Show(errors.First());
            }
            else
            {
                MessageBox.Show(MetaDrawSettings.ExportType + "(s) exported to: " + directoryPath);
            }
        }

        private void SequenceCoverageExportButton_Click(object sender, RoutedEventArgs e)
        {
            if (dataGridScanNums.SelectedItems.Count == 0 || dataGridScanNums.SelectedItems.Count > 1)
            {
                MessageBox.Show("Please select one psm to export sequence coverage");
                return;
            }

            string directoryPath = Path.Combine(Path.GetDirectoryName(MetaDrawLogic.PsmResultFilePaths.First()), "MetaDrawExport",
                    DateTime.Now.ToString("yyyy-MM-dd", CultureInfo.InvariantCulture));
            PsmFromTsv psm = (PsmFromTsv)dataGridScanNums.SelectedItem;
            MetaDrawLogic.ExportSequenceCoverage(sequenceText, map, directoryPath, psm);
            
            if (Directory.Exists(directoryPath))
            {
                MessageBox.Show(MetaDrawSettings.ExportType + " exported to: " + directoryPath);
            }
        }    

        private void SequenceAnnotationExportButton_Click(object sender, RoutedEventArgs e)
        {
            if (dataGridScanNums.SelectedItems.Count == 0 || dataGridScanNums.SelectedItems.Count > 1)
            {
                MessageBox.Show("Please select one psm to export sequence annotation");
                return;
            }

            string directoryPath = Path.Combine(Path.GetDirectoryName(MetaDrawLogic.PsmResultFilePaths.First()), "MetaDrawExport",
                    DateTime.Now.ToString("yyyy-MM-dd", CultureInfo.InvariantCulture));
            PsmFromTsv psm = (PsmFromTsv)dataGridScanNums.SelectedItem;

            int width = (int)SequenceAnnotationGrid.ActualWidth;
            MetaDrawLogic.ExportAnnotatedSequence(sequenceAnnotationCanvas, SequenceCoveragePtmLegendControl, psm, directoryPath, width);
            
            if (Directory.Exists(directoryPath))
            {
                MessageBox.Show(MetaDrawSettings.ExportType + " exported to: " + directoryPath);
            }
        }

        private void SetUpPlots()
        {
            foreach (var plot in PlotModelStat.PlotNames)
            {
                plotTypes.Add(plot);
            }
        }

        private void loadFilesButtonStat_Click(object sender, RoutedEventArgs e)
        {
            // check for validity
            if (!MetaDrawLogic.PsmResultFilePaths.Any())
            {
                MessageBox.Show("Please add a search result file.");
                return;
            }

            (sender as Button).IsEnabled = false;
            selectPsmFileButtonStat.IsEnabled = false;
            resetPsmFileButtonStat.IsEnabled = false;
            prgsFeedStat.IsOpen = true;

            // load the PSMs
            this.prgsTextStat.Content = "Loading data...";
            MetaDrawLogic.LoadFiles(loadSpectra: false, loadPsms: true);

            PsmStatPlotFiles.Clear();
            foreach (var item in MetaDrawLogic.PsmsGroupedByFile)
            {
                PsmStatPlotFiles.Add(item.Key);
            }

            // done loading - restore controls
            this.prgsFeedStat.IsOpen = false;
            (sender as Button).IsEnabled = true;
            selectPsmFileButtonStat.IsEnabled = true;
            resetPsmFileButtonStat.IsEnabled = true;
        }

        /// <summary>
        /// Export for Data Visualization tab
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void CreatePlotPdf_Click(object sender, RoutedEventArgs e)
        {
            var selectedItem = plotsListBox.SelectedItem;

            if (selectedItem == null)
            {
                MessageBox.Show("Select a plot type to export!");
                return;
            }

            if (!MetaDrawLogic.PsmResultFilePaths.Any())
            {
                MessageBox.Show("No PSMs are loaded!");
                return;
            }

            if (selectSourceFileListBox.SelectedItems.Count == 0)
            {
                MessageBox.Show("Please select a source file.");
                return;
            }

            var plotName = selectedItem as string;
            var fileDirectory = Path.Combine(Path.GetDirectoryName(MetaDrawLogic.PsmResultFilePaths.First()), "MetaDrawExport",
                    DateTime.Now.ToString("yyyy-MM-dd", CultureInfo.InvariantCulture));
            var fileName = String.Concat(plotName, ".pdf");

            // update font sizes to exported PDF's size
            double tmpW = plotViewStat.Width;
            double tmpH = plotViewStat.Height;
            plotViewStat.Width = 1000;
            plotViewStat.Height = 700;
            plotViewStat.UpdateLayout();
            PlotViewStat_SizeChanged(plotViewStat, null);

            if (!Directory.Exists(fileDirectory))
            {
                Directory.CreateDirectory(fileDirectory);
            }

            using (Stream writePDF = File.Create(Path.Combine(fileDirectory, fileName)))
            {
                PdfExporter.Export(plotViewStat.Model, writePDF, 1000, 700);
            }
            plotViewStat.Width = tmpW;
            plotViewStat.Height = tmpH;
            MessageBox.Show(MetaDrawSettings.ExportType + " Created at " + Path.Combine(fileDirectory, fileName) + "!");
        }

        private async void PlotSelected(object sender, SelectionChangedEventArgs e)
        {
            var listview = sender as ListView;
            var plotName = listview.SelectedItem as string;

            if (MetaDrawLogic.FilteredListOfPsms.Count == 0)
            {
                MessageBox.Show("There are no PSMs to analyze.\n\nLoad the current file or choose a new file.");
                return;
            }
            if (selectSourceFileListBox.SelectedItems.Count == 0)
            {
                MessageBox.Show("Please select a source file.");
                return;
            }

            // get psms from selected source files
            ObservableCollection<PsmFromTsv> psms = new ObservableCollection<PsmFromTsv>();
            Dictionary<string, ObservableCollection<PsmFromTsv>> psmsBSF = new Dictionary<string, ObservableCollection<PsmFromTsv>>();
            foreach (string fileName in selectSourceFileListBox.SelectedItems)
            {
                psmsBSF.Add(fileName, MetaDrawLogic.PsmsGroupedByFile[fileName]);
                foreach (PsmFromTsv psm in MetaDrawLogic.PsmsGroupedByFile[fileName])
                {
                    psms.Add(psm);
                }
            }
            PlotModelStat plot = await Task.Run(() => new PlotModelStat(plotName, psms, psmsBSF));
            plotViewStat.DataContext = plot;
            PlotViewStat_SizeChanged(plotViewStat, null);
        }

        private void selectSourceFileListBox_SelectionChanged(object sender, EventArgs e)
        {
            // refreshes the plot using the new source file
            if (plotsListBox.SelectedIndex > -1 && selectSourceFileListBox.SelectedItems.Count != 0)
            {
                PlotSelected(plotsListBox, null);
            }
        }

        private void selectAllSourceFiles_Click(object sender, RoutedEventArgs e)
        {
            selectSourceFileListBox.SelectAll();
        }

        private void deselectAllSourceFiles_Click(object sender, RoutedEventArgs e)
        {
            selectSourceFileListBox.SelectedIndex = -1;
        }

        // scales the font size down for the x axis labels of the PTM histogram when the window gets too small
        private void PlotViewStat_SizeChanged(object sender, SizeChangedEventArgs e)
        {
            if (plotsListBox.SelectedItem == null || !plotsListBox.SelectedItem.ToString().Equals("Histogram of PTM Spectral Counts"))
            {
                return;
            }
            OxyPlot.Wpf.PlotView plot = sender as OxyPlot.Wpf.PlotView;
            if (plot != null && plot.Model != null)
            {
                plot.Model.DefaultXAxis.TitleFontSize = plot.Model.DefaultFontSize; // stops the title from being scaled
                int count = (int)plot.Model.DefaultXAxis.ActualMaximum;
                int widthCountRatio = 23;   // maintains this width:number of PTM types ratio
                if (plot.ActualWidth / count < widthCountRatio)
                {
                    plot.Model.DefaultXAxis.FontSize = plot.Model.DefaultFontSize * (plot.ActualWidth / (count * widthCountRatio));
                }
                else
                {
                    plot.Model.DefaultXAxis.FontSize = plot.Model.DefaultFontSize;
                }
            }
        }

        /// <summary>
        /// Redraws the Stationary Sequence whenever the scrolling sequence is scrolled
        /// </summary>
        /// <param name="sender"> 
        /// <param name="e"></param>
        private void wholeSequenceCoverageHorizontalScroll_Scroll(object sender, ScrollChangedEventArgs e)
        {
            PsmFromTsv psm = (PsmFromTsv)dataGridScanNums.SelectedItem;
            if (AmbiguousSequenceOptionBox.Items.Count > 1 && AmbiguousSequenceOptionBox.SelectedItem != null)
            {
                psm = (PsmFromTsv)AmbiguousSequenceOptionBox.SelectedItem;

                // Draw the matched ions for the first ambiguous sequence only
                if (AmbiguousSequenceOptionBox.SelectedIndex == 0)
                {
                    MetaDrawSettings.DrawMatchedIons = true;
                }
                else
                {
                    MetaDrawSettings.DrawMatchedIons = false;
                }
            }
            SetSequenceDrawingPositionSettings();
            if (MetaDrawLogic.StationarySequence != null && !psm.FullSequence.Contains('|'))
            {
                DrawnSequence.DrawStationarySequence(psm, MetaDrawLogic.StationarySequence, 10);
            }
                
        }

        /// <summary>
        /// Redraws the Stationary and Scrollable sequences upon the selection of an Ambiguous PSM
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void AmbiguousSequenceOptionBox_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            if (AmbiguousWarningTextBlocks.Visibility != Visibility.Collapsed)
            {
                AmbiguousWarningTextBlocks.Visibility = Visibility.Collapsed;
            }

            wholeSequenceCoverageHorizontalScroll.ScrollToLeftEnd();
            PsmFromTsv psm = (PsmFromTsv)dataGridScanNums.SelectedItem;
            if (AmbiguousSequenceOptionBox.Items.Count > 1 && AmbiguousSequenceOptionBox.SelectedItem != null)
            {
                psm = (PsmFromTsv)AmbiguousSequenceOptionBox.SelectedItem;
                wholeSequenceCoverageHorizontalScroll.Visibility = Visibility.Visible;

                // Draw the matched ions for the first ambiguous sequence only
                if (AmbiguousSequenceOptionBox.SelectedIndex == 0)
                {
                    MetaDrawSettings.DrawMatchedIons = true;
                }
                else
                {
                    MetaDrawSettings.DrawMatchedIons = false;
                }
            }
            SetSequenceDrawingPositionSettings(true);
            object obj = new object();
            if (AmbiguousSequenceOptionBox.Items.Count > 0)
                dataGridScanNums_SelectedCellsChanged(obj, null);
            MetaDrawLogic.DisplaySequences(stationarySequenceCanvas, scrollableSequenceCanvas, sequenceAnnotationCanvas, psm);
        }

        /// <summary>
        /// Method to set the MetaDrawSettings fields FirstAAOnScreen and NumberofAAonScreen to the current scrolling sequence position
        /// </summary>
        private void SetSequenceDrawingPositionSettings(bool reset = false)
        {
            if (dataGridScanNums.SelectedItem == null)
                return;
            double width = SequenceAnnotationArea.ActualWidth;
            double offset = wholeSequenceCoverageHorizontalScroll.HorizontalOffset;
            if (reset)
            {
                offset = 0;
            }
            PsmFromTsv psm = (PsmFromTsv)dataGridScanNums.SelectedItem;
            if (AmbiguousSequenceOptionBox.Items.Count > 1 && AmbiguousSequenceOptionBox.SelectedItem != null)
            {
                psm = (PsmFromTsv)AmbiguousSequenceOptionBox.SelectedItem;
            }

            int lettersOnScreen = (int)Math.Round((width - 10) / MetaDrawSettings.AnnotatedSequenceTextSpacing, 0);
            int firstLetterOnScreen = (int)Math.Round((offset) / MetaDrawSettings.AnnotatedSequenceTextSpacing, 0);
            if ((firstLetterOnScreen + lettersOnScreen) > psm.BaseSeq.Length)
            {
                lettersOnScreen = psm.BaseSeq.Length - firstLetterOnScreen;
            }
            MetaDrawSettings.FirstAAonScreenIndex = firstLetterOnScreen;
            MetaDrawSettings.NumberOfAAOnScreen = lettersOnScreen;
        }

        /// <summary>
        /// Allows the color settings to load asynchronously, avoiding a minor delay in MetaDraw Launch
        /// </summary>
        private async void InitializeColorSettingsView()
        {
            SettingsViewModel view = new SettingsViewModel();
            await view.Initialization;
            SettingsView = view;
        }

        /// <summary>
        /// Fires the command in the PtmLegend to decrease the residues per segment by one
        /// </summary>
        /// <param name="sender">Button Object</param>
        /// <param name="e"></param>
        private void residuesPerSegmentcmdDown_Click(object sender, RoutedEventArgs e)
        {
            if (PtmLegend.ResiduesPerSegment == 1)
            {
                MessageBox.Show("Value must be greater than 0");
                return;
            }

            PtmLegend.DecreaseResiduesPerSegment();
            PsmFromTsv psm = (PsmFromTsv)dataGridScanNums.SelectedItem;
            MetaDrawLogic.DisplaySequences(null, null, sequenceAnnotationCanvas, psm);
        }

        /// <summary>
        /// Fires the command in the PtmLegend to increase the residues per segment of the annotated sequence by one
        /// </summary>
        /// <param name="sender">Button Object</param>
        /// <param name="e"></param>
        private void residuesPerSegmentcmdUp_Click(object sender, RoutedEventArgs e)
        {
            PtmLegend.IncreaseResiduesPerSegment();
            PsmFromTsv psm = (PsmFromTsv)dataGridScanNums.SelectedItem;
            MetaDrawLogic.DisplaySequences(null, null, sequenceAnnotationCanvas, psm);
        }

        /// <summary>
        /// Fires the command in the PtmLegend to decrease the segments per row of the annotated sequence by one
        /// </summary>
        /// <param name="sender">Button Object</param>
        /// <param name="e"></param>
        private void segmentsPerRowcmdDown_Click(object sender, RoutedEventArgs e)
        {
            if (PtmLegend.SegmentsPerRow == 1)
            {
                MessageBox.Show("Value must be greater than 0");
                return;
            }

            PtmLegend.DecreaseSegmentsPerRow();
            PsmFromTsv psm = (PsmFromTsv)dataGridScanNums.SelectedItem;
            MetaDrawLogic.DisplaySequences(null, null, sequenceAnnotationCanvas, psm);
        }

        /// <summary>
        /// Fires the command in the PtmLegend to increase the segments per row of the annotated sequence by one
        /// </summary>
        /// <param name="sender">Button Object</param>
        /// <param name="e"></param>
        private void segmentsPerRowcmdUp_Click(object sender, RoutedEventArgs e)
        {
            PtmLegend.IncreaseSegmentsPerRow();
            PsmFromTsv psm = (PsmFromTsv)dataGridScanNums.SelectedItem;
            MetaDrawLogic.DisplaySequences(null, null, sequenceAnnotationCanvas, psm);
        }


        private void MetaDrawTabControl_OnSelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            PsmFromTsv selectedPsm = (PsmFromTsv)dataGridScanNums.SelectedItem;

            // switch from chimera to other views
            if (e.RemovedItems.Count > 0 && ((TabItem)e.RemovedItems[0]).Name == "ChimeraScanPlot")
            {
                MetaDrawLogic.FilterPsms();

                // reselect what was selected
                if (selectedPsm != null && MetaDrawLogic.FilteredListOfPsms.Contains(selectedPsm))
                {
                    int psmIndex = MetaDrawLogic.FilteredListOfPsms.IndexOf(selectedPsm);
                    dataGridScanNums.SelectedIndex = psmIndex;
                    dataGridScanNums_SelectedCellsChanged(new object(), null);
                }
            }

            // switch from other view to chimera
            if (e.AddedItems.Count > 0 && ((TabItem)e.AddedItems[0]).Name == "ChimeraScanPlot")
            {
                MetaDrawLogic.FilterPsmsToChimerasOnly();

                // reselect what was selected
                if (selectedPsm == null || !MetaDrawLogic.FilteredListOfPsms.Contains(selectedPsm)) return;
                int psmIndex = MetaDrawLogic.FilteredListOfPsms.IndexOf(selectedPsm);
                dataGridScanNums.SelectedIndex = psmIndex;
                dataGridScanNums_SelectedCellsChanged(new object(), null);
            }
        }

        /// <summary>
        /// Clears and resets the presentation area
        /// </summary>
        private void ClearPresentationArea()
        {
            DrawnSequence.ClearCanvas(scrollableSequenceCanvas);
            DrawnSequence.ClearCanvas(stationarySequenceCanvas);
            DrawnSequence.ClearCanvas(map);
            DrawnSequence.ClearCanvas(sequenceText);
            DrawnSequence.ClearCanvas(sequenceAnnotationCanvas);
            plotView.Visibility = Visibility.Hidden;
            GrayBox.Opacity = 0;
            wholeSequenceCoverageHorizontalScroll.Visibility = Visibility.Collapsed;
            AmbiguousSequenceOptionBox.Items.Clear();
            plotView.Visibility = Visibility.Hidden;

            if (ChimeraLegend != null)
                ChimeraLegend.Visibility = false;
            if (PtmLegend != null)
                PtmLegend.Visibility = false;
        }

        /// <summary>
        /// Enables and disables the buttons on the main MetaDraw view
        /// </summary>
        /// <param name="value">true = enabled, false = disable</param>
        private void ToggleButtonsEnabled(bool value)
        {
            loadFiles.IsEnabled = value;
            selectSpectraFileButton.IsEnabled = value;
            selectPsmFileButton.IsEnabled = value;
            selectSpecLibraryButton.IsEnabled = value;
            resetPsmFileButton.IsEnabled = value;
            resetSpectraFileButton.IsEnabled = value;
            resetSpectraFileButton.IsEnabled = value;
            ExportButton.IsEnabled = value;
        }
    }
}