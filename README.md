### STORM analysis
The MATLAB codes for the 1D autocorrelation and cross-correlation analyses of the STORM images of the axon segments can be found in the folder: `./auto-correlation and cross-correlation analyses/`.

`./auto-correlation and cross-correlation analyses auto/autocorr_batch_analysis_bin.m` and `./auto-correlation and cross-correlation analyses auto/autocorr_batch_analysis_text.m`: These two codes calculate the average 1D auto-correlation of the imaged molecules projected to the longitudinal axis of the axon segments, the average spacing of the MPS, and average diameter of axons. The first code takes .bin as inputs, and the second code takes .txt as inputs.

`./auto-correlation and cross-correlation analyses auto/crosscorrelation_batch_analysis.m`: It calculates the average 1D cross-correlation of the imaged molecules projected to the longitudinal axis of the axon segments, the average spacing of the MPS, as well as the histograms of normalized localization counts of the two color STORM signals projected to the longitudinal axis of the axon segments.

### Conventional microscopy analysis
`./conventional microscopy analysis/` contains the MATLAB codes used for calculating results from conventional fluorescence microscopy.
`./conventional microscopy anatlysis/axon_bundling.m`:  it calculates the average axon bundle diameter.
`./conventional microscopy analysis/axon_dendrite_bundling.m`:  it calculates the fraction of dendrite length adhered to axons.
`./conventional microscopy analysis/colocalize_synapse.m`:  it calculates the average synapse density. 
`./conventional microscopy analysis/axon_intensity.m`: it calculates the average fluorescence intensity along axons (for the quantifications of knockdown efficiencies of various molecules and the average surface expression levels of the cell adhesion molecules).
`./conventional microscopy analysis/Dendrite_diameter_vs_distance.m`: it calculates the average diameter of dendrites versus the distance from soma. 

### Data IO and other utilities
`./data_IO/`: contains the MATLAB codes used for converting between different file formats for image processing (`./data_IO /bin2mat.m and ./data_IO /convert_molecule_list_txt2bin.m`), STORM visualization (`./data_IO/bin2png_2c.m` and `./data_IO /bin2png.m`), molecules list writer (`./data_IO/CreateDefaultParameters.m` and `./data_IO/CreateMoleculeList.m`) and neurite tracing across several fields of view (`./data_IO/mosaic_to_matlab.py`). We note that this directory is used to perform tasks such as image processing or format conversion. 

### Two color STORM image processing
`./two color STORM imaging processing/` contains the MATLAB codes used to differentiate the Alexa 647 and CF680 molecules by calculating the intensity ratios of single fluorophores in the two color-channels and to generate a pair of Alexa 647 and CF 680  STORM images.
