PCA analysis performed in Bruker software:
    - Lipidomics MALDI TOF image data consolidated by region and condition
    - PCA analysis performed (note that features were isolated to those expected to be significant from co-correlation analysis)
        - As a result, PCA plots are biased (included features are already co-correlated between conditions)
	- Organized by coefficient score (highest variance to lowest)
	- Index pertains to spectral coordinate in MALDI TOF scan.
    - Data exported to .csv (organized by ctrl/ko and brain set - see .py files for indices)
        - Both 12-month datasets combined, only a single dataset for 1-month mice

PCA data was plotted in matplotlib to improve visual elements of the plot (see .py file for library declarations).