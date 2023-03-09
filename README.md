# bios785_team_project

1. The paper is "Single-cell sequencing shows cellular heterogeneity of cutaneous lesions in lupus erythematosus" and can be accessed at [here](https://www.nature.com/articles/s41467-022-35209-1#Abs1
)
2. The GitHub repo is [here](https://github.com/zml314/skin/tree/v1.0.1)

## Data Import
1. The raw data file can be downloaded at [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179633).
2. The disease states are HC, DLE and SLE. The locations are epidermis (E) and Dermis (D). The number of samples are  
||HC|DLE|SLE|
|--|--|--|--|
|E|4|5|5|
|D|4|5|7|
3. Each sample has three files, which are barcodes, features, and matrix. The name of files are organized as "sampleID + (DLE, SLE, or HC) + subjectnumber by diseas state + (D or E)".
4. The data for HC\_2E and HC\_5E sameples are damaged to loss some rows, so I did not include those samples in the Seurat object. The total number of samples for HC and E is two.
5. The Seurat object for each disease state and location combination is stored as Rdata format in Rdata folder.
6. In analysis folder, I uploaded analysis R script which contains the code to load Rdata.

## Data Processing

## Data Analysis
