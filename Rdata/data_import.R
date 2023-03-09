#############################
### BIOS 785 Team Project ###
######## Data Import ########
#############################


# Set WD
setwd("D:/UNC-PhD/Spring 2023/BIOS 785/teamproject/data")


# Load library
library(Seurat)


# data import
## Seurat Object creating Function
somaker <- function(x, state, location){
  barcode <- x
  feature <- x + 1
  matrix <- x + 2
  
  dt <- ReadMtx(mtx = list.files(pattern = paste0(state, "_\\d+.", location))[matrix],
                cells = list.files(pattern = paste0(state, "_\\d+.", location))[barcode],
                features = list.files(pattern = paste0(state, "_\\d+.", location))[feature])
  
  so_dt <- CreateSeuratObject(counts=dt)
  
  return(so_dt)
}



####################
## DLE and Epidermal
list.files(pattern="DLE_\\d+.E_")
DLE_1E <- somaker(1, "DLE", "E")
DLE_5E <- somaker(4, "DLE", "E")
DLE_6E <- somaker(7, "DLE", "E")
DLE_7E <- somaker(10, "DLE", "E")
DLE_8E <- somaker(13, "DLE", "E")

### Merge DLE_E
DLE_E <- merge(DLE_1E, y = c(DLE_5E, DLE_6E, DLE_7E, DLE_8E), add.cell.ids = c("DLE_1E", "DLE_5E", "DLE_6E", "DLE_7E", "DLE_8E"))

### Check
GetAssayData(DLE_E, slot="counts")[1:20, 1:50]
rownames(GetAssayData(DLE_E, slot="counts"))

### Save as Rdata
save(DLE_E, file="DLE_E")



#################
## DLE and Dermal
list.files(pattern="DLE_\\d+.D_")
DLE_1D <- somaker(1, "DLE", "D")
DLE_2D <- somaker(4, "DLE", "D")
DLE_3D <- somaker(7, "DLE", "D")
DLE_4D <- somaker(10, "DLE", "D")
DLE_5D <- somaker(13, "DLE", "D")

### Merge DLE_E
DLE_D <- merge(DLE_1D, y = c(DLE_2D, DLE_3D, DLE_4D, DLE_5D), add.cell.ids = c("DLE_1D", "DLE_2D", "DLE_3D", "DLE_4D", "DLE_5D"))

### Check
GetAssayData(DLE_D, slot="counts")[1:20, 1:50]
rownames(GetAssayData(DLE_D, slot="counts"))

### Save as Rdata
save(DLE_D, file="DLE_D")



####################
## SLE and Epidermal
list.files(pattern="SLE_\\d+.E_")
SLE_6E <- somaker(1, "SLE", "E")
SLE_7E <- somaker(4, "SLE", "E")
SLE_8E <- somaker(7, "SLE", "E")
SLE_9E <- somaker(10, "SLE", "E")
SLE_10E <- somaker(13, "SLE", "E")

### Merge DLE_E
SLE_E <- merge(SLE_6E, y = c(SLE_7E, SLE_8E, SLE_9E, SLE_10E), add.cell.ids = c("SLE_6E", "SLE_7E", "SLE_8E", "SLE_9E", "SLE_10E"))

### Check
GetAssayData(SLE_E, slot="counts")[1:20, 1:50]
rownames(GetAssayData(SLE_E, slot="counts"))

### Save as Rdata
save(SLE_E, file="SLE_E")



#################
## SLE and Dermal
list.files(pattern="SLE_\\d+.D_")
SLE_1D <- somaker(1, "SLE", "D")
SLE_2D <- somaker(4, "SLE", "D")
SLE_3D <- somaker(7, "SLE", "D")
SLE_4D <- somaker(10, "SLE", "D")
SLE_5D <- somaker(13, "SLE", "D")
SLE_6D <- somaker(16, "SLE", "D")
SLE_7D <- somaker(19, "SLE", "D")

### Merge DLE_E
SLE_D <- merge(SLE_1D, y = c(SLE_2D, SLE_3D, SLE_4D, SLE_5D, SLE_6D, SLE_7D), add.cell.ids = c("SLE_1D", "SLE_2D", "SLE_3D", "SLE_4D", "SLE_5D", "SLE_6D", "SLE_7D"))

### Check
GetAssayData(SLE_D, slot="counts")[1:20, 1:50]
rownames(GetAssayData(SLE_D, slot="counts"))

### Save as Rdata
save(SLE_D, file="SLE_D")



########################
## Healthy and Epidermal
list.files(pattern="HC_\\d+.E_")
#HC_2E <- somaker(1, "HC", "E")
HC_3E <- somaker(4, "HC", "E")
HC_4E <- somaker(7, "HC", "E")
#HC_5E <- somaker(10, "HC", "E")

### Merge DLE_E
HC_E <- merge(HC_3E, y = c(HC_4E), add.cell.ids = c("HC_3E", "HC_4E"))

### Check
GetAssayData(HC_E, slot="counts")[1:20, 1:50]
rownames(GetAssayData(HC_E, slot="counts"))

### Save as Rdata
save(HC_E, file="HC_E")



#####################
## Healthy and Dermal
list.files(pattern="HC_\\d+.D_")
HC_1D <- somaker(1, "HC", "D")
HC_2D <- somaker(4, "HC", "D")
HC_4D <- somaker(7, "HC", "D")
HC_3D <- somaker(10, "HC", "D")

### Merge DLE_E
HC_D <- merge(HC_1D, y = c(HC_2D, HC_3D, HC_4D), add.cell.ids = c("HC_1D", "HC_2D", "HC_3D", "HC_4D"))

### Check
GetAssayData(HC_D, slot="counts")[1:20, 1:50]
rownames(GetAssayData(HC_D, slot="counts"))

### Save as Rdata
save(HC_D, file="HC_D")
