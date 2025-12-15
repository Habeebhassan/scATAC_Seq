# --- Prerequisites: Load Library and Set Genome (Task 1.1) ---
library(ArchR) 

# Set the reference genome to 'hg38'
addArchRGenome("hg38") 
addArchRThreads(threads = 6) # Set parallel threads (adjust as needed)

# --- Task 1.2: Read Data and Apply Lenient Filtering ---

# 1. Rerun the command to find your input files
# ***CRITICAL: Use the corrected path and pattern based on file inspection.***
inputFiles <- list.files(
  path = "/workspace/dataset/", 
  # Pattern is adjusted for the double-extension file names
  pattern = "fragments\\.tsv\\.gz_fragments\\.tsv\\.gz$", 
  full.names = TRUE
)

print(inputFiles)

# 2. Create the Arrow Files
# We omit 'outputDirectory' as it caused an error in your ArchR version.
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  # Adjust sample name extraction to remove the long suffix
  sampleNames = gsub("_fragments.tsv.gz_fragments.tsv.gz", "", basename(inputFiles)),
  minTSS = 4,      # Lenient filtering: filterTSS > 4
  minFrags = 1000,  # Lenient filtering: filterFrags > 500
  addTileMat = TRUE, 
  addGeneScoreMat = TRUE 
)
print("--- Arrow File Creation Complete. Files saved to default 'ArrowFiles/' directory. ---")

# --- Task 1.3: Identify Doublets ---

# 1. Manually recreate the ArrowFiles character vector by listing the generated .arrow files
# ***CORRECTION***: We search the current directory ('.') for the files, 
# since you confirmed they are not in a subfolder.
ArrowFiles <- list.files(
  path = ".", # Search the current working directory
  pattern = ".arrow$", 
  full.names = TRUE
)

# 1. Add Doublet Scores
# NOTE: The ArrowFiles object is modified in memory with the scores.
ArrowFiles <- addDoubletScores(
  input = ArrowFiles,
  k = 10
)
print("--- Doublet Scoring Complete. ---")

# --- Task 1.4: Collect All Samples into a Joint Data Structure ---

# 1. Manually recreate the ArrowFiles character vector for ArchRProject()
# This handles the object type change from addDoubletScores().
ArrowFiles <- list.files(
  path = ".", 
  pattern = ".arrow$", 
  full.names = TRUE
)

# 2. Create the ArchR Project object
print("--- Creating ArchR Project ---")
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  # This output directory is for the main project structure (Crucial for organizing results)
  outputDirectory = "scATAC_Project", 
  copyArrows = TRUE
)

# 3. Filter doublets from the project
print("--- Filtering Doublets ---")
proj <- filterDoublets(proj)

print("--- Project Creation and Doublet Filtering Complete ---")

# 4. Inspect and Report Metadata (Answering Task 1.4 Questions)
cat("\n--- Task 1.4 Answers (Post-Filtering) ---")

# How many cells are included in your project? 
n_cells <- nCells(proj)
cat("The number of cells included are: ",n_cells)

# What is the median TSS value and the median number of fragments?
median_tss <- median(proj$TSSEnrichment, na.rm = TRUE)
median_frags <- median(proj$nFrags, na.rm = TRUE)
cat(paste0("\n• Median TSS Enrichment Score: ", round(median_tss, 2)))
cat(paste0("\n• Median Number of Fragments: ", format(median_frags, big.mark = ",")))

# What is the dimension of the tile set your dataset?
tile_ranges <- getFeatures(proj, useMatrix = "TileMatrix")
n_tiles <- length(tile_ranges)
cat(paste0("\n• Dimension of the Tile Matrix:     ", format(n_tiles, big.mark = ","), " x ", format(n_cells, big.mark = ",")))
cat("\n-------------------------------------------")

# Save the Project object
saveArchRProject(ArchRProj = proj, outputDirectory = "scATAC_Project", load = TRUE)