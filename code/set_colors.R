## Purpose:  Set the global colors based on important factors
## Author: Marian L. Schmidt 
 

# Set global colors for the different taxonomic phyla
phylum_colors <- c( 
  Acidobacteria = "navy", 
  Actinobacteria = "darkslategray2", 
  Alphaproteobacteria = "tomato", 
  Aminicenantes = "cornflowerblue",
  Armatimonadetes = "wheat", 
  Bacteria_unclassified = "#508578", 
  Bacteroidetes = "gold", 
  Betaproteobacteria = "plum1", 
  "Candidate_division_OP3" = "slategray3",
  Chlamydiae = "#A20E42",
  Chlorobi="magenta", 
  Chloroflexi="black", 
  Cyanobacteria = "limegreen",
  "Deinococcus-Thermus" = "black",
  Deltaproteobacteria = "grey", 
  Firmicutes = "#3E9B96",
  Gammaproteobacteria = "cyan",
  Gemmatimonadetes = "yellow",
  Gracilibacteria = "#FD823F",
  JTB23 = "#B5D6AA",
  Latescibacteria = "salmon4",
  Lentisphaerae = "palevioletred1",
  Nitrospirae = "forestgreen",
  Omnitrophica = "violet",
  Microgenomates = "blue",
  Parcubacteria = "#531A4D",
  Planctomycetes = "darkorange", 
  Proteobacteria_unclassified = "greenyellow",
  Spirochaetae = "royalblue",
  TA06 = "peachpuff",
  TA18 = "burlywood", 
  TM6 = "olivedrab",
  Verrucomicrobia = "purple4",
  "WCHB1-60" = "green")

# Set global colors for the filter fraction sizes 
fraction_colors <- c(
  Particle = "#FF6600",
  Free = "skyblue",
  WholePart =  "firebrick3", 
  WholeFree = "cornflowerblue",
  Sediment = "#8A2667")

# Set global colors for the depth of the sample
depth_colors <- c(
  Top = "#1AB58A",
  Bottom = "#2F438A")

# Set global colors for the lakesite (i.e. sampling station)
lakesite_colors <- c(
  MIN = "#C24704",
    River = "#C24704",
  MBR = "#D9CC3C",
    Bear = "#D9CC3C",
  MDP = "#A0E0BA",
    Deep = "#A0E0BA",
  MOT = "#00ADA7",
    Outlet = "#00ADA7")

lakesite_shapes <- c(
  MIN = 24,
    River = 24,
  MBR = 23,
    Bear = 23,
  MDP = 22,
    Deep = 22,
  MOT = 21,
    Outlet = 21)

# Set global colors for the season
season_colors <- c(
  Spring = "#675B78",
  Summer = "#AC6C82",
  Fall = "#EB7B88")

season_shapes <- c(
  Spring = 25,
  Summer = 23,
  Fall = 21)


# Set global colors for the year
year_colors <- c(
  "2014" = "#F8BB0A", 
  "2015" = "#FF4E50")

# Set global colors for the pvalue for MPD and MNTD analyses
pd_colors <- c(
  "high_pval" = "red", 
  "insignificant" = "grey",
  "low_pval" = "blue")

# Set global colors for testing rare taxa
tons_colors <- c(
  "1-tons" = "#3B556A", 
  "5-tons" = "#427276", 
  "10-tons" = "#3F9A7A", 
  "30-tons" = "#94C660", 
  "60-tons" = "#F3E99F", 
  "90-tons" = "#BFCFBB", 
  "150-tons" = "#6F755F", 
  "225-tons" = "#403F33",
  "300-tons" = "#4D2B2F")
