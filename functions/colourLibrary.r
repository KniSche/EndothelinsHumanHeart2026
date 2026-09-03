
# trifurcation colours
  cambridge.sunset <- c("#57D8F7","#3488cf","#845F96", "#D148A8","#FEA1BA","#FEFDFD", "#FCEDAA","#FFD169","#FFA762","#E15C53", "#692F34")
  
  require(RColorBrewer)
  
  cambridge.sunrise = 
  colorRampPalette(
    c(
      colorRampPalette(c("darkblue", cambridge.sunset[c(2,1,6)]))(6)[1:5],
      colorRampPalette(c(cambridge.sunset[c(6,7,8,9)],"red"))(6)[2:6]
    )
 )(11)
 
# Some color scheme (similar to / inspired by Viridis) 
    alt.cols <- c("#332366", "#44599E", "#4AB38D", "#78CB4E", "#FFE458")

# plot schemes
  plot.dark <- function(){
    par(bg=rgb(0.1,0.12,0.12), col="white", col.lab="white", col.axis="white", col.main="white", col.sub="white", fg="white")
  }
  plot.light <- function(){
    par(bg="white", col="black", col.lab="black", col.axis="black", col.main="black", col.sub="black", fg="black")
  }

# aja.yellow
   aja.blue <- "#614fff"
   aja.yellow <- "#ffce1f"
   
# my custom green library:
  verdant = c(
    "#020f01",
    "#192819",
    "#365126",
    "#4f7836",
    "#6d9b40",
    "#8dbe40",
    "#aed135",
    "#c4da4a",
    "#d8e476",
    "#eceea2",
    "#f6f5c5",
    "#fbfcec"
  )



# Heat-like scale   
   Heat <- c(
   "#8e133b",
   "#b9534c",
   "#da8459",
   "#eeab65",
   "#f6c971",
   "#f1de81",
   "#e2e6bd"   
   )[7:1]
   Heat.pal=colorRampPalette(col=Heat)
   
# bifurcation colours
   Earth <- c(
    "#8c510a",
    "#bf812d",
    "#dfc27d",
    "#f6e8c3",
    "#f5f5f5",
    "#c7eae5",
    "#80cdc1",
    "#35978f",
    "#01665e"
    )

# discrete data...     
  Discrete <- c(
    "#000000",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E422",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
    "#999999"
    )
 
 
# Discrete but up to 12
  Discrete12 = c(
  '#a6cee3',
  '#1f78b4',
  '#b2df8a',
  '#33a02c',
  '#fb9a99',
  '#e31a1c',
  '#fdbf6f',
  '#ff7f00',
  '#cab2d6',
  '#6a3d9a',
  '#ffff99',
  '#b15928'
  )
  
  
# Discrete Ultimate edition
  Discrete.U = c(
  '#000000', # 1        
  '#606060', # 2
  '#DFDFDF', # 3

  '#780033', # 4   
  '#C00052', # 5 
  '#D86093', # 6 
 
  '#006978', # 7 
  '#00A8C0', # 8 
  '#60C9D8', # 9 

  '#833515', # 10   
  '#D25521', # 11
  '#E39574', # 12 

  '#600000', # 13
  '#C00000', # 14
  '#D86060', # 15
  
  '#786E05', # 16
  '#F1DD0B', # 17
  '#F8EE85', # 18
  
  '#003E60', # 19
  '#007CC0', # 20
  '#80BEE0', # 21
  
  '#007800', # 22
  '#00C000', # 23
  '#80E080', # 24
  
  '#380060', # 25
  '#7000C0', # 26
  '#B780E0' # 27
  )

 Discrete.U2 = Discrete.U[c(seq(2,27,3), seq(3,27,3), seq(1,27,3))]
  
  
# Project.J
  nc.yellow = "#FDD466"
  cm.red    = "#E06866"
  epi.blue  = "#6D9EEB"
  nccm.orange = "#EF9E66"
  epicm.purple = "#A582A9"
  
  
  
