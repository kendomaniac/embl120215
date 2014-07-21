###  define class
setClass("fSet",
     representation(
		 fusionInfo="list", 
		 fusionLoc="GRangesList", 
		 fusionRNA="DNAStringSet", 
		 fusionGA="GAlignments"
		 )
		)
#constructor		
newfSet <- function(fusionInfo = list(fusionTool = NA,
                                      UniqueCuttingPositionCount = 0,
                                      SeedCount = 0, 
                                      RescuedCount = 0, 
                                      SplicePattern = NA, 
                                      FusionGene = NA,
                                      frameShift = NA), 
                             fusionLoc = GRangesList(),
		                     fusionRNA = DNAStringSet(),
		                     fusionGA = GAlignments()){
					                new(Class = "fSet", fusionInfo = fusionInfo, fusionLoc = fusionLoc, fusionRNA = fusionRNA, fusionGA = fusionGA)
		                    }		
