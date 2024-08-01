
library(polyG)
library(here)

## run the poly-G data processing pipeline to generate the angular distance data
polyG(input_dir=here('original_data/polyG/science'), 
      results_dir=here('processed_data/polyG/science'), seed=42)

polyG(input_dir=here('original_data/polyG/natgen'), 
      results_dir=here('processed_data/polyG/natgen'), seed=42)

polyG(input_dir=here('original_data/polyG/lung'), 
      results_dir=here('processed_data/polyG/lung'), seed=42)

polyG(input_dir=here('original_data/polyG/peritoneal'), 
      results_dir=here('processed_data/polyG/peritoneal'), seed=42)




