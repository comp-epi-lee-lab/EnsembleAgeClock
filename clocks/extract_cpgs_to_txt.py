import pandas as pd 
from han2020_coefs import cpg_sites

#df = pd.read_csv("clocks/YingCausAge.csv")
df = pd.DataFrame(cpg_sites)
print(df)

#cpgs = df.iloc[:, 0]
#print(cpgs)

# Specify the output file path
output_file = '/Users/yashagarwal/Documents/GitHub/EnsembleMeAgingClock/tools/han2020_cpgs_list.txt'



# Save the selected values to the text file
#df.to_csv(output_file, header=False, index=False)

with open('tools/horvath_cpgs_list.txt', 'r') as file:
    lines = file.readlines()

lines = [line.strip() for line in lines]

print(lines)
print(len(lines))