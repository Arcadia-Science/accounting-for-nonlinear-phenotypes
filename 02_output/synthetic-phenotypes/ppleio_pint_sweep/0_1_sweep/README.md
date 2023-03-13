genotype and phenotype output from 2d sweep through p_pleio and p_int using the make_genotype <br /> 
function in tools_for_phen_gen_creatino.py <br />
<br />
Fixed parameters: <br />
n_important_loci=100 <br />
n_loci=300 <br />
n_phen=100 <br />
n_animals=500 <br />

p_pleio and p_int had values of: [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] <br />
The organization of the data is: <br />
data[p_pleio][p_int][phens][individuals] <br />
Such that data[2][3][40][200] is the: <br />
Value of the 40th phenotype for the 200th individual where p_pleio was 0.2 and p_int was 0.3 <br />
The preceeding statement assumes python indexing where the first item in a list is position 0<br />
