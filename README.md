This the nucleons coalescence code to generate the light nuclei, deuteron and triton, 
by the coalescence process. 

To compile code, please use the command:

      gfortran -o Triton cluster_trition_output_spectra.f 
      gfortran -o Deuteron cluster_deuteron_output_spectra.f 

Then use 

      ./Triton path_to_input_nucleon_file 

and

      ./Deuteron path_to_input_nucleon_file 

to output the pT-spectra and dN/dy of triton and deuteron. 
The format for the input nucleon is:

      number_of_nucleon
      id, px, py, pz, mass, x, y, z, t

The input.txt is the parameter file. You can also get the (anti)hypertriton, 
(anti)helium3 in this nucleon coalaescence code.

If you have any questions about it, please send message to email zhaowenb@pku.edu.cn

Wenbin Zhao
 



