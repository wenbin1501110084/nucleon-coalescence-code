This the nucleons coalaescenc code to generate the light nuclei, deuteron and triton, by the coalaescenc process. To compile code, please use the command:

      gfortran -o Triton cluster_trition_output_spectra.f 
      gfortran -o Deuteron cluster_deuteron_output_spectra.f 

Then use 
      ./Triton path_to_input_nucleon_file 
and
      ./Deuteron path_to_input_nucleon_file 
to output the pT-spectra and dN/dy of triton and deuteron. The format for the input nucleon is:

      number_of_nucleon
      id, px, py, pz, xm, x, y, z, t
      
If you have any questions about it, please let me know. My email is: 1501110084@pku.edu.cn
Wenbin Zhao
 



