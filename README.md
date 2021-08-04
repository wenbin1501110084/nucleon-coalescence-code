# nucleon-coalescence-code
nucleon coalescence 
This the nucleons coalaescenc code to generate the light nuclei, deuteron and triton, by the coalaescenc process. To compile code, please use the command:
      gfortran -o Triton cluster_trition_output_spectra.f 
      gfortran -o Deuteron cluster_deuteron_output_spectra.f 
Then use ./Triton path_to_nucleon_file, or ./Deuteron path_to_nucleon_file to output the pT-spectra and dN/dy of triton and deuteron. nucleon_test.txt is the example file of the nucleon input file.
 
When you use these code,please cite papers:

@article{Zhao:2021dka,
    author = "Zhao, Wenbin and Sun, Kai-jia and Ko, Che Ming and Luo, Xiaofeng",
    title = "{Multiplicity Scaling of Light Nuclei Production in Relativistic Heavy-Ion Collisions}",
    eprint = "2105.14204",
    archivePrefix = "arXiv",
    primaryClass = "nucl-th",
    month = "5",
    year = "2021"
}

@article{Zhao:2020irc,
    author = "Zhao, Wenbin and Shen, Chun and Ko, Che Ming and Liu, Quansheng and Song, Huichao",
    title = "{Beam-energy dependence of the production of light nuclei in Au + Au collisions}",
    eprint = "2009.06959",
    archivePrefix = "arXiv",
    primaryClass = "nucl-th",
    doi = "10.1103/PhysRevC.102.044912",
    journal = "Phys. Rev. C",
    volume = "102",
    number = "4",
    pages = "044912",
    year = "2020"
}

@article{Zhao:2018lyf,
    author = "Zhao, Wenbin and Zhu, Lilin and Zheng, Hua and Ko, Che Ming and Song, Huichao",
    title = "{Spectra and flow of light nuclei in relativistic heavy ion collisions at energies available at the BNL Relativistic Heavy Ion Collider and at the CERN Large Hadron Collider}",
    eprint = "1807.02813",
    archivePrefix = "arXiv",
    primaryClass = "nucl-th",
    doi = "10.1103/PhysRevC.98.054905",
    journal = "Phys. Rev. C",
    volume = "98",
    number = "5",
    pages = "054905",
    year = "2018"
}


