#######
Methods
#######

predict
=======

::
   predict(seq, struct=None):                                                  
    
what it does                                                                
------------                                                                
Assigns to each aminoacid of the sequence the probability of being in a     
AD region.                                                                  
                                                                            
parameters                                                                  
----------                                                                  
  - sequence: sequence of amino acids (aa). OR uniprot ID!                  
    ACHTUNG!! ==> Use a valid id or <name>_<specie> or else you might       
    get the wrong result.                                                   
  - second_struct: sequence of secondary elements of each aa in sequence.   
                                                                            
returns                                                                     
-------                                                                     
  - numpy array with probabilities of AD for each aa in the input sequence.


Saturated mutagenesis
=====================
::
   saturated_mutagenesis(sequence, second_struct, predictions, *args)         

what is does                                                                
------------                                                                
Uses mut_analysis to calculate saturated mutagenesis                        
                                                                            
parameters                                                                  
----------                                                                  
    - sequence: Protein sequence                                            
    - second_struct: string of protein secondary structure elements         
    - predictions: numpy array of adpred predictions                        
                                                                            
    - optionals                                                             
        - 'second_struct_on_each_mutant': string, second_struct will be computed
          on each mutant.                                                   
                                                                            
returns                                                                     
-------                                                                     
    - 2d array of positions and adpred probabilities for each of the 30 residues
      in the sequence. The order follows the order of the list utils.aa 



Calculate secondary structure
=============================
::
   calculate_psipred(fasta_name)
                           
what it does                                                                
------------                                                                
predict secondaary structure                                                
                                                                            
parameters                                                                  
----------                                                                  
    - fasta_name: filename of fasta file                                    
                                                                            
results                                                                     
-------                                                                     
    -string of secondary structure elements.    

