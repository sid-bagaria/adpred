########################
Class protein Attributes
########################

Attributes
==========
* prot_id                                
* sequence                                           
* second_struct                                           
* predictions                                            
* meta_data (usually header of fasta in uniprot) 

#####################
Class protein Methods
#####################

predict
=======
.. code-block:: python
    :linenos:

    predict()                                                 
    
what it does                                                                
------------                                                                
Assigns to each aminoacid of the sequence the probability of being in a     
AD region and populates the object attribute.                                                       

predict_second_struct
=====================
.. code-block:: python
    :linenos:

    predict_second_struct()                                            
    
what it does                                                                
------------                                                                
predict secondary structure of protein object and populates object attribute.

Saturated mutagenesis
=====================

.. code-block:: python
    :linenos:

    saturated_mutagenesis(sequence, second_struct, predictions, *args)

