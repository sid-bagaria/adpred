.. adpred's docs documentation master file, created by
   sphinx-quickstart on Sun May 10 10:07:39 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

####################
ADpred documentation
####################

Introduction
============

adpred is a tool to predict Transcription Activation Domains from protein sequences.  
It includes **methods** and a **class** to simplify the commands.  

In addition it includes an Web Application for those afraid of using terminal.
This application is called **adpred-app**

For single predictions use to the `adpred-web server <adpred.fredhutch.org>`_ 

**Adpred** uses psipred to predict secondary structure of the sequences, which, 
together with the sequence conforms the input of the CNN model.

By default, adpred uses psipred server to calculate the secondary structure. However, 
the user can define (after loading the module)::

   local_psipred = <path/to/psipred_run>  # e,g, '~/some-path/run_psipred'

or, if run_psipred is loaded into the environment::

   local_psipred_env = <command>  # depending on the version, In my case is 'run_psipred'

and adpred will use the local psipred installation.



The simple and quick use is::

   # import the module  
   import adpred
   
   # with own sequence:
   sequence = "ATREFERTATREFERTAADDWLNDCWATREFERTA"

   # with uniprot identifier (example gcn4)
   gcn4 = adpred.protein('GCN4_YEAST')

   # if secondary structure is not known:
   gcn4.predict()   # This will predict secondary structure

   # If you wish to first get the secondary structure
   gcn4.predict_second_struct()

   # do saturated mutagenesis analysis to 30mer between residue positions 108 and 138.
   gcn4.saturated_mutagenesis(start=108, end=138)

   # By default the WT structure is used for all mutants, however, 
   # If you wish to recalculate the secondary structure for each mutant
   gcn4.saturated_mutagenesis(start=108, end=138, 'second_struct_on_each_mutant')


Contents
========

.. toctree::
   :maxdepth: 6
   :titlesonly:
   :glob:

   Installation <installation>
   Module Methods <methods>
   Class protein <protein-class>
   Application <application>

Credits
=======

* Ariel Erijman
