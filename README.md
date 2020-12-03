# End-of-epidemic-declarations

Matlab and R code for computing the probability of elimination, z_s, of an infectious disease epidemic at time s. Adopts an event-triggered
approach providing end-of-epidemic declaration times with 95% confidence when z_s > 0.95.

Approach examines several reproduction number profiles and computes z_s under renewal models featuring serial interval distributions
from several infectious diseases. Current implementation assumes perfect surveillance and compares biases resulting from under-reporting
or case importation. Compensating for these biases will improve judgments of the end of the epidemic.


Code author: Kris V Parag

Cite as: An exact method for quantifying the reliability of end-of-epidemic declarations in real time
         by Kris V Parag, Christl A. Donnelly, Rahul Jha & Robin N Thompson at
         PLoS Comput Biol: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008478
           

