This repository contains all files used in <ins>**On the nature of structure in open empirical food webs**</ins>.

Code to reproduce all numerics in the manuscript is in **Brimacombe_et_al_2023.r**.

Code to reproduce the multidimensional plot (Figure 3) is in **mds.py**.

All metadata for the associated networks is in **metadata.csv**.

All computed pairwise graphlet correlated distance-11s between networks is in **gcd11.csv**.

<ul>
  <li>All networks used in the manuscript are in the folder <ins>networks</ins>.
  <ul>
  <li>Networks are categorized and organized into folders based on their respective source, which can be traced back to one of three distinct repositories: (1) https://www.canberra.edu.au/globalwebdb/, (2) http://www.ecologia.ib.usp.br/iwdb/resources.html, and (3) https://www.web-of-life.es/.
  <ul>
  <li>Within each source's folder, networks are organized into two folders: (1) raw_networks which contain the unaltered networks downloaded from the repository, and (2) edgeList_binary_white_space_removed_no_cannibalism_no_duplicate: which contain the edge list of networks used in the manuscript, and have been modified such that there is no white space in species names, no cannibalism, and no duplicate nodes.</li>
  </ul>
  </li>
  </ul>
  </li>
</ul>
 
