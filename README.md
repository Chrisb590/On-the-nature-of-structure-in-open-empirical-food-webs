 <font size="+3">This repository contains all files used in <ins>On the nature of structure in collections of freely available food webs</ins>. </font>

<ul>
  <li>Code to reproduce all numerics in the main manuscript (i.e., mean pairwise GCD-11) is in <ins>brimacombe_et_al_2024.R</ins>.</li>
  <li>Code to reproduce all median numerics in the appendix manuscript is in <ins>brimacombe_et_al_2024_median.R</ins>.</li>
  <li>Code to reproduce the multidimensional plot (Figure 3) is in <ins>mds.py</ins>.</li>
  <li>All metadata for the associated networks is in <ins>metadata.csv</ins>.</li>
  <li><ins>metadata_for_mds.csv</ins> is the metadata file without bibtex citations, used in <ins>mds.py</ins> (Python doesn't like bibtex code).</li>
  <li>All computed pairwise graphlet correlated distance-11s between networks is in <ins>gcd11.csv</ins>.</li>
  <li>All networks used in the manuscript are in the folder <ins>networks</ins>.
  <ul>
    <li>Networks are categorized and organized into folders based on their respective source, which can be traced back to one of three distinct repositories: (1) https://www.canberra.edu.au/globalwebdb/, (2) http://www.ecologia.ib.usp.br/iwdb/resources.html, and (3) https://www.web-of-life.es/.
    <ul>
      <li>Within each source's folder, networks are organized into two folders: 
        <ul>
        <li><ins>raw_networks</ins> which contain the unaltered networks downloaded from the repository.</li> 
        <li><ins>edgeList_binary_white_space_removed_no_cannibalism_no_duplicate</ins> which contain the edge list of networks used in the manuscript, and have been modified such that there is no white spaces in node names, no cannibalism edges, and no duplicate nodes. Moreover, networks were required to have at least five unique predator and prey nodes, and at least 10 nodes in total.</li>
        </ul>
    </ul>
    </li>
  </ul>
  </li>
</ul>
 
