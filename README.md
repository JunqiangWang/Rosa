# Rosa (version 0.777)
<img src="/Illustrative_Figures/Rosea_logo2.png" width="200" height="200">
Rosa is an algorithm that implements Regulon set enrichment analysis

version 0.777

### The illustrative framework of Rosea

Rosa—RegulOn Structure-based Activity inference
Rosa performs regulon structure-based protein activity inference, which is a statistical model for quantitative inference of the protein activities using reversed-engineered gene regulatory networks from scRNA-seq.
There are several elements in the model: first, candidate regulators including TFs, coTFs, surface proteins and signaling transduction proteins, and associated targets; second, the regulatory modes (mod_reg) including positive regulation or the negative regulation; third, the weight of regulatory strength (w_reg); fourth, the normalized gene expression Z scores; fifth, the weight of the gene expression variation (w_var, optional). 
Rosea calculates a normalized enrichment score (NES) representing the relative protein activities of the candidate regulators. The NES calculation in Rosea is equivalent to the weighted Stouffer integration of the Z scores. 


### Author 
Junqiang Wang

Email: junqiangwang333@gmail.com



