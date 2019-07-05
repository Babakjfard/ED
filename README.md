# Ecological Dynamics
Studying food webs and their changes and effects is very important. It can lead into solutions for maintain the food resources for the human beings, saving species from extinction, planning for optimum harvesting, and many other applications which makes it very important. The important problem of climate change and its effect of food chains and species make it even more important to research. Through the times different approaches and methods have been used to solve the complexity of food webs. Introduction of the graph theory, or network science, and the advent of metacommunity models were good advancements in modeling the behaviors of food webs; leading to an important conclusion that complexity maintains stability, in opposite of what was suggested earlier by some scientists. This research is trying to model the interaction between food webs as models of networks of networks
# Methods and Materials
In the proposed model we have used combinations of different methods to make a more realistic model of connected food webs. In the local level, the predator-prey connections in each food web are modeled by an adjacency matrix of that web.

prey’s growth rate:      ![This equation](https://latex.codecogs.com/gif.latex?dN/dt%3DrN%281-N/K%29-aNP/%281&plus;ahN%29) . (1) .  

growth rate of predator:      ![Second equation](https://latex.codecogs.com/gif.latex?dP/dt%3DabNP/%281&plus;ahN%29-mP) . (2) . 

The system also includes basal species which compete for one resource. The competition is modeled by mechanistic model as in equations (3) and (4). 

![3rd eq](https://latex.codecogs.com/gif.latex?%28dN_i%29/dt%3DN_i%5C%20%28%5Cmu_i%5C%20%5C%20R/%28K_i&plus;R%29-M_i%5C%20%29) . (3)

<img src="https://github.com/Babakjfard/ED/blob/master/eq_04.png" width="300"> (4)

The metacommunity model is then used to model the regional level of interactions, which are represented by dispersals between the neighboring food webs as patches for dispersal model. Figure 1 depicts a simple model with two food webs and their resources, which are governed by dispersal rules as follow:


