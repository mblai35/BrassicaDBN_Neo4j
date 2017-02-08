# Neo4j_visNetwork.R
# R version 3.3.1 (2016-06-21)
# February 5, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Create interactive visualization for  pheno network for 
# Brassica data using bnlearn package. Data taken from Brassica control
# and droughted conditions. 


# ERRORS!!!!!!!!


#-----------------------------------------------------------------------
library(igraph)
library(visNetwork)
library(RNeo4j)
#-----------------------------------------------------------------------

# Start Neo4j graph. 
graph = startGraph("http://localhost:7474/db/data/", username = "neo4j",
                   password = "plantanalytics")

# Query for Intervention effects nodes. 
query = "
MATCH (n)-[]->(m)
WHERE n.name < {x}
RETURN n.name AS from, m.name AS to
"
  
    # Return query and store in 'edges'. 
    edges <- cypher(graph, query, x = "Stem")
    
    # Create unique nodes. 
    nodes = data.frame(id=unique(c(edges$from, edges$to)))
    nodes$label = nodes$id
    
    # Create an igraph object.
    ig = graph_from_data_frame(edges, directed=F)
    
    # Make the sizes of the nodes a function of their 
    # betweenness centrality.
    nodes$value = betweenness(ig)
    
    color <- nodes$value
    color[color > 16] <- 50
    color[color == 0 ] <- 16
    color[color < 16] <- 5
    
    # Color by edge_betweenness. 
    # clusters = cluster_edge_betweenness(ig)
    nodes$group = color

    # Plot network. 
    visNetwork(nodes, edges, 
               main = "Connectedness of Nodes Across Intervention")
    
    
# Query for predawn and midday nodes. 
query = "
MATCH (n:Midday)-[:INFLUENCES]->(m) 
RETURN n.name AS from, m.name AS to
"

    # Return query and store in 'edges'. 
    edges2 <- cypher(graph, query)
    
    # Create unique nodes. 
    nodes = data.frame(id=unique(c(edges2$from, edges2$to)))
    nodes$label = nodes$id
    
    # Create an igraph object.
    ig = graph_from_data_frame(edges2, directed=F)
    
    # Make the sizes of the nodes a function of their 
    # betweenness centrality.
    nodes$value = betweenness(ig)
    
    color <- nodes$value
    color[color > 16] <- 50
    color[color == 0 ] <- 0
    color[color < 5 & color > 0] <- 2
    
    
    nodes$group = color
    # Plot network. 
    visNetwork(nodes, edges2)
    
    
# Query for Photosynthesis nodes. 
query = "
    MATCH (n:Predawn)-[:INFLUENCES]->(m:Predawn) 
    RETURN n.name AS from, m.name AS to
    "
    
    # Return query and store in 'edges'. 
    edges3 <- cypher(graph, query)
    
    # Create unique nodes. 
    nodes = data.frame(id=unique(c(edges3$from, edges3$to)))
    nodes$label = nodes$id
    
    # Create an igraph object.
    ig = graph_from_data_frame(edges3, directed=T)
    
    # Make the sizes of the nodes a function of their 
    # betweenness centrality.
    nodes$value = betweenness(ig)
    
    color <- nodes$value
    
    nodes$group = color
    # Plot network. 
    visNetwork(nodes, edges3)


