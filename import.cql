// Create indexes for faster lookup
CREATE INDEX ON :ORGANISM(sciName);
CREATE INDEX ON :GENE(symbol);

// Create constraints
CREATE CONSTRAINT ON (g:GENE) ASSERT g.entrezId IS UNIQUE;
CREATE CONSTRAINT ON (o:ORGANISM) ASSERT o.taxId IS UNIQUE;
CREATE CONSTRAINT ON (p:PLATFORM) ASSERT p.platAccess IS UNIQUE;
CREATE CONSTRAINT ON (e:EXPERIMENT) ASSERT e.nameExp IS UNIQUE;
CREATE CONSTRAINT ON (c:CLUSTER) ASSERT c.clustInfo IS UNIQUE;
CREATE CONSTRAINT ON (b:BP) ASSERT b.GOTerm IS UNIQUE;

// Create Genes 
USING PERIODIC COMMIT
LOAD CSV WITH HEADERS FROM "file:///data.csv" as row
CREATE (GENE:GENE {entrezId: row.GeneID, synonyms: row.Aliases, symbol: row.Symbol, description: row.description, mapLocation: row.map_location, chromosome: row.chromosome,
currentId: row.CurrentID, status: row.Status, otherDesignations: row.other_designations,genomicNuclVers: row.genomic_nucleotide_accession_version,
Start:row.start_position_on_the_genomic_accession,End:row.end_position_on_the_genomic_accession,orientation:row.orientation,exon_count:row.exon_count,OMIM: row.OMIM});

// Create Organism
USING PERIODIC COMMIT
LOAD CSV WITH HEADERS FROM "file:///taxon.csv" as row
CREATE (ORGANISM:ORGANISM {taxId: row.tax_id, sciName: row.Org_name});

// Create relationships: Genes to Organisms
USING PERIODIC COMMIT
LOAD CSV WITH HEADERS FROM "file:///data.csv" AS row
MATCH (GENE:GENE {entrezId: row.GeneID})
MATCH (ORGANISM:ORGANISM {taxId: row.tax_id})
MERGE (ORGANISM)-[:Has_a]->(GENE); 

//Create relationship: PPI
USING PERIODIC COMMIT
LOAD CSV WITH HEADERS FROM "file:///PPI.csv" AS row
MATCH (GENE1:GENE {entrezId: row.protein1})
MATCH (GENE2:GENE {entrezId: row.protein2})
MERGE (GENE1)-[rel:PPI_interaction]->(GENE2)
ON CREATE SET rel.neighborhood = toInt(row.neighborhood), rel.fusion = toInt(row.fusion),
rel.cooccurence = toInt(row.cooccurence), rel.coexpression = toInt(row.coexpression), rel.experimental = toInt(row.experimental),
rel.database = toInt(row.database),rel.textmining = toInt(row.textmining),rel.combined_score = toInt(row.combined_score); 
