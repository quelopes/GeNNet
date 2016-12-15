// Create indexes for faster lookup
//CREATE INDEX ON :Category(categoryName);
//CREATE INDEX ON :SubCategory(subCategoryName);
CREATE INDEX ON :ORGANISM(sciName);
CREATE INDEX ON :GENE(symbol);

// Create constraints
//CREATE CONSTRAINT ON (t:TF) ASSERT t.geneType IS UNIQUE;
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

// Create Transcription Factors
//USING PERIODIC COMMIT
//LOAD CSV WITH HEADERS FROM "file:///TFs.csv" AS row
//CREATE (TF:TF {geneType: "Transcription Factor"});

// Create relationships: TFs to Genes 
//USING PERIODIC COMMIT
//LOAD CSV WITH HEADERS FROM "file:///TFs.csv" AS row
//MATCH (GENE:GENE {entrezId: row.entrez_human})
//MATCH (TF:TF {geneType: "Transcription Factor"})
//MERGE (GENE)-[it:Is_a]->(TYPE) 
//ON CREATE SET it.AnimalTFDB = toInt(row.AnimalTFDB), it.JASPAR = toInt(row.JASPAR),
//it.GOC = toInt(row.GOC), it.DbTF = toInt(row.DbTF), it.Ravasi = toInt(row.Ravasi_et_al),
//it.TcoF = toInt(row.TcoF_DB),it.TFCat = toInt(row.TFCat),it.DBD = toInt(row.DBD),
//it.Vaquerizas = toInt(row.Vaquerizas_et_al),it.ORFeome = toInt(row.ORFeome);
// [1] "gene_symbol"      "entrez_human"     "entrez_rat"       "entrez_mouse"     "Unirprot_Human"  
// [6] "Uniprot_Rat"      "Uniprot_mouse"    "gene_name"        "synonym"          "DbTF"            
//[11] "pubmed_id"        "evidence.type"    "AnimalTFDB"       "DBD"              "GOC"             
//[16] "JASPAR"           "ORFeome"          "Ravasi.et.al"     "TcoF.DB"          "TFCat"           
//[21] "TFClass_human"    "TFclass_mouse"    "Vaquerizas.et.al" "pubmed_id.1"      "taxon"           
//[26] "pubmed_id.2"      "taxon.1"          "pubmed_id.3"      "taxon.2"          "pubmed_id.4"     
//[31] "taxon.3"          "pubmed_id.5"      "taxon.4"          "pubmed_id.6"      "taxon.5"         
//[36] "GOC981.select_EC"


// Create relationships: Genes to Genes Orthologous
//USING PERIODIC COMMIT
//LOAD CSV WITH HEADERS FROM "file:///home/quelopes/Dropbox/GeNNet-Background-v.2.0/Background/Orthologous/OrthologousHsMmu.csv" AS row
//MATCH (GENE1:GENE {entrezId: row.gene1})
//MATCH (GENE2:GENE {entrezId: row.gene2})
//MERGE (GENE1)-[rela:Are_Orthologos]->(GENE2)
//ON CREATE SET rela.Source="IParanoid";

// Create Platform
//USING PERIODIC COMMIT
//LOAD CSV WITH HEADERS FROM "file:///home/quelopes/Dropbox/GeNNet-Background-v.2.0/Background/Platform/Platforms.csv" as row
//CREATE (PLATFORM:PLATFORM {platAccess: row.PlatAccession, platName: row.PlatformName,Organism:row.Organism,id:row.ID,DataSets:row.DataSets,series:row.Series,samples:row.Samples});

// Create relationships: Platform to Organism
//USING PERIODIC COMMIT
//LOAD CSV WITH HEADERS FROM "file:///home/quelopes/Dropbox/GeNNet-Background-v.2.0/Background/Platform/Platforms.csv" AS row
//MATCH (ORGANISM:ORGANISM {taxId: row.Org_Id})
//MATCH (PLATFORM:PLATFORM {platAccess: row.PlatAccession})
//MERGE (PLATFORM)-[v:Was_based_on]->(ORGANISM);

//HIV-1_GeneID,HIV-1_Prot_Acc,HIV-1_Prot_Name,Keyword,Human_GeneId,Human_GeneSymbol,Human_Prot_Acc,Human_Prot_Name,PMID(s),Modified_Date,Interaction_Desc

