-- 1. Create the samples table for your first CSV file metadata
CREATE TABLE IF NOT EXISTS samples (
    sample_id VARCHAR(50) PRIMARY KEY,
    dex VARCHAR(20) NOT NULL,          -- e.g., 'treated' or 'control'
    celltype VARCHAR(50) NOT NULL,     -- e.g., 'N61311', 'N052611', etc.
    geo_id VARCHAR(50) UNIQUE NOT NULL -- e.g., GEO series sample accession
);

-- 2. Create the gene expression table for your second CSV file
CREATE TABLE IF NOT EXISTS gene_expression (
    id SERIAL PRIMARY KEY,
    sample_id VARCHAR(50) NOT NULL,
    ensg_id VARCHAR(50) NOT NULL,      -- Ensembl Gene ID (e.g., ENSG00000...)
    expression_value NUMERIC(12, 4),   -- Counts, TPM, or FPKM expression value
    -- Ensures expression rows link strictly to a valid sample in the samples table
    CONSTRAINT fk_sample 
        FOREIGN KEY (sample_id) 
        REFERENCES samples (sample_id) 
        ON DELETE CASCADE

    CONSTRAINT uq_gene_expression_sample_gene UNIQUE (sample_id, ensg_id)
);

-- 3. Create an index to speed up lookups by gene or sample ID during analysis
CREATE INDEX IF NOT EXISTS idx_gene_sample ON gene_expression(ensg_id, sample_id);