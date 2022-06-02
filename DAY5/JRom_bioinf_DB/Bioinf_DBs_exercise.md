Bioinformatics services - exercises
================
Julia Romanowska
2022-06-02

-   [Read the data](#read-the-data)
-   [Write out results](#write-out-results)
-   [Web services](#web-services)
    -   [ensembl BioMart - genes](#ensembl-biomart---genes)
        -   [biomart on the web](#biomart_web)
        -   [using `biomaRt` package](#biomart_pack)
        -   [using results from biomart](#using-results-from-biomart)
    -   [STRING - database of protein-protein interactions
        (PPI)](#string---database-of-protein-protein-interactions-ppi)
        -   [ensembl BioMart - regulatory
            regions](#ensembl-biomart---regulatory-regions)
    -   [biomart on the web](#biomart_regul_web)
        -   [using `biomaRt` package](#biomart_regul_pack)
        -   [using results from biomart](#using-results-from-biomart-1)
    -   [GenEnhancer by GeneCards](#genenhancer-by-genecards)

# Read the data

This was done for you, since the entire dataset is too big.

``` r
diff_methyl <- read_csv(
  file.path(data_dir, "diffMethTable_site_cmp13.csv")
)
diff_methyl %>% glimpse()

qvals <- qvalue(diff_methyl$diffmeth.p.val)

diff_methyl <- diff_methyl %>%
  add_column(qval = qvals$qvalues)

# setting here to very low just to get not too many results
signif_qval <- 0.00001
diff_methyl_signif <- diff_methyl %>%
  filter(qval < signif_qval)
```

You have access to the file that contains only the most significant
results

``` r
diff_methyl_signif <- read_csv(
  file.path(data_dir, "diffMethTable.csv")
)
```

    ## Rows: 4492 Columns: 29
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (3): cgid, Chromosome, Strand
    ## dbl (26): id, Start, mean.Primary solid Tumor, mean.Recurrent Solid Tumor, m...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
diff_methyl_signif %>% glimpse()
```

    ## Rows: 4,492
    ## Columns: 29
    ## $ id                                           <dbl> 24, 199, 408, 461, 507, 1…
    ## $ cgid                                         <chr> "cg18761878", "cg17505193…
    ## $ Chromosome                                   <chr> "chr1", "chr1", "chr1", "…
    ## $ Start                                        <dbl> 568475, 896310, 995013, 1…
    ## $ Strand                                       <chr> "-", "-", "-", "+", "+", …
    ## $ `mean.Primary solid Tumor`                   <dbl> 0.52846308, 0.03686572, 0…
    ## $ `mean.Recurrent Solid Tumor`                 <dbl> 0.63372874, 0.04625797, 0…
    ## $ mean.diff                                    <dbl> -0.105265660, -0.00939224…
    ## $ mean.quot.log2                               <dbl> -0.25760545, -0.26352420,…
    ## $ diffmeth.p.val                               <dbl> 3.726879e-05, 1.054755e-0…
    ## $ `max.Primary solid Tumor`                    <dbl> 0.70135511, 0.10776697, 0…
    ## $ `min.Primary solid Tumor`                    <dbl> 0.22790224, 0.02145792, 0…
    ## $ `sd.Primary solid Tumor`                     <dbl> 0.087267475, 0.009237890,…
    ## $ `max.Recurrent Solid Tumor`                  <dbl> 0.67628324, 0.05845909, 0…
    ## $ `min.Recurrent Solid Tumor`                  <dbl> 0.52956157, 0.03632786, 0…
    ## $ `sd.Recurrent Solid Tumor`                   <dbl> 0.037504125, 0.006828548,…
    ## $ min.diff                                     <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ diffmeth.p.adj.fdr                           <dbl> 6.034336e-03, 1.204082e-0…
    ## $ combinedRank                                 <dbl> 56149, 285598, 244933, 28…
    ## $ `num.na.Primary solid Tumor`                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ `num.na.Recurrent Solid Tumor`               <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ `mean.covg.Primary solid Tumor`              <dbl> 11.108108, 12.432432, 11.…
    ## $ `mean.covg.Recurrent Solid Tumor`            <dbl> 13.000000, 13.384615, 12.…
    ## $ `min.covg.Primary solid Tumor`               <dbl> 4, 5, 3, 1, 4, 5, 4, 2, 5…
    ## $ `min.covg.Recurrent Solid Tumor`             <dbl> 6, 7, 8, 5, 6, 7, 8, 4, 7…
    ## $ `max.covg.Primary solid Tumor`               <dbl> 22, 21, 20, 22, 20, 21, 2…
    ## $ `max.covg.Recurrent Solid Tumor`             <dbl> 18, 18, 17, 20, 17, 18, 2…
    ## $ `nsamples.covg.thres5.Primary solid Tumor`   <dbl> 109, 111, 110, 109, 110, …
    ## $ `nsamples.covg.thres5.Recurrent Solid Tumor` <dbl> 13, 13, 13, 13, 13, 13, 1…

Now, `diff_methyl_signif` contains only the selected data.

# Write out results

Some services accept user data, but they require it to be in specific
format. It’s easy to write various formats with R.

``` r
# write out BED file
write_delim(
  diff_methyl_signif %>%
    mutate(End = Start + 1) %>%
    dplyr::select(Chromosome, Start, End, Name = cgid, Score = diffmeth.p.val, Strand),
  file = file.path(data_dir, "signif_CpGs.BED"),
  delim = " ",
  col_names = FALSE
)

# write out region data
write_delim(
  diff_methyl_signif %>%
    mutate(strand_number = ifelse(Strand == "-", -1, 1)) %>%
    mutate(Chr = substr(Chromosome, 4, nchar(Chromosome))) %>%
    mutate(reg = paste0(Chr, ":", Start, ":", Start + 1, ":", strand_number)) %>%
    dplyr::select(reg),
  file = file.path(data_dir, "signif_CpGs.regions"),
  delim = " ",
  col_names = FALSE
)

write_delim(
  diff_methyl_signif %>%
    mutate(Chr = substr(Chromosome, 4, nchar(Chromosome))) %>%
    mutate(reg = paste0(Chr, ":", Start, ":", Start + 1)) %>%
    dplyr::select(reg),
  file = file.path(data_dir, "signif_CpGs_noStrand.regions"),
  delim = " ",
  col_names = FALSE
)
```

------------------------------------------------------------------------

# Web services

We will use the results above to perform several searches online.

## ensembl BioMart - genes

The biomart on [ensembl](http://www.ensembl.org/biomart/martview/) is a
web service that allows to send batch queries to ensembl and get results
in a text format.

We will and use `signif_CpGs.regions` file to extract the genes that lay
near CpGs.

Here, you can choose how to do that:

-   [through the web browser](#biomart_web) or by [using `biomaRt`
    package](#biomart_pack).

### biomart on the web

-   Choose database: Ensembl Genes

-   Choose dataset: Human Genes

-   On the left (blue menu), choose `Filters`, expand `REGION` and
    select `Multiple regions`, then Browse to the file

-   Next, choose `Attributes` (below `Filters`) and expand `GENE`, then
    choose the following:

    -   Gene stable ID
    -   Protein stable ID
    -   Gene description
    -   Chromosome/scaffold name
    -   Gene start (bp)
    -   Gene end (bp)
    -   Gene name

-   Click on `Results` (above `Filters`) and choose to “Export all
    results to \| File \| TSV”, tick off “Unique results only” and click
    “Go”

-   Save the file to the disk, we will read it below

### using `biomaRt` package

``` r
# check all the servers by running biomaRt::listEnsembl()
# check all the datasets by running biomaRt::listDatasets(ensembl_mart)
ensembl_genes <- useEnsembl("genes",
  dataset = "hsapiens_gene_ensembl",
  GRCh = 37 # don't use this argument if you want the newest version
)

# check all the attributes by running biomaRt::listAttributes(ensembl_mart)
my_attribs <- c(
  "ensembl_gene_id",
  "ensembl_peptide_id",
  "description",
  "chromosome_name",
  "start_position",
  "end_position",
  "external_gene_name"
)
# check all the attributes by running biomaRt::listFilters(ensembl_mart)
my_filters <- c(
  "chromosome_name",
  "start",
  "end",
  "strand"
)

# we need to have an input data that is a data.frame with columns
#  as in 'my_filters'
data_in <- diff_methyl_signif %>%
    mutate(End = Start + 1) %>%
    dplyr::select(
      chromosome_name = Chromosome,
      start = Start,
      end = End,
      strand = Strand
    )
data_in
```

    ## # A tibble: 4,492 × 4
    ##    chromosome_name   start     end strand
    ##    <chr>             <dbl>   <dbl> <chr> 
    ##  1 chr1             568475  568476 -     
    ##  2 chr1             896310  896311 -     
    ##  3 chr1             995013  995014 -     
    ##  4 chr1            1020738 1020739 +     
    ##  5 chr1            1051364 1051365 +     
    ##  6 chr1            1334803 1334804 +     
    ##  7 chr1            1543828 1543829 -     
    ##  8 chr1            1710208 1710209 +     
    ##  9 chr1            1822399 1822400 +     
    ## 10 chr1            2126425 2126426 -     
    ## # … with 4,482 more rows

``` r
out_data_colnames <- c("Gene_ID", "Protein_ID", "Chr", "Start", "End",
                       "Gene_name", "Gene_descr")

# since this requires internet connection, we will save it in a file after
#   downloading the first time and read from a file next time;
#   if we want to update the information, we will need to delete this file
#   before running the code
genes_file_out <- file.path(data_dir, "mart_export.txt")

if(!file.exists(genes_file_out)){
# this might take some time...
  genes_near_signif_cpgs <- purrr::map(seq_len(nrow(data_in)), function(row){
    cur_cpg <- as.list(data_in[row, ])
    cur_out <- biomaRt::getBM(
      attributes = my_attribs,
      filters = my_filters,
      values = cur_cpg,
      mart = ensembl_genes
    )
    if(nrow(cur_out) == 0){
      return(NULL)
    }
    return(cur_out)
  })
  genes_near_signif_cpgs <- dplyr::bind_rows(genes_near_signif_cpgs)
  names(genes_near_signif_cpgs) <- out_data_colnames
  
  write_delim(
    genes_near_signif_cpgs,
    file = file.path(data_dir, "mart_export.txt"),
    delim = "\t"
  )
}
```

### using results from biomart

``` r
# after getting results from ensembl BioMart
genes_near_signif_cpgs <- read_delim(
  file = file.path(data_dir, "mart_export.txt"),
  delim = "\t",
  col_names = out_data_colnames,
  skip = 1
)
```

    ## Rows: 759 Columns: 7
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (4): Gene_ID, Protein_ID, Gene_name, Gene_descr
    ## dbl (3): Chr, Start, End
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
genes_near_signif_cpgs
```

    ## # A tibble: 759 × 7
    ##    Gene_ID         Protein_ID        Chr    Start      End Gene_name  Gene_descr
    ##    <chr>           <chr>           <dbl>    <dbl>    <dbl> <chr>      <chr>     
    ##  1 ENSG00000120949 ENSP00000263932     1 12063303 12144207 TNFRSF8    TNF recep…
    ##  2 ENSG00000120949 ENSP00000390650     1 12063303 12144207 TNFRSF8    TNF recep…
    ##  3 ENSG00000120949 ENSP00000421938     1 12063303 12144207 TNFRSF8    TNF recep…
    ##  4 ENSG00000120949 ENSP00000398337     1 12063303 12144207 TNFRSF8    TNF recep…
    ##  5 ENSG00000120949 <NA>                1 12063303 12144207 TNFRSF8    TNF recep…
    ##  6 ENSG00000261025 <NA>                1 24968423 24970865 AL445471.1 novel tra…
    ##  7 ENSG00000084070 ENSP00000389895     1 40344850 40423326 SMAP2      small Arf…
    ##  8 ENSG00000084070 ENSP00000361803     1 40344850 40423326 SMAP2      small Arf…
    ##  9 ENSG00000084070 ENSP00000361793     1 40344850 40423326 SMAP2      small Arf…
    ## 10 ENSG00000084070 ENSP00000479285     1 40344850 40423326 SMAP2      small Arf…
    ## # … with 749 more rows

``` r
# let's see how many results per chromosome
knitr::kable(
  genes_near_signif_cpgs %>%
    distinct(Gene_ID, Chr) %>%
    count(Chr)
)
```

| Chr |   n |
|----:|----:|
|   1 |  11 |
|   2 |  14 |
|   3 |   6 |
|   4 |   6 |
|   5 |   9 |
|   6 |   3 |
|   7 |   1 |
|   8 |   2 |
|   9 |   1 |
|  10 |   6 |
|  11 |   6 |
|  12 |   8 |
|  13 |   3 |
|  14 |   4 |
|  15 |   3 |
|  16 |   3 |
|  17 |  10 |
|  18 |   2 |
|  19 |   8 |
|  20 |   5 |
|  21 |   1 |
|  22 |   1 |

To use the results in the next step, we need only a list of the gene
names.

``` r
write_delim(
  genes_near_signif_cpgs %>%
    distinct(Gene_name),
  file = file.path("data", "genes_near_signif_cpgs.txt"),
  delim = " ",
  col_names = FALSE
)
```

## STRING - database of protein-protein interactions (PPI)

We can check whether our set of genes is enriched in something specific
or has other common features or functions by e.g., checking it in
[STRING](https://string-db.org/)

-   Click on SEARCH

-   Choose `Multiple proteins` in the list to the left

-   Click on the “Browse” button and navigate to
    `genes_near_signif_cpgs.txt` file

-   Enter `Homo sapiens` in the “Organism” field and click “Search”

-   You will be presented with possible matches to the identifiers that
    we provided - accept all and continue

-   Explore the network, e.g.,

    -   in the “Analysis” tab, check enrichment
    -   in the “Clusters” tab, play with k-means clustering
    -   what is the protein “in the middle”?

### ensembl BioMart - regulatory regions

ensembl browser has many databases - one contains regulatory features of
the genome (such as promoters and enhancers). Export the regulatory
regions that are in vicinity of our selected CpGs.

Here, you can choose how to do that:

-   [through the web browser](#biomart_regul_web) or by [using `biomaRt`
    package](#biomart_regul_pack).

## biomart on the web

-   Choose database: Ensembl Regulation

-   Choose dataset: Human Regulatory Features

-   On the left (blue menu), choose `Filters`, expand `REGION` and
    select `Multiple regions`, then Browse to the file
    `signif_CpGs_noStrand.regions`

-   Next, choose `Attributes` (below `Filters`) and expand `GENE`, then
    choose the following:

    -   Chromosome/scaffold name
    -   Start (bp)
    -   End (bp)
    -   Feature type
    -   Regulatory stable ID
    -   SO term accession

-   Click on `Results` (above `Filters`) and choose to “Export all
    results to \| File \| TSV”, tick off “Unique results only” and click
    “Go”

-   Save the file to the disk as `mart_export_regulatory_feat.txt`, we
    will read it below

### using `biomaRt` package

``` r
# check all the servers by running biomaRt::listEnsembl()
# check all the datasets by running biomaRt::listDatasets(ensembl_mart)
ensembl_reg <- useEnsembl("regulation",
  dataset = "hsapiens_regulatory_feature",
  GRCh = 37 # don't use this argument if you want the newest version
)

# check all the attributes by running biomaRt::listAttributes(ensembl_mart)
my_attribs <- c(
  "chromosome_name",
  "chromosome_start",
  "chromosome_end",
  "feature_type_name",
  "regulatory_stable_id",
  "so_accession"
)
# check all the attributes by running biomaRt::listFilters(ensembl_mart)
my_filters <- c(
  "chromosome_name",
  "start",
  "end"
)

data_in <- diff_methyl_signif %>%
    mutate(End = Start + 1) %>%
    dplyr::select(
      chromosome_name = Chromosome,
      start = Start,
      end = End
    )
data_in
```

    ## # A tibble: 4,492 × 3
    ##    chromosome_name   start     end
    ##    <chr>             <dbl>   <dbl>
    ##  1 chr1             568475  568476
    ##  2 chr1             896310  896311
    ##  3 chr1             995013  995014
    ##  4 chr1            1020738 1020739
    ##  5 chr1            1051364 1051365
    ##  6 chr1            1334803 1334804
    ##  7 chr1            1543828 1543829
    ##  8 chr1            1710208 1710209
    ##  9 chr1            1822399 1822400
    ## 10 chr1            2126425 2126426
    ## # … with 4,482 more rows

``` r
reg_data_colnames <- c("chr", "start", "end", "feat_type", "reg_ID", "SO_term")

# since this requires internet connection, we will save it in a file after
#   downloading the first time and read from a file next time;
#   if we want to update the information, we will need to delete this file
#   before running the code
regul_file_out <- file.path(data_dir, "mart_export_regulatory_feat.txt")

if(!file.exists(regul_file_out)){
# this might take some time...
  regulatory_near_signif_cpgs <- purrr::map(seq_len(nrow(data_in)), function(row){
    cur_cpg <- as.list(data_in[row, ])
    cur_out <- biomaRt::getBM(
      attributes = my_attribs,
      filters = my_filters,
      values = cur_cpg,
      mart = ensembl_reg
    )
    if(nrow(cur_out) == 0){
      return(NULL)
    }
    return(cur_out)
  })
  regulatory_near_signif_cpgs <- dplyr::bind_rows(regulatory_near_signif_cpgs)
  names(regulatory_near_signif_cpgs) <- reg_data_colnames
  write_delim(
    regulatory_near_signif_cpgs,
    file = regul_file_out,
    delim = "\t"
  )
}
```

### using results from biomart

``` r
# after checking regulatory regions
regulatory_near_signif_cpgs <- read_delim(
  file = regul_file_out,
  delim = "\t",
  skip = 1,
  col_names = reg_data_colnames
)
```

    ## Rows: 66 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): feat_type, reg_ID, SO_term
    ## dbl (3): chr, start, end
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
regulatory_near_signif_cpgs
```

    ## # A tibble: 66 × 6
    ##      chr     start       end feat_type                reg_ID          SO_term   
    ##    <dbl>     <dbl>     <dbl> <chr>                    <chr>           <chr>     
    ##  1     1   2158602   2161400 Promoter Flanking Region ENSR00000344819 SO:0001952
    ##  2     1 232764855 232769000 Promoter Flanking Region ENSR00000390326 SO:0001952
    ##  3     1  24969401  24969600 CTCF Binding Site        ENSR00000351034 SO:0001974
    ##  4     1  12079201  12080600 CTCF Binding Site        ENSR00000250247 SO:0001974
    ##  5     1  12078600  12079601 Promoter                 ENSR00000919835 SO:0000167
    ##  6    10  45450002  45455400 Promoter Flanking Region ENSR00000402406 SO:0001952
    ##  7    11  96122001  96126199 Promoter Flanking Region ENSR00000439752 SO:0001952
    ##  8    11  67164802  67171799 Promoter Flanking Region ENSR00000433249 SO:0001952
    ##  9    11   1410402   1411800 Promoter Flanking Region ENSR00000952322 SO:0001952
    ## 10    11    840800    846401 Promoter                 ENSR00000035700 SO:0000167
    ## # … with 56 more rows

``` r
regulatory_near_signif_cpgs %>%
  count(SO_term) %>%
  left_join(regulatory_near_signif_cpgs %>%
              distinct(SO_term, feat_type))
```

    ## Joining, by = "SO_term"

    ## # A tibble: 6 × 3
    ##   SO_term        n feat_type               
    ##   <chr>      <int> <chr>                   
    ## 1 SO:0000165     2 Enhancer                
    ## 2 SO:0000167    23 Promoter                
    ## 3 SO:0000235     2 TF binding site         
    ## 4 SO:0001747     2 Open chromatin          
    ## 5 SO:0001952    29 Promoter Flanking Region
    ## 6 SO:0001974     8 CTCF Binding Site

``` r
regulatory_near_signif_cpgs %>%
  filter(feat_type == "Enhancer")
```

    ## # A tibble: 2 × 6
    ##     chr     start       end feat_type reg_ID          SO_term   
    ##   <dbl>     <dbl>     <dbl> <chr>     <chr>           <chr>     
    ## 1    12  98908801  98910200 Enhancer  ENSR00000467225 SO:0000165
    ## 2     3 197475801 197477000 Enhancer  ENSR00000714253 SO:0000165

How many different regulatory features were found?

## GenEnhancer by GeneCards

Choose one enhancer (by `reg_ID` - this is the ensembl ID of the region)
and enter it in the “Keywords” search field at
<https://www.genecards.org/>. *(NOTE: this is not the main search field,
which only searches for gene names)*

What genes were found? Check their Cards and see where are they located
on the genome in relation to the enhancer you’ve chosen.
