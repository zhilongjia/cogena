#! /usr/bin/Rscript

gmt2anno <- function(gmt_fns) {
	require(GSEABase)

    for (i in gmt_fns) {
        assign(i, geneIds(getGmt(i)))
        #save(i, file=paste0(i, ".RData"))
    }
}
gmt_fns <- c("c2.cgp.v4.0.symbols.gmt", "c2.cp.kegg.v4.0.symbols.gmt", "c2.cp.reactome.v4.0.symbols.gmt", 
	"c2.cp.v4.0.symbols.gmt", "c5.bp.v4.0.symbols.gmt", "c6.all.v4.0.symbols.gmt")
gmt2anno(gmt_fns)
save.image("annotation.RData")