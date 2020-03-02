library(ggplot2)
library(Biostrings)
library(hues)

# Setting up file names

## I SUSPECT YOUR ERRORS WILL BE IN THIS SECTION
# folder where your files are  
base = "/home/moritz/people/0023_anoxicencyclo/course"
sample_name = "Loc090907-8-6m"
# You will need to run jgi_summarize_bam_contig_depthsfrom "module load MetaBat" with the --outputDepth option  on your sorted BAM  to get make that file 
mapping_file = file.path(base, "bins", sample_name, "final.contigs.fa.depth.txt")
# your assembly file
assembly_name = file.path(base, "assemblies", sample_name, "final.contigs.fa")

binning_file = NA #file.path(base, "bins", sample_name, "final.contigs.fa.depth.txt") # we will not use that in the first part of the tutorial ##FIX THIS WHEN WE DO BINNING
out_table_file = file.path(base, "script_out", "assembly_table.csv")
out_plot_file = file.path(base, "script_out", "assembly_plot.pdf")

# loading and cleaning up the coverage table
coverage_data = read.table(mapping_file, sep="\t", h=T, as.is = TRUE)
coverage_data$contig_name = sapply(strsplit(coverage_data$contigName, " "), "[",1)
row.names(coverage_data) = coverage_data$contig_name
coverage_data = coverage_data[,c("contigLen","totalAvgDepth")]

# loading assembly
assembly = readDNAStringSet(assembly_name)

# computing GC-content of assembly
gc_content = as.vector(letterFrequency(assembly, "GC")/width(assembly))
names(gc_content) = sapply(strsplit(names(assembly), " "), "[",1)
coverage_data$gc_content = gc_content[row.names(coverage_data)]
coverage_data$bin = NA

# if you have a some bins, parse them here to find out which contigs are in which bin # we will not use that in the first part of the tutorial
if(!is.na(binning_file))
{
  prokkas_folder = file.path(base, "all_bins")
  bin_folders = grep( sample_name, list.files(prokkas_folder), value = TRUE)
  all_bins = sapply(bin_folders, function(x) readDNAStringSet(file.path(base, "all_bins", x, paste0(x, ".fna"))))
  contig2bin = unlist(lapply(names(all_bins), function(x) {
    contigs = names(all_bins[[x]]);
    mapi = rep(x,length(contigs))
    names(mapi) = contigs
    mapi
  }))
  coverage_data$bin = contig2bin[row.names(coverage_data)]
}

# fix colours for the plot
coverage_data$bin[is.na(coverage_data$bin)] = "unbinned"
colmap = as.vector(iwanthue(length(levels(factor(coverage_data$bin)))))
colmap[length(colmap)] = "#999999"

# make a pretty plot
p1 = ggplot(coverage_data[coverage_data$contigLen > 2500,], aes(x=gc_content, y=totalAvgDepth, size=contigLen, col=bin))+geom_point()+scale_y_log10()+scale_x_log10()+scale_color_manual(values = colmap)+theme_minimal()

# prep output data
out_data = list(
  total_length = sum(width(assembly)),
  mean_length = mean(width(assembly)),
  median_length = median(width(assembly)),
  gc_content = sum(letterFrequency(assembly, "GC"))/sum(width(assembly)),
  nb_contigs = length(assembly),
  N50 = sum(cumsum(sort(width(assembly))) > (sum(width(assembly))/2)),
  N50_len = sort(width(assembly))[cumsum(sort(width(assembly))) > (sum(width(assembly))/2)][1]
)
out_data = data.frame(out_data)
row.names(out_data) = c(sample_name)

#write output
write.table(out_data, file = out_table_file)
ggsave(file = out_plot_file, plot = p1)
