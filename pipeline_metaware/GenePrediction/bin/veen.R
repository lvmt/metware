library(optparse)
library(ggplot2)
library(ggsci)
library(sf)
library(ggVennDiagram)

#自定义颜色；
color1 <- alpha("#f8766d",0.9)
color2 <- alpha("#FF99CC",0.7)
color3 <- alpha("#c77cff",0.5)
color4 <- alpha("#99CC00",0.5)

print('aaaa')
option_list <- list(
    make_option(c('-g', '--gene_table'), help='基因文件'),
    make_option(c('-c', '--groups'), help='分组信息a_b_c,b_c_d等')
)
opt = parse_args(OptionParser(option_list = option_list, usage = "This Script is a test for arguments!"))

print('bbb')
print(opt)

draw_veen <- function(df, single_group) {
    png("aaa.png")
    gene_names = rownames(df)
    stat_list = list()
    samples = unlist(strsplit(single_group, '_'))
    for(sample in samples) {
        stat_list[[sample]] = gene_names[as.numeric(df[, sample]) > 0]
    }

    ggVennDiagram(stat_list, label_alpha=0) + scale_fill_gradient(low="white",high =color4 ,guide="none")
    dev.off()
}


############################################### 
# df = read.csv(opt$gene_table,header = TRUE, row=1, sep='\t')
# groups_list = unlist(strsplit(opt$groups, ','))
# gene_names = rownames((df))

# for(single_group in groups_list) {
#     print(single_group)
#     draw_veen(df, single_group)
# }




